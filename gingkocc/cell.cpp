///////////////////////////////////////////////////////////////////////////////
//
// GINGKO Biogeographical Evolution Simulator.
//
// Copyright 2009 Jeet Sukumaran and Mark T. Holder.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License along
// with this program. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#include "cell.hpp"
#include "landscape.hpp"
#include "randgen.hpp"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <set>

using namespace gingko;

///////////////////////////////////////////////////////////////////////////////	
// Cell

std::vector<const Organism*> Cell::breeding_female_ptrs; // scratch space for breeding
std::vector<const Organism*> Cell::breeding_male_ptrs;   // scratch space for breeding
OrganismVector Cell::previous_gen;                       // scratch space to hold previous generation during reproduction

// --- lifecycle and assignment ---

Cell::Cell(CellIndexType index,
           CellIndexType x,
           CellIndexType y,
           unsigned num_fitness_factors,
           Landscape& landscape, 
           const SpeciesByLabel& species, 
           RandomNumberGenerator& rng)     
        : index_(index),
          x_(x),
          y_(y),
          carrying_capacity_(0),
          num_fitness_factors_(num_fitness_factors),
          landscape_(landscape),
          species_(species),
          rng_(rng) {   
    memset(this->environment_, 0, 
    MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));         
}

// --- basic biotics ---

void Cell::generate_new_organisms(Species * sp, CellIndexType num) {
    this->organisms_.reserve(this->organisms_.size() + num);
    for ( ; num > 0; --num) {
        this->organisms_.push_back(sp->new_organism());
    }
}

// --- primary biogeographical and evolutionary processes ---

void Cell::reproduction() {
	
	Cell::previous_gen.clear();
	Cell::previous_gen.swap(this->organisms_);
		
    for (SpeciesByLabel::const_iterator spi = this->species_.begin(); spi != this->species_.end(); ++spi) {
        Species * sp = spi->second;
		Cell::breeding_female_ptrs.clear();
		Cell::breeding_male_ptrs.clear();

        // species-level reproduction rate for now: later this will be at the 
        // organism level and subject to evolution
        unsigned num_offspring = sp->get_mean_reproductive_rate();
        
        this->extract_breeding_groups(sp,
            Cell::previous_gen,
            Cell::breeding_female_ptrs, 
            Cell::breeding_male_ptrs);
        if ( (Cell::breeding_female_ptrs.size() > 0) and (Cell::breeding_male_ptrs.size() > 0)) {
            for (std::vector<const Organism*>::iterator fptr = Cell::breeding_female_ptrs.begin();
                    fptr != Cell::breeding_female_ptrs.end();
                    ++fptr) {
                for (unsigned n = 0; n <= num_offspring; ++n) {                    
                    const Organism* male = this->rng_.select(breeding_male_ptrs);
                    const Organism* female = *fptr;
                    this->organisms_.push_back(sp->new_organism(*female, *male));
                } // for each offspring     
            } // for each female
        } // if females > 0 and males > 0
    }  // for each species

    // post reproduction shuffle, to destroy correlation between adjacent
    // organisms
    RandomPointer rp(this->rng_);
    std::random_shuffle(this->organisms_.begin(), this->organisms_.end(), rp);    
    
}

void Cell::migration() {

    for (OrganismVector::iterator og = this->organisms_.begin(); og != this->organisms_.end(); ++og) {             
        assert(!og->is_expired());                
        Species& sp = og->species();
        int movement = sp.get_movement_capacity();
        CellIndexType curr_idx = this->index_;
                        
        while (movement > 0) {
            CellIndexType dest_idx = this->landscape_.random_neighbor(curr_idx);
            movement -= sp.movement_cost(dest_idx);
            if (movement >= 0) {
                curr_idx = dest_idx;
            }
        } 
        
        if (curr_idx != this->index_) {
            this->landscape_.add_migrant(*og, curr_idx);
            og->set_expired(true);            
        }
    }
    this->purge_expired_organisms();
}

void Cell::survival() {
    for (OrganismVector::iterator og = this->organisms_.begin(); og != this->organisms_.end(); ++og) {        
        assert(!og->is_expired());                
        Species& sp = og->species();
        float fitness = sp.calc_fitness(*og, this->environment_);       
        og->set_fitness(fitness);
        if (this->rng_.uniform_01() > fitness) {
            og->set_expired(true);
        }
    }
    this->purge_expired_organisms();
}

void Cell::competition() {
    // NOTE: ASSUMES THAT FITNESS HAS BEEN CALCULATED FOR THE ORGANISM IN THIS CELL!
    // This would have been done during the survival phase.
    if (this->organisms_.size() > this->carrying_capacity_) {
    
        // shuffle vector, so that if no selection is operating, random 
        // selection will determine competition winners
//         RandomPointer rp(this->rng_);
//         std::random_shuffle(this->organisms_.begin(), this->organisms_.end(), rp);        
        
        // build set of organisms sorted by fitness
        std::multiset<Organism *, CompareOrganismFitness> optrs;
        for (OrganismVector::iterator oi = this->organisms_.begin();
                oi != this->organisms_.end();
                ++oi) {
            optrs.insert( &(*oi) );
        }
        assert(optrs.size() == this->organisms_.size());
        
        // find winners
        OrganismVector winners;
        winners.reserve(this->carrying_capacity_);
        unsigned long count = 0;
        for (std::multiset<Organism *, CompareOrganismFitness>::iterator opi = optrs.begin();
                count < this->carrying_capacity_;
                ++opi, ++count) {
            winners.push_back(**opi);                
        }
        
        this->organisms_.swap(winners);
    }
    assert(this->organisms_.size() <= this->carrying_capacity_);
//     std::cout << " / " << this->organisms_.size() << ", " << this->carrying_capacity_ << std::endl;
}

// --- for trees etc ---

void Cell::sample_organisms(Species * sp_ptr,
    std::vector<const Organism *>& samples, unsigned long num_organisms) {
    std::vector<const Organism *> species_organisms;
    species_organisms.reserve(this->organisms_.size());
    for (OrganismVector::const_iterator og = this->organisms_.begin(); og != this->organisms_.end(); ++og) {
        if (&og->species() == sp_ptr) {
            sp_ptr->get_organism_label(*og, this->x_, this->y_);
            species_organisms.push_back(&(*og));
        }
    }
    if (num_organisms == 0 or num_organisms >= species_organisms.size()) {
        samples.reserve(samples.size() + species_organisms.size());
        std::copy(species_organisms.begin(), species_organisms.end(), std::back_inserter(samples));
    } else {
        samples.reserve(samples.size() + num_organisms);
        RandomPointer rp(this->rng_);
        std::random_shuffle(species_organisms.begin(), species_organisms.end(), rp);
        std::copy(species_organisms.begin(), species_organisms.begin() + num_organisms, std::back_inserter(samples)); 
    }
}

unsigned long Cell::num_organisms(Species * sp_ptr) const {
    unsigned long count = 0;
    for (OrganismVector::const_iterator og = this->organisms_.begin(); og != this->organisms_.end(); ++og) {
        if (&og->species() == sp_ptr) {
            ++count;
        }
    }
    return count;
}   

// --- supporting methods ---

//! Extracts pointers to male and female organisms of a particular species from 
//! a vector of organisms passed to it.
void Cell::extract_breeding_groups(Species * sp_ptr,
        const OrganismVector& organisms,
        std::vector<const Organism*>& female_ptrs,
        std::vector<const Organism*>& male_ptrs) const { 
    for (OrganismVector::const_iterator og = organisms.begin(); og != organisms.end(); ++og) {
        if (&og->species() == sp_ptr) {
            if (og->is_female()) {
                female_ptrs.push_back(&(*og));
            } else {
                male_ptrs.push_back(&(*og));
            }
        }
    }
}                                         

// Cell
///////////////////////////////////////////////////////////////////////////////	
