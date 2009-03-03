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

#include "cell.h"
#include "landscape.h"

#include <iostream>
#include <iomanip>

using namespace gingko;

///////////////////////////////////////////////////////////////////////////////	
// Cell

std::vector<const Organism*> Cell::breeding_female_ptrs; // scratch space for breeding
std::vector<const Organism*> Cell::breeding_male_ptrs;   // scratch space for breeding
OrganismVector Cell::previous_gen;                       // scratch space to hold previous generation during reproduction

// --- lifecycle and assignment ---

Cell::Cell(CellIndexType index, 
           unsigned num_fitness_factors,
           Landscape& landscape, 
           const SpeciesPointerVector& species, 
           RandomNumberGenerator& rng)     
    : index_(index),
      num_fitness_factors_(num_fitness_factors),
      landscape_(landscape),
      species_(species),
      rng_(rng) {
    this->carrying_capacity_ = 0;
    memset(this->environment_, 0, this->num_fitness_factors_*sizeof(FitnessFactorType));    
}

// --- primary biogeographical and evolutionary processes ---

void Cell::reproduction() {
	
	Cell::previous_gen.clear();
	Cell::previous_gen.swap(this->organisms_);
		
    for (SpeciesPointerVector::const_iterator sp = this->species_.begin(); sp != this->species_.end(); ++sp) {
		Cell::breeding_female_ptrs.clear();
		Cell::breeding_male_ptrs.clear();

        // species-level reproduction rate for now: later this will be at the 
        // organism level and subject to evolution
        unsigned num_offspring = (*sp)->get_mean_reproductive_rate();
        
        this->extract_breeding_groups((*sp)->get_index(), Cell::breeding_female_ptrs, Cell::breeding_male_ptrs);
        if ( (Cell::breeding_female_ptrs.size() > 0) and (Cell::breeding_male_ptrs.size() > 0)) {
            for (std::vector<const Organism*>::iterator fptr = Cell::breeding_female_ptrs.begin();
                    fptr != Cell::breeding_female_ptrs.end();
                    ++fptr) {
                for (unsigned n = 0; n <= num_offspring; ++n) {                    
                    const Organism* male = this->rng_.select(breeding_male_ptrs);
                    const Organism* female = *fptr;
                    this->organisms_.push_back((*sp)->new_organism(*female, *male));
                } // for each offspring     
            } // for each female
        } // if females > 0 and males > 0
    }  // for each species  
}

void Cell::migration() {

    for (OrganismVector::iterator og = this->organisms_.begin(); og != this->organisms_.end(); ++og) {
    
        assert(og->species_index() < this->species_.size());          
        assert(!og->is_expired());        
        
        Species& sp = *this->species_[og->species_index()];
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
        assert(og->species_index() < this->species_.size());          
        assert(!og->is_expired());                
        Species& sp = *this->species_[og->species_index()];
        float fitness = sp.calc_fitness(*og, this->environment_);       
        og->set_fitness(fitness);
        if (this->rng_.uniform_real() > fitness) {
            og->set_expired(true);
        }
    }
    this->purge_expired_organisms();
}

void Cell::competition() {
    // NOTE: ASSUMES THAT FITNESS HAS BEEN CALCULATED FOR THE ORGANISM IN THIS CELL!
//     std::cout << this->organisms_.size() << ", " << this->carrying_capacity_;
    if (this->organisms_.size() > this->carrying_capacity_) {
        // defaults to using Organism::operator<(), which also checks that 
        // fitness has been set (i.e., >= 0)
        std::sort(this->organisms_.begin(), this->organisms_.end());
        this->organisms_.erase(this->organisms_.begin()+this->carrying_capacity_, 
            this->organisms_.end());
        assert(this->organisms_.size() == this->carrying_capacity_);
    }
    assert(this->organisms_.size() <= this->carrying_capacity_);
//     std::cout << " / " << this->organisms_.size() << ", " << this->carrying_capacity_ << std::endl;
}


// --- supporting biogeographical and evolutionary processes ---

//! Extracts pointers to male and female organisms of a particular species
void Cell::extract_breeding_groups(unsigned species_index, 
                                        std::vector<const Organism*>& female_ptrs,
                                        std::vector<const Organism*>& male_ptrs) const {
    assert(species_index < this->species_.size());    
    for (OrganismVector::const_iterator og = this->organisms_.begin(); og != this->organisms_.end(); ++og) {
        if (og->species_index() == species_index) {
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
