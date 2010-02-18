///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Biogeographical Evolution Simulator.
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

using namespace ginkgo;

///////////////////////////////////////////////////////////////////////////////
// Cell

OrganismVector Cell::next_gen_females;
OrganismVector Cell::next_gen_males;

// --- lifecycle and assignment ---

Cell::Cell(CellIndexType index,
           CellIndexType x,
           CellIndexType y,
           unsigned num_fitness_traits,
           Landscape& landscape,
           const SpeciesByLabel& species,
           RandomNumberGenerator& rng)
        : index_(index),
          x_(x),
          y_(y),
          carrying_capacity_(0),
          num_fitness_traits_(num_fitness_traits),
          landscape_(landscape),
          species_(species),
          populations_(species, rng),
          rng_(rng) {
    memset(this->fitness_trait_optimum_, 0,
    MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));
}

// --- basic biotics ---

void Cell::generate_new_organisms(Species * sp, PopulationCountType num) {
    BreedingPopulation& pop = this->populations_[sp];
    for ( ; num > 0; --num) {
        pop.add(sp->new_organism(this->index_));
    }
}

void Cell::generate_new_population(Species * sp,
        PopulationCountType final_pop_size,
        PopulationCountType ancestral_pop_size,
        GenerationCountType ancestral_generations) {
    if (ancestral_pop_size == 0) {
        ancestral_pop_size = final_pop_size;
    }
    if (ancestral_generations == 0) {
        ancestral_generations = ancestral_pop_size * 10;
    }
    Cell temp_cell(this->index_, this->x_, this->y_, this->num_fitness_traits_, this->landscape_, this->species_, this->rng_);
    temp_cell.generate_new_organisms(sp, ancestral_pop_size);
    for (GenerationCountType g = 0; g != ancestral_generations; ++g) {
        temp_cell.reproduction(false); // reproduce without evolving fitness component traits
        temp_cell.populations_.retain(ancestral_pop_size);
        // std::random_shuffle(temp_cell.organisms_.begin(), temp_cell.organisms_.end(), rp);
    }
    temp_cell.populations_.retain(final_pop_size);
    this->populations_ = temp_cell.populations_;
}

// --- primary biogeographical and evolutionary processes ---

void Cell::reproduction(bool evolve_fitness_components) {

    for (SpeciesByLabel::const_iterator spi = this->species_.begin(); spi != this->species_.end(); ++spi) {
        Species * sp = spi->second;
        BreedingPopulation pop = this->populations_[sp];
        OrganismVector& females = pop.females();
        OrganismVector& males = pop.males();
        if ( (females.size() > 0) and (males.size() > 0)) {
            BreedingPopulation next_gen;
            unsigned num_offspring = sp->get_mean_reproductive_rate();
            for (OrganismVector::iterator fi = females.begin();
                    fi != females.end();
                    ++fi) {
                for (unsigned n = 0; n <= num_offspring; ++n) {
                    const Organism* male = this->rng_.select_ptr(males);
                    next_gen.add(sp->new_organism(*fi, *male, this->index_, evolve_fitness_components));
                } // for each offspring
            } // for each female
            this->populations_[sp].swap(next_gen);

            // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // TODO: check if this is neccessary
            this->populations_[sp].shuffle(this->rng_);
            // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        } else {
            this->populations_[sp].clear();
        } // if females == 0 or males == 0
    }  // for each species
}

void Cell::migration() {
    for (SpeciesByLabel::const_iterator spi = this->species_.begin();
            spi != this->species_.end();
            ++spi) {
        Species * sp = spi->second;
        BreedingPopulation& pop = this->populations_[sp];
    }
//
//    for (OrganismVector::iterator og = this->organisms_.begin(); og != this->organisms_.end(); ++og) {
//        assert(!og->is_expired());
//        Species& sp = og->species();
//
//        if (this->rng_.uniform_01() <= sp.movement_probability(this->index_)) {
//            MovementCountType movement = sp.get_movement_capacity();
//            CellIndexType curr_idx = this->index_;
//            MovementCountType movement_cost = 0;
//            while (movement > 0) {
//                CellIndexType dest_idx = this->landscape_.random_neighbor(curr_idx);
//                assert(dest_idx < this->landscape_.size());
//                movement_cost = sp.movement_cost(dest_idx);
//                if ( movement >= movement_cost) {
//                    movement -= movement_cost;
//                    curr_idx = dest_idx;
//                } else {
//                    movement -= sp.movement_cost(curr_idx);
//                }
//            }
//
//            if (curr_idx != this->index_) {
//                this->landscape_.add_migrant(*og, curr_idx);
//                og->set_expired();
//            }
//        }
//    }
//    this->purge_expired_organisms();
}

void Cell::survival() {
//    for (OrganismVector::iterator og = this->organisms_.begin(); og != this->organisms_.end(); ++og) {
//        assert(!og->is_expired());
//        Species& sp = og->species();
//        float fitness = sp.calc_fitness(*og, this->fitness_trait_optimum_);
//        og->set_fitness(fitness);
//        if (this->rng_.uniform_01() > fitness) {
//            og->set_expired();
//        }
//    }
//    this->purge_expired_organisms();
}

void Cell::competition() {
//    // NOTE: ASSUMES THAT FITNESS HAS BEEN CALCULATED FOR THE ORGANISM IN THIS CELL!
//    // This would have been done during the survival phase.
//    if (this->organisms_.size() > this->carrying_capacity_) {
//
//        // build set of organisms sorted by fitness
//        std::multiset<Organism *, CompareOrganismFitness> optrs;
//        for (OrganismVector::iterator oi = this->organisms_.begin();
//                oi != this->organisms_.end();
//                ++oi) {
//            optrs.insert( &(*oi) );
//        }
//        assert(optrs.size() == this->organisms_.size());
//
//        // find winners
//        OrganismVector winners;
//        winners.reserve(this->carrying_capacity_);
//        unsigned long count = 0;
//        for (std::multiset<Organism *, CompareOrganismFitness>::iterator opi = optrs.begin();
//                count < this->carrying_capacity_;
//                ++opi, ++count) {
//            winners.push_back(**opi);
//        }
//
//        this->organisms_.swap(winners);
//    }
//    assert(this->organisms_.size() <= this->carrying_capacity_);

}

// --- for trees etc ---

void Cell::sample_organisms(Species * sp_ptr,
    std::vector<const Organism *>& samples, PopulationCountType num_organisms) {
//    std::vector<const Organism *> available_organisms;
//    available_organisms.reserve(this->organisms_.size());
//    for (OrganismVector::const_iterator og = this->organisms_.begin(); og != this->organisms_.end(); ++og) {
//        if (sp_ptr == NULL) {
//            (*og).species().set_organism_label(*og, this->index_, this->x_, this->y_);
//            available_organisms.push_back(&(*og));
//        } else if (&og->species() == sp_ptr) {
//            sp_ptr->set_organism_label(*og, this->index_, this->x_, this->y_);
//            available_organisms.push_back(&(*og));
//        }
//    }
//    if (num_organisms == 0 or num_organisms >= available_organisms.size()) {
//        samples.reserve(samples.size() + available_organisms.size());
//        std::copy(available_organisms.begin(), available_organisms.end(), std::back_inserter(samples));
//    } else {
//        samples.reserve(samples.size() + num_organisms);
//        RandomPointer rp(this->rng_);
//        std::random_shuffle(available_organisms.begin(), available_organisms.end(), rp);
//        std::copy(available_organisms.begin(), available_organisms.begin() + num_organisms, std::back_inserter(samples));
//    }
}

void Cell::num_organisms(Species * species_ptr, PopulationCountType& num_females, PopulationCountType& num_males) const {
/// TODO !!! ///
}

PopulationCountType Cell::num_organisms(Species * sp_ptr) const {
    PopulationCountType  count = 0;
//    for (OrganismVector::const_iterator og = this->organisms_.begin(); og != this->organisms_.end(); ++og) {
//        if (&og->species() == sp_ptr) {
//            ++count;
//        }
//    }
    return count;
}

// Cell
///////////////////////////////////////////////////////////////////////////////
