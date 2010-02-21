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

OrganismMemoryManager& ORGANISM_MM = OrganismMemoryManager::get_instance();

///////////////////////////////////////////////////////////////////////////////
// Cell

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
          rng_(rng),
          organism_memory_manager_(ORGANISM_MM) {
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
        BreedingPopulation& pop = this->populations_[sp];

        // TODO: instead of handles to the male and female vectors,
        //       BreedingPopulation should return iterators

        OrganismPointers& female_ptrs = pop.females();
        OrganismPointers& male_ptrs = pop.males();
        if ( (female_ptrs.size() > 0) and (male_ptrs.size() > 0)) {
            BreedingPopulation next_gen;
            unsigned num_offspring = sp->get_mean_reproductive_rate();
            for (OrganismPointers::iterator fi = female_ptrs.begin();
                    fi != female_ptrs.end();
                    ++fi) {
                const Organism* female_ptr = *fi;
                for (unsigned n = 0; n <= num_offspring; ++n) {
                    const Organism* male_ptr = *(this->rng_.select_ptr(male_ptrs));
                    next_gen.add(sp->new_organism(female_ptr, male_ptr, this->index_, evolve_fitness_components));
                } // for each offspring
            } // for each female
            this->populations_[sp].clear();
            this->populations_[sp].swap(next_gen);
            // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // TODO: check if this is neccessary
            //this->populations_[sp].shuffle(this->rng_);
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
        BreedingPopulation& current_pop = this->populations_[sp];
        BreedingPopulation remaining_pop;
        for (BreedingPopulation::iterator oi = current_pop.begin(); oi != current_pop.end(); ++oi) {
            Organism * og_ptr = *oi;
            if (this->rng_.uniform_01() <= sp->movement_probability(this->index_)) {
                MovementCountType movement = sp->get_movement_capacity();
                CellIndexType curr_idx = this->index_;
                MovementCountType movement_cost = 0;
                while (movement > 0) {
                    CellIndexType dest_idx = this->landscape_.random_neighbor(curr_idx);
                    assert(dest_idx < this->landscape_.size());
                    movement_cost = sp->movement_cost(dest_idx);
                    if ( movement >= movement_cost) {
                        movement -= movement_cost;
                        curr_idx = dest_idx;
                    } else {
                        movement -= sp->movement_cost(curr_idx);
                    }
                }
                if (curr_idx != this->index_) {
                    this->landscape_.add_migrant(og_ptr, curr_idx);
                } else {
                    remaining_pop.add(og_ptr);
                }
            }
        }
        this->populations_[sp] = remaining_pop;
    }
}

PopulationCountType Cell::survival() {
    PopulationCountType num_survivors = 0;
    for (SpeciesByLabel::const_iterator spi = this->species_.begin();
            spi != this->species_.end();
            ++spi) {
        Species * sp = spi->second;
        BreedingPopulation& pop = this->populations_[sp];
        for (BreedingPopulation::iterator oi = pop.begin(); oi != pop.end(); ++oi) {
            Organism * og_ptr = *oi;
            float fitness = sp->calc_fitness(og_ptr, this->fitness_trait_optimum_);
            og_ptr->set_fitness(fitness);
            if (this->rng_.uniform_01() > fitness) {
                og_ptr->set_expired();
            } else {
                ++num_survivors;
            }
        }
    }
    this->purge_expired_organisms();
    return num_survivors;
}

void Cell::competition() {
    // NOTE: ASSUMES THAT FITNESS HAS BEEN CALCULATED FOR THE ORGANISM IN THIS CELL!
    // This would have been done during the survival phase.
    if (this->populations_.size() > this->carrying_capacity_) {

        CompareOrganismFitnessFuncPtrType comp_fitness_fptr = &compare_organism_fitness;

        // build set of organisms sorted by fitness
        std::multiset<Organism *, CompareOrganismFitnessFuncPtrType> organism_fitness_map(comp_fitness_fptr);
        for (SpeciesByLabel::const_iterator spi = this->species_.begin();
                spi != this->species_.end();
                ++spi) {
            Species * sp = spi->second;
            BreedingPopulation& pop = this->populations_[sp];
            for (BreedingPopulation::iterator oi = pop.begin(); oi != pop.end(); ++oi) {
                organism_fitness_map.insert(*oi);

                // we set the organism's expired flag here, and unset it if
                // it is in the top K organisms later; we do it this way
                // b/c of the failure in the original approach noted below
                (*oi)->set_expired();

            }
        }
        assert(organism_fitness_map.size() == this->populations_.size());

        // alternate approach: "rescue" the top K individuals from expiration
        unsigned long count = 0;
        for (std::multiset<Organism *, CompareOrganismFitnessFuncPtrType>::iterator opi = organism_fitness_map.begin();
                count < this->carrying_capacity_ && opi != organism_fitness_map.end();
                ++opi, ++count) {
            (*opi)->set_unexpired();
        }
        this->purge_expired_organisms();

/******************************************************************************
Original Approach
------------------
Following does not work if many organisms have the same fitness. Specifically,
the "opi != organism_fitness_map.end()" seems to evaluate to True even if it
has not actually reached the end.
******************************************************************************/

//        std::multiset<Organism *, CompareOrganismFitnessFuncPtrType>::iterator opi = organism_fitness_map.begin();
//        PopulationCountType retained = 0;
//
//        // --FOR DEBUGGING-- //
//        PopulationCountType original_size = this->populations_.size();
//        PopulationCountType removed = 0;
//        // --FOR DEBUGGING-- //
//
//        // advance past the top K individuals, where K == carrying capacity
//        while ( (retained < (this->carrying_capacity_-1) ) && opi != organism_fitness_map.end() ) {
//            ++retained;
//            ++opi;
//        }
//
//        // mark remaining individuals for expiration
//        while (opi != organism_fitness_map.end())  {
//            (*opi)->set_expired();
//            ++opi;
//
//            // --FOR DEBUGGING-- //
//            ++removed;
//            // --FOR DEBUGGING-- //
//
//        }
//        assert(retained+removed == original_size);
//
//        // expire
//        this->purge_expired_organisms();
//
//        // --FOR DEBUGGING-- //
//        if (!(this->populations_.size() <= this->carrying_capacity_)) {
//            std::cout << " Original: " << original_size << std::endl;
//            std::cout << " Capacity: " << this->carrying_capacity_ << std::endl;
//            std::cout << " Retained: " << retained << std::endl;
//            std::cout << "  Removed: " << removed << std::endl;
//            std::cout << "Remaining: " << this->populations_.size() << std::endl;
//            OrganismPointers orgs = this->populations_.get_organism_ptrs();
//            unsigned long expired = 0;
//            for (OrganismPointers::iterator oi = orgs.begin(); oi != orgs.end(); ++oi) {
//                if ((*oi)->is_expired()) {
//                    ++expired;
//                }
//            }
//            std::cout << expired << " expired organisms still in population" << std::endl;
//        }
//        // --FOR DEBUGGING-- //
    }
    assert(this->populations_.size() <= this->carrying_capacity_);
}

// --- for trees etc ---

void Cell::sample_organisms(Species * sp_ptr,
        std::vector<const Organism *>& samples,
        PopulationCountType num_organisms) {
    assert(sp_ptr);
    OrganismPointers s = this->populations_[sp_ptr].sample_organism_ptrs(num_organisms, this->rng_);
    for (OrganismPointers::iterator og = s.begin(); og != s.end(); ++og) {
        sp_ptr->set_organism_label(*og, this->index_, this->x_, this->y_);
    }
    samples.insert(samples.end(), s.begin(), s.end());
}

void Cell::num_organisms(Species * species_ptr, PopulationCountType& num_females, PopulationCountType& num_males) const {
    this->populations_[species_ptr].num_organisms(num_females, num_males);
}

PopulationCountType Cell::num_organisms(Species * sp_ptr) const {
    PopulationCountType num_females = 0;
    PopulationCountType num_males = 0;
    this->num_organisms(sp_ptr, num_females, num_males);
    return num_females + num_males;
}

// Cell
///////////////////////////////////////////////////////////////////////////////
