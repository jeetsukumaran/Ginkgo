///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Phylogeographical Evolution Simulator.
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

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <list>

#include "cell.hpp"
#include "landscape.hpp"
#include "randgen.hpp"
#include "world.hpp"

using namespace ginkgo;

///////////////////////////////////////////////////////////////////////////////
// Cell

// --- lifecycle and assignment ---

Cell::Cell(CellIndexType index,
           CellIndexType x,
           CellIndexType y,
           unsigned num_fitness_traits,
           Landscape& landscape)
        : index_(index),
          x_(x),
          y_(y),
          carrying_capacity_(0),
          num_fitness_traits_(num_fitness_traits),
          landscape_(landscape),
          species_registry_(SpeciesRegistry::get_instance()),
          rng_(RandomNumberGenerator::get_instance()),
          organism_memory_manager_(OrganismMemoryManager::get_instance()) {
    memset(this->fitness_trait_optimum_, 0, MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));
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
    Cell temp_cell(this->index_, this->x_, this->y_, this->num_fitness_traits_, this->landscape_);
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

    for (SpeciesRegistry::iterator spi = this->species_registry_.begin(); spi != this->species_registry_.end(); ++spi) {
        Species * sp = *spi;
        BreedingPopulation& pop = this->populations_[sp];

        // TODO: instead of handles to the male and female vectors,
        //       BreedingPopulation should return iterators

        OrganismPointers& female_ptrs = pop.females();
        OrganismPointers& male_ptrs = pop.males();
        if ( (female_ptrs.size() > 0) && (male_ptrs.size() > 0)) {
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
    for (SpeciesRegistry::const_iterator spi = this->species_registry_.begin();
            spi != this->species_registry_.end();
            ++spi) {
        Species * sp = *spi;
        BreedingPopulation& current_pop = this->populations_[sp];
        BreedingPopulation remaining_pop;
        for (BreedingPopulation::iterator oi = current_pop.begin(); oi != current_pop.end(); ++oi) {
            Organism * og_ptr = *oi;
//            if (this->rng_.uniform_01() <= sp->movement_probability(this->index_)) {
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
//            } // movement probability: if (this->rng_.uniform_01() <= sp->movement_probability(this->index_))
        }
        this->populations_[sp] = remaining_pop;
    }
}

PopulationCountType Cell::survival() {
    PopulationCountType num_survivors = 0;
    for (SpeciesRegistry::const_iterator spi = this->species_registry_.begin();
            spi != this->species_registry_.end();
            ++spi) {
        Species * sp = *spi;
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

template <class T>
void set_unexpired_(T& organism_ptrs) {
    for (typename T::iterator i = organism_ptrs.begin();
            i != organism_ptrs.end();
            ++i) {
        (*i)->set_unexpired();
    }
}

void Cell::competition() {

    if (this->populations_.size() <=  this->carrying_capacity_) {
        return;
    }
    if (this->carrying_capacity_ == 0) {
        this->populations_.clear_organisms();
        return;
    }
    OrganismPointers original_pop;
    this->populations_.get_organism_ptrs(original_pop);

    std::map<double, std::list<Organism *> > fitness_organisms_map;
    OrganismPointers::const_iterator pop_iter = original_pop.begin();

    // add initial batch
    for (PopulationCountType num_added=0;
            num_added < this->carrying_capacity_;
            ++pop_iter, ++num_added) {
        (*pop_iter)->set_expired();
        fitness_organisms_map[(*pop_iter)->get_fitness()].push_back(*pop_iter);
    }

    PopulationCountType num_in_map = this->carrying_capacity_;
    std::map<double, std::list<Organism *> >::iterator lowest_fitness_iter = fitness_organisms_map.begin();
    double lowest_stored_fitness_score = lowest_fitness_iter->first;
    PopulationCountType num_in_lowest_fitness = lowest_fitness_iter->second.size();

    for (; pop_iter != original_pop.end(); ++pop_iter) {
        (*pop_iter)->set_expired();
        const double curr_fitness = (*pop_iter)->get_fitness();
        if (curr_fitness > lowest_stored_fitness_score) {
            fitness_organisms_map[curr_fitness].push_back(*pop_iter);
            if (num_in_map - num_in_lowest_fitness == this->carrying_capacity_ - 1) {
                num_in_map = this->carrying_capacity_;
                fitness_organisms_map.erase(lowest_fitness_iter);
                lowest_fitness_iter = fitness_organisms_map.begin();
                lowest_stored_fitness_score = lowest_fitness_iter->first;
                num_in_lowest_fitness = lowest_fitness_iter->second.size();
            } else {
                ++num_in_map;
            }
        } else if (curr_fitness == lowest_stored_fitness_score) {
            ++num_in_lowest_fitness;
            ++num_in_map;
            lowest_fitness_iter->second.push_back(*pop_iter);
        }
    }

    PopulationCountType total_unexpired = 0;
    if (num_in_map > this->carrying_capacity_) {
        std::map<double, std::list<Organism *> >::iterator mi = fitness_organisms_map.begin();
        ++mi;
        for (; mi != fitness_organisms_map.end();
                ++mi) {
            set_unexpired_(mi->second);
            total_unexpired += mi->second.size();
        }

        // TODO: !!!REFACTOR FOR EFFICIENCY!!!
        // randomly sample list
        // not crazy about current approach ...
        lowest_fitness_iter = fitness_organisms_map.begin();
        std::vector<Organism *> final(lowest_fitness_iter->second.begin(), lowest_fitness_iter->second.end());
        RandomPointer rp(this->rng_);
        std::random_shuffle(final.begin(), final.end(), rp);
        std::vector<Organism *>::iterator i = final.begin();
        while ( (total_unexpired < (this->carrying_capacity_) ) && i != final.end()) {
            (*i)->set_unexpired();
            ++total_unexpired;
            ++i;
        }
    } else {
        for (std::map<double, std::list<Organism *> >::iterator mi = fitness_organisms_map.begin();
                mi != fitness_organisms_map.end();
                ++mi) {
            set_unexpired_(mi->second);
        }
    }
    this->purge_expired_organisms();
    assert(this->populations_.size() <= this->carrying_capacity_);
}

// --- for trees etc ---

void Cell::sample_organisms(Species * sp_ptr,
        std::vector<const Organism *>& samples,
        PopulationCountType num_organisms) {
    assert(sp_ptr);
    OrganismPointers s;
    if (num_organisms > 0) {
        s = this->populations_[sp_ptr].sample_organism_ptrs(num_organisms, this->rng_);
    } else {
        this->populations_[sp_ptr].get_organism_ptrs(s);
    }
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
