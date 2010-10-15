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

#if !defined(GINKGO_POPULATION_H)
#define GINKGO_POPULATION_H

#include <vector>
#include <iterator>
#include <map>

#include "ginkgo_defs.hpp"
#include "randgen.hpp"
#include "organism.hpp"
#include "species.hpp"

namespace ginkgo {

// collection of Organism pointers
typedef Organism *                     OrganismPointer;
typedef std::vector<OrganismPointer>   OrganismPointers;

///////////////////////////////////////////////////////////////////////////////
// Census

/**
 * Count of individuals based on cell of origin (in turn, based on the
 * geographical location of the first (diploid) allele.
 */
class PopulationCensus {

    public:
        PopulationCensus();
        void log(const OrganismPointer optr);

    private:
        std::map<CellIndexType, PopulationCountType>  counts_;

};
// Census
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Breeding Population

/**
 * Manages a breeding pool.
 */
class BreedingPopulation {

    public:

        /**
         * Returns total number of individuals in this population.
         */
        PopulationCountType size() const {
            return this->males_.size() + this->females_.size();
        }

        void num_organisms(PopulationCountType& num_females, PopulationCountType& num_males) const {
            num_females = this->females_.size();
            num_males = this->males_.size();
        }

        /**
         * Returns reference to females.
         */
        OrganismPointers& females() {
            return this->females_;
        }

        /**
         * Returns reference to males.
         */
        OrganismPointers& males() {
            return this->males_;
        }

        /**
         * Adds a new organism to this breeding population.
         */
        void add(Organism* organism) {
            if (organism->is_male()) {
                this->males_.push_back(organism);
            } else {
                this->females_.push_back(organism);
            }
        }

        /**
         * Returns pointers to all organisms.
         */
        OrganismPointers get_organism_ptrs();

        /**
         * Inserts pointers to all organisms across all species into vector
         * of pointers passed as argument.
         */
        OrganismPointers& get_organism_ptrs(OrganismPointers& optrs);

        /**
         * Selects pointers to random organisms across all species.
         * @param   num_organisms   number of organisms
         */
        OrganismPointers sample_organism_ptrs(PopulationCountType num_organisms, RandomNumberGenerator& rng);

        /**
         * Removes and deallocates all organisms from population.
         */
        void clear();

        /**
         * Randomizes order of male and female organism vectors.
         */
        void shuffle(RandomNumberGenerator& rng);

        /**
         * Exchanges breeding vectors.
         */
        void swap(BreedingPopulation& other);

        /**
         * Removes and deallocates all organisms marked for expiration.
         */
        PopulationCountType purge_expired_organisms();

        /**
         * Counts organisms in this population, recording their cell
         * of origin.
         */
        PopulationCensus get_census();


        ///////////////////////////////////////////////////////////////////////
        // iteration support

        class iterator
        {
            public:
                typedef iterator self_type;
                typedef Organism* value_type;
                typedef std::forward_iterator_tag iterator_category;
                typedef int difference_type;
                iterator(OrganismPointers::iterator f_begin,
                         OrganismPointers::iterator f_end,
                         OrganismPointers::iterator m_begin,
                         OrganismPointers::iterator m_end)
                        : f_current_(f_begin),
                          f_end_(f_end),
                          m_current_(m_begin),
                          m_end_(m_end) {  }
                self_type operator++() {
                    self_type i = *this;
                    if (this->f_current_ != this->f_end_) {
                        ++(this->f_current_);
                    } else {
                        ++(this->m_current_);
                    }
                    return *this;
                }
                self_type operator++(int) {
                    ++(*this);
                    return *this;
                }
                value_type operator*() {
                    if (f_current_ != f_end_) {
                       return *(this->f_current_);
                    } else {
                       return *(this->m_current_);
                    }
                }
                value_type operator->() {
                    if (f_current_ != f_end_) {
                       return *(this->f_current_);
                    } else {
                       return *(this->m_current_);
                    }
                }
                bool operator==(const self_type& rhs) {
                    return this->f_current_ == rhs.f_current_
                        && this->f_end_ == rhs.f_end_
                        && this->m_current_ == rhs.m_current_
                        && this->m_end_ == rhs.m_end_;
                }
                bool operator!=(const self_type& rhs) {
                    return !(*this == rhs);
                }
            private:
                OrganismPointers::iterator       f_current_;
                OrganismPointers::iterator       f_end_;
                OrganismPointers::iterator       m_current_;
                OrganismPointers::iterator       m_end_;
        };

        iterator begin()
        {
            return iterator(this->females_.begin(), this->females_.end(), this->males_.begin(), this->males_.end());
        }

        iterator end()
        {
            return iterator(this->females_.end(), this->females_.end(), this->males_.end(), this->males_.end());
        }

        // iteration support
        ///////////////////////////////////////////////////////////////////////

    private:
        /** all females in the population */
        OrganismPointers      females_;
        /** all males in the population */
        OrganismPointers      males_;
};

///////////////////////////////////////////////////////////////////////////////
// BreedingPopulations
/**
 * Manages multiple breeding pools.
 */
class BreedingPopulations {

    public:

        /**
         * Constructor.
         */
        BreedingPopulations();

        /**
         * Returns a reference to the breeding population referenced by the
         * Species pointer 'sp'. Behavior is similar to [] of std::map.
         */
        BreedingPopulation& operator[](const Species * sp) {
            return this->species_populations_[sp];
        }

        /**
         * Returns a reference to the breeding population referenced by the
         * Species pointer 'sp'. Behavior is similar to [] of std::map.
         */
        const BreedingPopulation& operator[](const Species * sp) const {
            std::map<const Species *, BreedingPopulation >::const_iterator bpi = this->species_populations_.find(sp);
            return bpi->second;
        }

        /**
         * Adds an organism to the mix.
         */
        void add(Organism * organism_ptr) {
            this->species_populations_[organism_ptr->species_ptr()].add(organism_ptr);
        }

        /**
         * Returns total number of individuals across all populations.
         */
        PopulationCountType size() const {
            PopulationCountType s = 0;
            for (std::map<const Species *, BreedingPopulation >::const_iterator spi = this->species_populations_.begin();
                    spi != this->species_populations_.end();
                    ++spi) {
                s += (*spi).second.size();
            }
            return s;
        }

        /**
         * Removes all organisms from the population (but retains species slots).
         */
        void clear_organisms() {
            for (std::map<const Species *, BreedingPopulation >::iterator spi = this->species_populations_.begin();
                    spi != this->species_populations_.end();
                    ++spi) {
                (*spi).second.clear();
            }
        }


        /**
         * Removes all individuals marked for expiration.
         */
        PopulationCountType purge_expired_organisms();

        /**
         * Returns pointers to all organisms across all species.
         */
        OrganismPointers get_organism_ptrs();

        /**
         * Inserts pointers to all organisms across all species into vector
         * of pointers passed as argument.
         */
        OrganismPointers& get_organism_ptrs(OrganismPointers& optrs);

        /**
         * Selects pointers to random organisms across all species.
         * @param   num_organisms   number of organisms
         */
        OrganismPointers sample_organism_ptrs(PopulationCountType num_organisms);

        /**
         * Removes random organisms such that the total number of organisms is equal to num_organisms.
         * @param   num_organisms   number of organisms to retain
         */
        void retain(PopulationCountType num_organisms);

    private:
        /** the breeding pools */
        std::map<const Species *, BreedingPopulation >    species_populations_;

};

} // ginkgo namespace

#endif
