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

#if !defined(GINKGO_BIOSYS_H)
#define GINKGO_BIOSYS_H

#include <algorithm>
#include <cassert>
#include <vector>
#include <iterator>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <cstring>

#include "ginkgo_defs.hpp"
#include "randgen.hpp"
#include "organism.hpp"
#include "species.hpp"

namespace ginkgo {

///////////////////////////////////////////////////////////////////////////////
// Breeding Population
/**
 * Manages a breeding pool.
 */
class BreedingPopulation {

    public:

        /**
         * Adds a new organism to this breeding population.
         */
        void add(const Organism& organism) {
            if (organism.is_male()) {
                this->males_.push_back(organism);
            } else {
                this->females_.push_back(organism);
            }
        }

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
        OrganismVector& females() {
            return this->females_;
        }

        /**
         * Returns reference to males.
         */
        OrganismVector& males() {
            return this->males_;
        }

        /**
         * Returns pointers to all organisms.
         */
        std::vector<const Organism *> organism_ptrs() {
            std::vector<const Organism *> optrs;
            optrs.reserve(this->size());
            for (OrganismVector::const_iterator ov = this->females_.begin();
                    ov != this->females_.end();
                    ++ov) {
                optrs.push_back(&*(ov));
            }
            for (OrganismVector::const_iterator ov = this->males_.begin();
                    ov != this->males_.end();
                    ++ov) {
                optrs.push_back(&*(ov));
            }
            return optrs;
        }

        /**
         * Selects pointers to random organisms across all species.
         * @param   num_organisms   number of organisms
         */
        std::vector<const Organism *> sample_organism_ptrs(PopulationCountType num_organisms) {
            std::vector<const Organism *> source = this->organism_ptrs();
            if (num_organisms <= source.size()) {
                std::vector<const Organism *> samples;
                samples.insert(samples.end(), source.begin(), source.begin() + num_organisms);
                return samples;
            } else {
                return source;
            }
        }

        void clear() {
            this->females_.clear();
            this->males_.clear();
        }

        void shuffle(RandomNumberGenerator& rng) {
            RandomPointer rp(rng);
            std::random_shuffle(females_.begin(), females_.end(), rp);
            std::random_shuffle(males_.begin(), males_.end(), rp);
        }

        void swap(BreedingPopulation& other) {
            this->females_.swap(other.females_);
            this->males_.swap(other.males_);
        }

        /**
         * Removes all individuals marked for expiration.
         */
        void purge_expired_organisms() {
            OrganismVector::iterator end_unexpired_females = std::remove_if(this->females_.begin(),
                this->females_.end(),
                std::mem_fun_ref(&Organism::is_expired));
            this->females_.erase(end_unexpired_females, this->females_.end());
            OrganismVector::iterator end_unexpired_males = std::remove_if(this->males_.begin(),
                this->males_.end(),
                std::mem_fun_ref(&Organism::is_expired));
            this->males_.erase(end_unexpired_males, this->males_.end());
        }

        ///////////////////////////////////////////////////////////////////////
        // iterators

        class iterator
        {
            public:
                typedef iterator self_type;
                typedef Organism value_type;
                typedef Organism& reference;
                typedef Organism* pointer;
                typedef std::forward_iterator_tag iterator_category;
                typedef int difference_type;
                iterator(OrganismVector::iterator f_begin,
                         OrganismVector::iterator f_end,
                         OrganismVector::iterator m_begin,
                         OrganismVector::iterator m_end)
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
                reference operator*() {
                    if (f_current_ != f_end_) {
                       return *(this->f_current_);
                    } else {
                       return *(this->m_current_);
                    }
                }
                pointer operator->() {
                    if (f_current_ != f_end_) {
                        return &(*(this->f_current_));
                    } else {
                        return &(*(this->m_current_));
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
                OrganismVector::iterator       f_current_;
                OrganismVector::iterator       f_end_;
                OrganismVector::iterator       m_current_;
                OrganismVector::iterator       m_end_;
        };

        iterator begin()
        {
            return iterator(this->females_.begin(), this->females_.end(), this->males_.begin(), this->males_.end());
        }

        iterator end()
        {
            return iterator(this->females_.end(), this->females_.end(), this->males_.end(), this->males_.end());
        }

        // iterators
        ///////////////////////////////////////////////////////////////////////

    private:
        /** all females in the population */
        OrganismVector      females_;
        /** all males in the population */
        OrganismVector      males_;
};

///////////////////////////////////////////////////////////////////////////////
// Breeding Population
/**
 * Manages multiple breeding pools.
 */
class BreedingPopulations {

    public:

        /**
         * Constructor, takes reference to species pool map and rng.
         */
        BreedingPopulations(const SpeciesByLabel& species_pool, RandomNumberGenerator& rng):
                rng_ptr_(&rng) {
            for (SpeciesByLabel::const_iterator spi = species_pool.begin();
                    spi != species_pool.end();
                    ++spi) {
                BreedingPopulation pop;
                this->species_populations_.insert(std::make_pair((*spi).second, pop));
            }
        }

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

//        const BreedingPopulations& operator=(const BreedingPopulations bp) {
//            return *this;
//        }

        /**
         * Adds an organism to the mix.
         */
        void add(const Organism& organism) {
            this->species_populations_[organism.species_ptr()].add(organism);
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
         * Removes all individuals marked for expiration.
         */
        void purge_expired_organisms() {
            for (std::map<const Species *, BreedingPopulation >::iterator spi = this->species_populations_.begin();
                    spi != this->species_populations_.end();
                    ++spi) {
                (*spi).second.purge_expired_organisms();
            }
        }

        /**
         * Returns pointers to all organisms across all species.
         */
        std::vector<const Organism *> organism_ptrs() {
            std::vector<const Organism *> optrs;
            optrs.reserve(this->size());
            for (std::map<const Species *, BreedingPopulation >::iterator spi = this->species_populations_.begin();
                    spi != this->species_populations_.end();
                    ++spi) {
//                const Species * sp = (*spi).first;
                BreedingPopulation& bp = (*spi).second;
                const std::vector<const Organism *>& poptrs = bp.organism_ptrs();
                optrs.insert(optrs.end(), poptrs.begin(), poptrs.end());
            }
            return optrs;
        }

        /**
         * Selects pointers to random organisms across all species.
         * @param   num_organisms   number of organisms
         */
        std::vector<const Organism *> sample_organism_ptrs(PopulationCountType num_organisms) {
            std::vector<const Organism *> source = this->organism_ptrs();
            RandomPointer rp(*(this->rng_ptr_));
            std::random_shuffle(source.begin(), source.end(), rp);
            if (num_organisms <= source.size()) {
                std::vector<const Organism *> samples;
                samples.insert(samples.end(), source.begin(), source.begin() + num_organisms);
                return samples;
            } else {
                return source;
            }
        }

        /**
         * Removes random organisms such that the total number of organisms is equal to num_organisms.
         * @param   num_organisms   number of organisms to retain
         */
        void retain(PopulationCountType num_organisms) {
            if (num_organisms <= this->size()) {
                return;
            }
            std::vector<const Organism *> source = this->organism_ptrs();
            RandomPointer rp(*(this->rng_ptr_));
            std::random_shuffle(source.begin(), source.end(), rp);
            std::map< const Species *, BreedingPopulation > new_pop;
            for (std::vector<const Organism *>::const_iterator oi = source.begin();
                    oi <= source.end();
                    ++oi) {
                new_pop[(**oi).species_ptr()].add(**oi);
            }
            this->species_populations_ = new_pop;
        }

    private:
        /** the breeding pools */
        std::map<const Species *, BreedingPopulation >    species_populations_;
        /** random number genereator */
        RandomNumberGenerator *                           rng_ptr_;

};

} // ginkgo namespace

#endif
