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

#include "organism.hpp"
#include "species.hpp"
#include "population.hpp"

using namespace ginkgo;

///////////////////////////////////////////////////////////////////////////////
// BreedingPopulation

// Returns pointers to all organisms.
std::vector<const Organism *> BreedingPopulation::organism_ptrs() {
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

// Selects pointers to random organisms across all species.
std::vector<const Organism *> BreedingPopulation::sample_organism_ptrs(PopulationCountType num_organisms) {
    std::vector<const Organism *> source = this->organism_ptrs();
    if (num_organisms <= source.size()) {
        std::vector<const Organism *> samples;
        samples.insert(samples.end(), source.begin(), source.begin() + num_organisms);
        return samples;
    } else {
        return source;
    }
}

void BreedingPopulation::clear() {
    this->females_.clear();
    this->males_.clear();
}

void BreedingPopulation::shuffle(RandomNumberGenerator& rng) {
    RandomPointer rp(rng);
    std::random_shuffle(females_.begin(), females_.end(), rp);
    std::random_shuffle(males_.begin(), males_.end(), rp);
}

void BreedingPopulation::swap(BreedingPopulation& other) {
    this->females_.swap(other.females_);
    this->males_.swap(other.males_);
}

// Removes all individuals marked for expiration.
void BreedingPopulation::purge_expired_organisms() {
    OrganismVector::iterator end_unexpired_females = std::remove_if(this->females_.begin(),
        this->females_.end(),
        std::mem_fun_ref(&Organism::is_expired));
    this->females_.erase(end_unexpired_females, this->females_.end());
    OrganismVector::iterator end_unexpired_males = std::remove_if(this->males_.begin(),
        this->males_.end(),
        std::mem_fun_ref(&Organism::is_expired));
    this->males_.erase(end_unexpired_males, this->males_.end());
}

// BreedingPopulation
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// BreedingPopulations

BreedingPopulations::BreedingPopulations(const SpeciesByLabel& species_pool, RandomNumberGenerator& rng):
        rng_ptr_(&rng) {
    for (SpeciesByLabel::const_iterator spi = species_pool.begin();
            spi != species_pool.end();
            ++spi) {
        BreedingPopulation pop;
        this->species_populations_.insert(std::make_pair((*spi).second, pop));
    }
}

void BreedingPopulations::purge_expired_organisms() {
    for (std::map<const Species *, BreedingPopulation >::iterator spi = this->species_populations_.begin();
            spi != this->species_populations_.end();
            ++spi) {
        (*spi).second.purge_expired_organisms();
    }
}

std::vector<const Organism *> BreedingPopulations::organism_ptrs() {
    std::vector<const Organism *> optrs;
    optrs.reserve(this->size());
    for (std::map<const Species *, BreedingPopulation >::iterator spi = this->species_populations_.begin();
            spi != this->species_populations_.end();
            ++spi) {
        BreedingPopulation& bp = (*spi).second;
        const std::vector<const Organism *>& poptrs = bp.organism_ptrs();
        optrs.insert(optrs.end(), poptrs.begin(), poptrs.end());
    }
    return optrs;
}

std::vector<const Organism *> BreedingPopulations::sample_organism_ptrs(PopulationCountType num_organisms) {
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

void BreedingPopulations::retain(PopulationCountType num_organisms) {
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

// BreedingPopulations
///////////////////////////////////////////////////////////////////////////////
