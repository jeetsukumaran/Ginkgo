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

#include "organism.hpp"
#include "species.hpp"
#include "population.hpp"

using namespace ginkgo;

// local helper functions

PopulationCountType purge_expired_organisms_(OrganismPointers& organism_ptrs) {
    PopulationCountType num_purged = 0;
	OrganismMemoryManager& organism_memory_manager = OrganismMemoryManager::get_instance();
	std::vector<Organism *> to_keep_ptrs;
	to_keep_ptrs.reserve(organism_ptrs.size());
	for (OrganismPointers::iterator i = organism_ptrs.begin(); i != organism_ptrs.end(); ++i) {
		Organism * org = *i;
		if (org->is_expired()) {
		    ++num_purged;
            organism_memory_manager.deallocate(org);
		} else {
			to_keep_ptrs.push_back(org);
		}
	}
	organism_ptrs.clear();
	to_keep_ptrs.swap(organism_ptrs);
	return num_purged;
}

void deallocate_and_clear_organisms_(OrganismPointers& organism_ptrs) {
	OrganismMemoryManager& organism_memory_manager = OrganismMemoryManager::get_instance();
    for (OrganismPointers::iterator oi = organism_ptrs.begin();
            oi != organism_ptrs.end();
            ++oi) {
        organism_memory_manager.deallocate(*oi);
    }
    organism_ptrs.clear();
}

///////////////////////////////////////////////////////////////////////////////
// Census

OrganismProvenances::OrganismProvenances() { }

void OrganismProvenances::log(const OrganismPointer optr) {
    // TODO: make this safe (i.e. if no diploid nodes, default to haploid node.
    GenealogyNode * allele1 = optr->get_diploid_node1(0);
    CellIndexType origin_cell_idx =  allele1->get_cell_index();
    this->counts_[origin_cell_idx] += 1;
}

PopulationCountType OrganismProvenances::get_count(CellIndexType cell_index) const {
    std::map<CellIndexType, PopulationCountType>::const_iterator ci = this->counts_.find(cell_index);
    if (ci == this->counts_.end()) {
        return 0;
    } else {
        return ci->second;
    }
}

// Census
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// BreedingPopulation

// Returns pointers to all organisms.
OrganismPointers BreedingPopulation::get_organism_ptrs() {
    OrganismPointers o;
    return this->get_organism_ptrs(o);
}

OrganismPointers& BreedingPopulation::get_organism_ptrs(OrganismPointers& optrs) {
    optrs.reserve(optrs.size() + this->females_.size() + this->males_.size());
    optrs.insert(optrs.end(), this->females_.begin(), this->females_.end());
    optrs.insert(optrs.end(), this->males_.begin(), this->males_.end());
    return optrs;
}

OrganismPointers BreedingPopulation::sample_organism_ptrs(PopulationCountType num_organisms, RandomNumberGenerator& rng) {
    OrganismPointers source;
    this->get_organism_ptrs(source);
    if (num_organisms <= source.size()) {
        OrganismPointers samples;
        RandomPointer rp(rng);
        std::random_shuffle(source.begin(), source.end(), rp);
        samples.insert(samples.end(), source.begin(), source.begin() + num_organisms);
        return samples;
    } else {
        return source;
    }
}

void BreedingPopulation::clear() {
    deallocate_and_clear_organisms_(this->females_);
    deallocate_and_clear_organisms_(this->males_);
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
PopulationCountType BreedingPopulation::purge_expired_organisms() {
    PopulationCountType num_purged = 0;
    num_purged += purge_expired_organisms_(this->females_);
    num_purged += purge_expired_organisms_(this->males_);
    return num_purged;
}

OrganismProvenances BreedingPopulation::get_organism_provenances() {
    OrganismProvenances provenances;
    OrganismPointers source;
    this->get_organism_ptrs(source);
    for (OrganismPointers::const_iterator optr = source.begin(); \
            optr != source.end(); \
            ++optr) {
        provenances.log(*optr);
    }
    return provenances;
}
// BreedingPopulation
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// BreedingPopulations

// constructor
BreedingPopulations::BreedingPopulations() {
    SpeciesRegistry& species_registry = SpeciesRegistry::get_instance();
    for (SpeciesRegistry::iterator spi = species_registry.begin();
            spi != species_registry.end();
            ++spi) {
        BreedingPopulation pop;
        this->species_populations_.insert(std::make_pair(*spi, pop));
    }
}

// removes organisms flagged for expiration
PopulationCountType BreedingPopulations::purge_expired_organisms() {
    PopulationCountType num_purged = 0;
    for (std::map<const Species *, BreedingPopulation >::iterator spi = this->species_populations_.begin();
            spi != this->species_populations_.end();
            ++spi) {
        num_purged += (*spi).second.purge_expired_organisms();
    }
    return num_purged;
}

// returns pointers to all organisms
OrganismPointers BreedingPopulations::get_organism_ptrs() {
    OrganismPointers o;
    return this->get_organism_ptrs(o);
}

OrganismPointers& BreedingPopulations::get_organism_ptrs(OrganismPointers& optrs) {
    optrs.reserve(optrs.size() + this->size());
    for (std::map<const Species *, BreedingPopulation >::iterator spi = this->species_populations_.begin();
            spi != this->species_populations_.end();
            ++spi) {
        BreedingPopulation& bp = (*spi).second;
        const OrganismPointers& poptrs = bp.get_organism_ptrs();
        optrs.insert(optrs.end(), poptrs.begin(), poptrs.end());
    }
    return optrs;
}

// returns random sample of pointers to organisms
OrganismPointers BreedingPopulations::sample_organism_ptrs(PopulationCountType num_organisms) {
    OrganismPointers source;
    this->get_organism_ptrs(source);
    RandomPointer rp(RandomNumberGenerator::get_instance());
    std::random_shuffle(source.begin(), source.end(), rp);
    if (num_organisms <= source.size()) {
        OrganismPointers samples;
        samples.insert(samples.end(), source.begin(), source.begin() + num_organisms);
        return samples;
    } else {
        return source;
    }
}

// removes all but specified number of organisms randomly
void BreedingPopulations::retain(PopulationCountType num_organisms) {
    long surplus = this->size() - num_organisms;
    if (surplus <= 0 ) {
        return;
    }
    OrganismPointers to_remove = this->sample_organism_ptrs(surplus);
    for (OrganismPointers::iterator i = to_remove.begin();
            i != to_remove.end();
            ++i) {
        (*i)->set_expired();
    }
    this->purge_expired_organisms();
}

// BreedingPopulations
///////////////////////////////////////////////////////////////////////////////
