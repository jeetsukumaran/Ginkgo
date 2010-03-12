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
#include "randgen.hpp"

using namespace ginkgo;

///////////////////////////////////////////////////////////////////////////////
// Species

Species::Species(const std::string& label)
    : label_(label),
      num_fitness_traits_(MAX_FITNESS_TRAITS),
      is_fixed_movement_capacity_(true),
      movement_capacity_(1),
      global_selection_strength_(1.0),
      rng_(RandomNumberGenerator::get_instance()),
      organism_label_index_(0),
      organism_memory_manager_(OrganismMemoryManager::get_instance()) {
    this->mean_reproductive_rate_ = 6;
    this->selection_weights_.assign(this->num_fitness_traits_, 1);
    this->fitness_trait_inheritance_sd_.assign(this->num_fitness_traits_, 0.7071068);
    memset(this->default_fitness_trait_genotypes_, 0,
        MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));
    this->movement_capacity_ = 1;
}

Species::Species(const std::string& label,
                 unsigned num_fitness_traits,
                 float global_selection_strength)
    : label_(label),
      num_fitness_traits_(num_fitness_traits),
      is_fixed_movement_capacity_(true),
      movement_capacity_(1),
      global_selection_strength_(global_selection_strength),
      rng_(RandomNumberGenerator::get_instance()),
      organism_label_index_(0),
      organism_memory_manager_(OrganismMemoryManager::get_instance()) {
    this->mean_reproductive_rate_ = 6;
    this->selection_weights_.assign(this->num_fitness_traits_, 1);
    this->fitness_trait_inheritance_sd_.assign(this->num_fitness_traits_, 0.7071068);
    memset(this->default_fitness_trait_genotypes_, 0,
        MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));
    this->movement_capacity_ = 1;
}
// Species
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// SpeciesRegistry

SpeciesRegistry SpeciesRegistry::instance_;

SpeciesRegistry::~SpeciesRegistry() {
    for (std::list<Species *>::iterator spi = this->species_list_.begin();
            spi != this->species_list_.end();
            ++spi) {
        assert(*spi != NULL);
        delete *spi;
    }
}

Species * SpeciesRegistry::new_species(const std::string& label) {
    if (this->label_species_map_.find(label) != this->label_species_map_.end()) {
        throw std::runtime_error("Species with label '" + label + "' already defined");
    }
    Species *  sp = new Species(label);
    this->species_list_.push_back(sp);
    this->label_species_map_.insert(std::make_pair(label, sp));
    return sp;
}

// SpeciesRegistry
///////////////////////////////////////////////////////////////////////////////
