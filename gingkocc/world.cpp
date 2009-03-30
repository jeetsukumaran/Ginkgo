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

#include "world.hpp"
#include <iostream>
#include <string>
#include <map>
#include <utility>

using namespace gingko;

// constructor
World::World() 
    : species_(),
      rng_(),
      landscape_(species_, rng_) {
    this->current_generation_ = 0;    
}

// constructor
World::World(unsigned long seed) 
    : species_(),
      rng_(seed),
      landscape_(species_, rng_) {
    this->current_generation_ = 0;    
}    

// clean up species pool
World::~World() {
    for (SpeciesByLabel::iterator sp = this->species_.begin();
            sp != this->species_.end();
            ++sp) {
        assert (sp->second != NULL);            
        delete sp->second;
        sp->second = NULL;
        this->species_.erase(sp);
    }            
}

// --- initialization and set up ---

// Creates a new landscape.
void World::generate_landscape(CellIndexType size_x, CellIndexType size_y, unsigned num_fitness_factors) {
    this->num_fitness_factors_ = num_fitness_factors;
    this->landscape_.generate(size_x, size_y, num_fitness_factors);
}

// Adds a new species definition to this world.
Species& World::new_species(const std::string& label) {
    Species* sp = new Species(label, 
                              this->num_fitness_factors_, 
                              this->rng_);
    this->species_.insert(std::make_pair(std::string(label), sp));
    return *sp;
}

// Populates the cell at (x,y) with organisms of the given species.
void World::seed_population(CellIndexType x, CellIndexType y, const std::string& species_label, unsigned long size) {
    assert(this->species_.find(species_label) != this->species_.end());
    this->landscape_.at(x, y).generate_new_organisms(this->species_[species_label], size);
}

// Populates the cell cell_index with organisms of the given species.
void World::seed_population(CellIndexType cell_index, const std::string& species_label, unsigned long size) {
    assert(this->species_.find(species_label) != this->species_.end());
    this->landscape_.at(cell_index).generate_new_organisms(this->species_[species_label], size);
}

// --- event handlers ---

void World::add_world_settings(unsigned long generation, const WorldSettings& world_settings) {
    this->world_settings_[generation] = world_settings;
}


// --- simulation cycles ---
void World::cycle() {

// Results in inflated population: the migrants get distributed after the 
// competition phase, resulting in a artificially (>> carrying capacity) 
// boosted population when entering the next generation's reproduction phase.
// This leads to a standing population at the end of each generation sometimes
// an order or more of magnitude above the carrying capacity of the cell.
//     for (CellIndexType i = this->landscape_.size()-1; i >= 0; --i) {
//         this->landscape_[i].reproduction(); 
//         this->landscape_[i].migration();
//         this->landscape_[i].survival();
//         this->landscape_[i].competition();
//     }
//     this->landscape_.process_migrants();
    ++this->current_generation_;

    std::cerr << "GENERATION " << this->current_generation_ << std::endl;

    for (CellIndexType i = this->landscape_.size()-1; i >= 0; --i) {
        this->landscape_[i].reproduction(); 
        this->landscape_[i].migration();
    }
    this->landscape_.process_migrants();    
    for (CellIndexType i = this->landscape_.size()-1; i >= 0; --i) {    
        this->landscape_[i].survival();
        this->landscape_[i].competition();        
    }    
}

void World::run(unsigned long num_generations) {    
    for (unsigned long g=0; g < num_generations; ++g) {
        this->cycle();        
    }
}
