///////////////////////////////////////////////////////////////////////////////
//
// GINGKO Biogeographical and Evolution Simulator.
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

#include "world.h"

using namespace gingko;

//! constructor: calls
World::World(unsigned long seed) 
    : species_pool_(),
      rng_(seed),
      landscape_(_species_pool, rng_) {
    this->current_generation_ = 0;    
}    

//! clean up species pool
World::~World() {
    for (std::vector<Species*>::iterator sp = this->species_pool_.begin();
            sp != this->species_pool_.end();
            ++sp) {
        delete *sp;            
    }            
}

// --- initialization and set up ---

//! Creates a new landscape.
void World::generate_landscape(CellIndexType size_x, CellIndexType size_y, unsigned num_environmental_factors) {
    this->num_fitness_factors_ = num_environmental_factors;
    this->landscape_.generate(size_x, size_y, num_environmental_factors);
}

//! Adds a new species definition to this world.
Species& World::new_species(const char* label) {
    Species* sp = new Species(this->species_pool_.size(),
                              label, 
                              this->num_fitness_factors_, 
                              this->rng_);
    this->species_pool_.push_back(sp);
    return *sp;
}

//! Populates the cell at (x,y) with organisms of the given species.
void World::seed_population(CellIndexType x, CellIndexType y, unsigned species_index, CellIndexType size) {
    this->landscape_.at(x, y).generate_new_organisms(species_index, size);
}

//! Populates the cell cell_index with organisms of the given species.
void World::seed_population(CellIndexType cell_index, unsigned species_index, CellIndexType size) {
    this->landscape_.at(cell_index).generate_new_organisms(species_index, size);
}

// --- species configuration ---


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

void World::run(int num_generations) {
    for ( ; num_generations >= 0; --num_generations, ++(this->current_generation_)) {
//##DEBUG##
DEBUG_BLOCK( std::cerr << "\n#### GENERATION " << this->current_generation_ << " ####\n\n"; )
DEBUG_BLOCK( this->landscape_.dump(std::cout); )        
        this->cycle();        
    }
}
