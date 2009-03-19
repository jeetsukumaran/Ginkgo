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

#include "landscape.h"

#include <iostream>
#include <iomanip>

using namespace gingko;

///////////////////////////////////////////////////////////////////////////////	
// Landscape 

// --- lifecycle and assignment --- 

Landscape::Landscape(const SpeciesPointerVector& species, RandomNumberGenerator& rng)
    : species_(species),
      rng_(rng) {
    this->size_x_ = 0;
    this->size_y_ = 0;
    this->size_ = 0;
}

// clean up cells
Landscape::~Landscape() {
    for (std::vector<Cell*>::iterator cell = this->cells_.begin();
            cell != this->cells_.end();
            ++cell) {
        delete *cell; 
    }  
}

// --- initialization and set up ---

void Landscape::generate(CellIndexType size_x, CellIndexType size_y, unsigned num_fitness_factors) {
    this->size_x_ = size_x;
    this->size_y_ = size_y;
    this->size_ = size_x * size_y;
    this->cells_.reserve(this->size_);
    for (CellIndexType x = 0, index = 0; x < size_x; ++x) {
        for (CellIndexType y = 0; y < size_y; ++y, ++index) {
            Cell* cell = new Cell(index, num_fitness_factors, *this, this->species_, this->rng_);
            this->cells_.push_back(cell);
        }
    }
}
        
// --- migration and movement ---
void Landscape::clear_migrants() {
    this->migrants_.clear();
}


void Landscape::process_migrants() {
    for (MigrationEvents::iterator m = this->migrants_.begin();
            m != this->migrants_.end();
            ++m) {
        m->first.set_fitness(-1);   // invalidate cached fitness                    
        this->cells_.at(m->second)->insert_organism(m->first);                
    }
    this->migrants_.clear();
}

// --- debugging ---
unsigned long Landscape::dump(std::ostream& output) {
    unsigned long num = 0;
    unsigned long total = 0;
    for (CellIndexType y = 0; y < this->size_y_; ++y) {
        for (CellIndexType x = 0; x < this->size_x_; ++x) {
            num = this->operator()(x,y).num_organisms();
            total += num;
            output << std::setw(4) << num << " ";
        }
        output << std::endl;
    }
    output << "---\nTotal organisms: " << total << std::endl;            
    return total;
}

// Landscape
///////////////////////////////////////////////////////////////////////////////	
