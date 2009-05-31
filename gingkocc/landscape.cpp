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

#include "landscape.hpp"

#include <iostream>
#include <iomanip>
#include "textutil.hpp"

using namespace gingko;

///////////////////////////////////////////////////////////////////////////////	
// Landscape 

// --- lifecycle and assignment --- 

Landscape::Landscape(const SpeciesByLabel& species, RandomNumberGenerator& rng)
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
    for (CellIndexType y = 0, index = 0; y < size_y; ++y) {
        for (CellIndexType x = 0; x < size_x; ++x, ++index) {
            Cell* cell = new Cell(index, x, y, num_fitness_factors, *this, this->species_, this->rng_);
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

// --- sampling for tree building ---

void Landscape::sample_organisms(Species * sp_ptr, 
        unsigned long num_organisms_per_cell,
        const std::set<CellIndexType>& cell_indexes,
        std::vector<const Organism *>& samples) {

    ////////////////
    // DEBUG CODE //
    ////////////////
    unsigned long old_size = 0;
    std::vector<unsigned> debug_results;
    std::vector<std::string> test_sample;
    ////////////////

    for (std::set<CellIndexType>::const_iterator ci = cell_indexes.begin();
            ci != cell_indexes.end();
            ++ci) {        
        assert(static_cast<unsigned long>(*ci) < this->cells_.size());
        this->cells_[*ci]->sample_organisms(sp_ptr, samples, num_organisms_per_cell);
        
        ////////////////
        // DEBUG CODE //
        ////////////////
        debug_results.push_back(samples.size() - old_size);
        if (samples.size() - old_size > 0) {
            std::string label = sp_ptr->get_organism_label(*samples.back());
            std::vector<std::string> parts = textutil::split(label, "_");
            test_sample.push_back(parts[1] + "," + parts[2]);
        } else {
            test_sample.push_back("--,--");
        }
        old_size = samples.size();
        ////////////////        
        
    }

    ////////////////
    // DEBUG CODE //
    ////////////////    
    std::cout << std::endl << "*** SAMPLING SCHEME ***" << std::endl;
    CellIndexType k = 0;    
    CellIndexType i = 0;
    CellIndexType x = 0;
    CellIndexType y = 0;
    for (y = 0; y < this->size_y_; ++y) {
        for (x = 0; x < this->size_x_; ++x, ++i) {
            std::cout << std::setw(8);
            if ( cell_indexes.find(i) != cell_indexes.end() ) {
                // std::cout << debug_results.at(k);
                std::cout << test_sample.at(k);
                k += 1;
            } else {
                std::cout << "-";
            }
        }
        std::cout << std::endl;
    }
    ////////////////    
    
}

void Landscape::count_organisms(Species * sp_ptr, std::vector<long>& counts) const {
    counts.reserve(this->cells_.size());
    for (std::vector<Cell *>::const_iterator ci = this->cells_.begin(); ci != this->cells_.end(); ++ci) {
        counts.push_back((*ci)->num_organisms(sp_ptr));
    }
}

// --- debugging ---

 // dump structure to std::cerr
void Landscape::debug_dump_structure(std::ostream& out) {
    out << std::endl << "Landscape Structure:";
    CellIndexType i = 0;
    unsigned long max_width_x = log10(static_cast<double>(this->size_x_)) + 1;
    unsigned long max_width_y = log10(static_cast<double>(this->size_y_)) + 1;
    for (std::vector<Cell*>::iterator ci = this->cells_.begin();
            ci != this->cells_.end();
            ++ci, ++i) {
        CellIndexType x = this->index_to_x(i);
        CellIndexType y = this->index_to_y(i);
        if (x == 0) {
            out << std::endl;
        }
        out << std::setfill('0') << std::setw(max_width_x) << (*ci)->get_x();
        out << ",";
        out << std::setfill('0') << std::setw(max_width_y) << (*ci)->get_y();
        out << " ";
    }
    out << std::endl << std::endl;
}

// Landscape
///////////////////////////////////////////////////////////////////////////////	
