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

#if !defined(GINGKO_LANDSCAPE_H)
#define GINGKO_LANDSCAPE_H

#include "biosys.h"
#include "cell.h"
#include <iostream>

namespace gingko {


///////////////////////////////////////////////////////////////////////////////	
//! The landscape.
class Landscape {

    public:
    
        typedef std::pair<Organism, CellIndexType>   MigrationEvent;
        typedef std::vector<MigrationEvent>          MigrationEvents;
    
        // --- lifecycle and assignment ---        

        Landscape(const SpeciesPointerVector& species, RandomNumberGenerator& rng);
        ~Landscape();
        
        // --- initialization and set up ---

        void generate(CellIndexType size_x, CellIndexType size_y, unsigned num_environmental_factors); 
 
        // --- landscape access, control and mutation ---                                      

        void set_cell_carrying_capacity(unsigned long carrying_capacity);
        
        // --- cell access and spatial mapping ---
        
        Cell& operator()(CellIndexType x, CellIndexType y) {
            return *this->cells_[this->xy_to_index(x, y)];
        }
        Cell& operator[](CellIndexType index) {
            return *this->cells_[index];
        }            
        Cell& at(CellIndexType x, CellIndexType y) {
            return *this->cells_.at(this->xy_to_index(x, y));
        }
        Cell& at(CellIndexType index) {
            return *this->cells_.at(index);
        }
        CellIndexType index_to_x(CellIndexType index) const {
            return index % this->size_x_;          
        }
        CellIndexType index_to_y(CellIndexType index) const {
            return static_cast<CellIndexType>(index / this->size_x_);            
        }
        CellIndexType xy_to_index(CellIndexType x, CellIndexType y) const {
            return (y * this->size_x_) + x;
        }
        CellIndexType size() const {
            return this->size_;
        }
        CellIndexType size_x() const {
            return this->size_x_;
        }
        CellIndexType size_y() const {
            return this->size_y_;
        }
        CellIndexType random_neighbor(CellIndexType i) {
            static CellIndexType x = 0;
            static CellIndexType y = 0;
            
            x = this->index_to_x(i) + this->rng_.uniform_int(-1, 1); // to reflect: % this->size_x_; 
            y = this->index_to_y(i) + this->rng_.uniform_int(-1, 1); // to reflect: % this->size_y_; 
            if (x >= this->size_x_) {
                x = this->size_x_ - 1;
            } else if (x < 0) {
                x = 0;
            }
            if (y >= this->size_y_) {
                y = this->size_y_ - 1;
            } else if (y < 0) {
                y = 0;
            }
            return this->xy_to_index(x, y);
        }
        
        // --- migration and movement ---

        void add_migrant(const Organism& organism, CellIndexType dest_cell) {
            this->migrants_.push_back(std::make_pair(organism, dest_cell));
        }
        void clear_migrants();               
        void process_migrants();
        
        // --- debugging ---
        
        unsigned long dump(std::ostream& output = std::cout);       
    
    private:
        CellIndexType                   size_x_;                // size of the landscape in the x dimension
        CellIndexType                   size_y_;                // size of the landscape in the y dimension        
        CellIndexType                   size_;                  // == x * y, cached here
        std::vector<Cell*>              cells_;                 // cells of the landscape
        
        // cache to collect migrating organisms over a round of migration
        MigrationEvents                 migrants_;
        
        const SpeciesPointerVector&     species_;               // species pool
        RandomNumberGenerator&          rng_;                   // random number generator        
};
// Landscape
//////////////////////////////////////////////////////////////////////////////

} // gingko namespace

#endif
