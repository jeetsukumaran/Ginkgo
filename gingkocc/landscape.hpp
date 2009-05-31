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

#include "biosys.hpp"
#include "cell.hpp"
#include <iostream>
#include <set>

namespace gingko {


///////////////////////////////////////////////////////////////////////////////	
// Landscape
/**
 * The geospatial framework.
 */
class Landscape {

    public:
    
        /** To store information on a single migration event: organism, dest */
        typedef std::pair<Organism, CellIndexType>   MigrationEvent;
        /** To store collection of single migration events */
        typedef std::vector<MigrationEvent>          MigrationEvents;
    
        // --- lifecycle and assignment ---

        /**
         * Constructs a landscape object, binding a species pool and a random
         * number generator.
         *
         * @param species   reference to map of Species objects
         * @param rng       reference to a RandomNumberGenerator
         */
        Landscape(const SpeciesByLabel& species, RandomNumberGenerator& rng);
        
        /** Destructor */
        ~Landscape();
        
        // --- initialization and set up ---

        /**
         * Create landscape data structure with specified dimensions.
         *
         * @param size_x                    x-dimension
         * @param size_y                    y-dimension
         * @param num_fitness_factors number of active fitness factors
         */
        void generate(CellIndexType size_x, CellIndexType size_y, unsigned num_fitness_factors); 
 
        // --- landscape access, control and mutation ---                                      

        /**
         * Globally fixes the maximum number of organisms that can occupy a 
         * cell.
         * Iterates through all cells, setting their carrying capacity.
         * @param carrying_capacity        maximum occupancy of each cell
         */
        void set_global_cell_carrying_capacity(unsigned long carrying_capacity);
        
        /**
         * Individually set the carrying capacity of cells.
         * @param cell_carrying_capacities   vector of numbers representing
         *                                   maximum cell occupancy on a per cell
         *                                   basis.
         */
        void set_carrying_capacities(const std::vector<unsigned long>& cell_carrying_capacities) {
            assert(cell_carrying_capacities.size() >= this->cells_.size());
            for (CellIndexType i = 0; i < this->cells_.size(); ++i) {
                (*this->cells_[i]).set_carrying_capacity(cell_carrying_capacities[i]);
            }
        }
        
        /**
         * Individually set the carrying capacity of cells.
         * @param cell_carrying_capacities   vector of numbers representing
         *                                   maximum cell occupancy on a per cell
         *                                   basis.
         */
        void set_carrying_capacities(const std::vector<long>& cell_carrying_capacities) {
            assert(cell_carrying_capacities.size() >= this->cells_.size());
            for (CellIndexType i = 0; i < this->cells_.size(); ++i) {
                (*this->cells_[i]).set_carrying_capacity(cell_carrying_capacities[i]);
            }
        }
        
        /**
         * Individually set the environment of cells.
         * @param factor_index       fitness factor index being set
         * @param cell_environments  vector of numbers representing the 
         *                           environmental value of cells
         */
        void set_environment(unsigned index, std::vector<long> cell_environments) {
            assert(cell_environments.size() >= this->cells_.size());
            for (CellIndexType i = 0; i < this->cells_.size(); ++i) {
                (*this->cells_[i]).set_environment_factor(index, cell_environments[i]);
            }        
        }
        
        // --- sampling and tree-building ---
        
        /** 
         * Loads vector of pointers to organisms of a particular species from
         * the given vector of cells. If num_organisms is 0 all organisms of 
         * that species are added, otherwise limited to num_organisms, sampled 
         * at random. If num_organisms exceeds the number of organisms of 
         * given species in the cell, then all the organisms are returned.
         * Organisms will be assigned labels based on their position.
         * @param sp                     pointer to Species object
         * @param num_organisms_per_cell number of organisms (0=all)
         * @param cell_indexes           indexes of cells from which to sample
         */
        void sample_organisms(Species * sp_ptr, 
                    unsigned long num_organisms_per_cell, 
                    const std::set<CellIndexType>& cell_indexes,
                    std::vector<const Organism *>& samples);
         
        /**
         * Count numbers of individuals of specified species in each cell.
         * @param sp_ptr    pointer to species
         * @param counts    vector of counts to populate
         */
        void count_organisms(Species * sp_ptr, std::vector<long>& counts) const;
           
        // --- cell access and spatial mapping ---
        
        /**
         * Returns cell at geospatial coordinates (x,y).
         *
         * @param   x   geospatial x-coordinate of the Cell
         * @param   y   geospatial y-coordinate of the Cell
         * @return      Cell object at (x,y)
         */
        Cell& operator()(CellIndexType x, CellIndexType y) {
            return *this->cells_[this->xy_to_index(x, y)];
        }
        
        /**
         * Returns cell at geospatial coordinates (x,y) with bounds checking.
         *
         * @param   x   geospatial x-coordinate of the Cell
         * @param   y   geospatial y-coordinate of the Cell
         * @return      Cell object at (x,y)
         */        
        Cell& at(CellIndexType x, CellIndexType y) {
            return *this->cells_.at(this->xy_to_index(x, y));
        }        
        
        /**
         * Returns cell at given vector index.
         *
         * @param   index   vector index of the Cell
         * @return          Cell object at given vector index
         */        
        Cell& operator[](CellIndexType index) {
            return *this->cells_[index];
        }         

        /**
         * Returns cell at given vector index with bounds checking.
         *
         * @param   index   vector index of the Cell
         * @return          Cell object at given vector index
         */                  
        Cell& at(CellIndexType index) {
            return *this->cells_.at(index);
        }
        
        /**
         * Returns geospatial x-coordinate for a given vector index.
         *
         * @param   index   index of element in vector
         * @return          geospatial x-coordinate corresponding to vector
         *                  index
         */
        CellIndexType index_to_x(CellIndexType index) const {
            return index % this->size_x_;          
        }
        
        /**
         * Returns geospatial y-coordinate for a given vector index.
         *
         * @param   index   index of element in vector
         * @return          geospatial y-coordinate corresponding to vector
         *                  index
         */        
        CellIndexType index_to_y(CellIndexType index) const {
            return static_cast<CellIndexType>(index / this->size_x_);            
        }
        
        /**
         * Returns geospatial (x, y) coordinates for a given vector index.
         *
         * @param    x  geospatial x-coordinate of the Cell
         * @param    y  geospatial y-coordinate of the Cell
         * @return      index of element in vector         
         */        
        CellIndexType xy_to_index(CellIndexType x, CellIndexType y) const {
            return (y * this->size_x_) + x;
        }
        
        /**
         * Returns length of vector of Cell objects.
         *
         * @return  number of Cell objects in this Landscape        
         */         
        CellIndexType size() const {
            return this->size_;
        }
        
        /**
         * Returns length of x-dimension of the geospatial framework 
         * superimposed on the vector of Cell objects.
         *
         * @return  length of x-dimension       
         */         
        CellIndexType size_x() const {
            return this->size_x_;
        }
        
        /**
         * Returns length of y-dimension of the geospatial framework 
         * superimposed on the vector of Cell objects.
         *
         * @return  length of y-dimension       
         */          
        CellIndexType size_y() const {
            return this->size_y_;
        }
        
        /**
         * Returns a Cell object that is adjacent to the cell of the given 
         * vector index when considered within the geospatial framework.         
         *
         * @param   i   index of the origin cell in the vector
         * @return      a random cell that is "next to" the origin cell in
         *              the geospatial framework
         */             
        CellIndexType random_neighbor(CellIndexType i) {
            static long x = 0;
            static long y = 0;
            
            x = this->index_to_x(i) + this->rng_.uniform_int(-1, 1); // to reflect: % this->size_x_; 
            y = this->index_to_y(i) + this->rng_.uniform_int(-1, 1); // to reflect: % this->size_y_; 
            if (x < 0) {
                x = 1;
            } else if (static_cast<CellIndexType>(x) >= this->size_x_) {
                x = this->size_x_ - 2;
            }
            if (y < 0) {
                y = 1;
            } else if (static_cast<CellIndexType>(y) >= this->size_y_) {
                y = this->size_y_ - 2;
            }
            return this->xy_to_index(x, y);
        }
        
        // --- migration and movement ---
        
        /**
         * Adds a migration event to be processed.        
         *
         * When simulating migration, we build up a collection of migrants
         * to be actually inserted into the destination cells at the end of the 
         * round.
         *
         * @param   organism    reference to organism that will be moved
         * @param   dest_cell   vector index of destination cell
         */ 
        void add_migrant(const Organism& organism, CellIndexType dest_cell) {
            this->migrants_.push_back(std::make_pair(organism, dest_cell));
        }
        
        /** Clear the collection of migration events. */
        void clear_migrants();               
        
        /** Actually process the migration events (move the organisms). */
        void process_migrants();
        
        // --- debugging ---
        
        // dump structure to std::cerr
        void debug_dump_structure(std::ostream& out = std::cerr);
        
    private:
    
        /** Size of the Landscape in the x-dimension in the geospatial framework */
        CellIndexType                   size_x_;
        /** Size of the Landscape in the y-dimension in the geospatial framework */
        CellIndexType                   size_y_;
        /** Total number of cells in this Landscape. */
        CellIndexType                   size_;
        /** Collection of cells in this Landscape. */
        std::vector<Cell*>              cells_;
        /** Collection migration events. */
        MigrationEvents                 migrants_;
        /** Reference to collection of pointers to Species objects in the World. */
        const SpeciesByLabel&           species_;
        /** Reference to RandomNumberGenerator of the World. */
        RandomNumberGenerator&          rng_;
};
// Landscape
//////////////////////////////////////////////////////////////////////////////


} // gingko namespace

#endif
