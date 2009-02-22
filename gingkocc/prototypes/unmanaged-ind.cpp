///////////////////////////////////////////////////////////////////////////////
//
// unmanaged-ind.cpp
//
// Performance evaluation of population evolution under non-managed allocation
// of individuals.
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

#include <vector>
#include <list>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
// POPULATION ECOLOGY AND GENETICS

class Species;
class Population;
class Cell;

typedef std::vector<float> Genotype;

/// A single organism of a population of a particular species.
/// Responsible for tracking (non-neutral) genotype and neutral marker 
/// histories. 
class Individual {
    public:
        
        /// instantiates individual with given genotype
        Individual(const Genotype& g)
            : genotype(g) {
        }            
                
    private:
        Genotype genotype;   /// non-neutral genotype: maps to fitness phenotype
        
        /// returns individual's genotype; private to disable copying
        const Genotype& get_genotype() const {
            return this->genotype;
        }

}; // Individual

/// A single population of a particular species.
/// Responsible for managing collections of individuals, and relating them to 
/// their species.
class Population {
    public:
        Population(const Species& sp)
            : species(sp) {
        }
    
    private:
        const Species&    species; // the species to which this population belongs
        Cell*             cell;    // the current location of this population
        std::vector<Individual> individuals; // members of this population        
}; // Population

/// A collection of Populations sharing the same ecologies (e.g. movement, 
/// fitness/survival functions, breeding pool)
class Species {
    public:
    
//     private:
//         std::list<Population*> populations;
}; // Species

typedef std::vector<Species> SpeciesContainer;
typedef std::vector<Species>::iterator SpeciesIterator;

///////////////////////////////////////////////////////////////////////////////
// SPATIAL FRAMEWORK

/// The fundamental atomic spatial unit of the world.
class Cell {
    public:
    
        static void set_carrying_capacity(Cell& cell) {
            cell.set_carrying_capacity(4);
        }
        
        void set_carrying_capacity(int cc) {
            this->carrying_capacity = cc;
        }
    
    private:
        int     carrying_capacity;
//         std::vector<Population> populations;    

}; // Cell

typedef std::vector<Cell> Cells;
typedef std::vector<Cell>::iterator CellIterator;


/// The world.
class World {

    public:
        World(int dim_x, int dim_y, const SpeciesContainer& spp);
        void initialize(int dim_x, int dim_y);
        void set_carrying_capacity(int carrying_capacity);
        
    private:
        int                 dim_x;
        int                 dim_y;
        std::vector<Cell>   cells;
        SpeciesContainer    species;   
        
}; // World

// constructor: calls landscape initializer
World::World(int dim_x, int dim_y, const SpeciesContainer& spp)
    : species(spp) {   
    this->initialize(dim_x, dim_y);
}

void World::initialize(int dim_x, int dim_y) {
    this->dim_x = dim_x;
    this->dim_y = dim_y;
    int num_cells = dim_x * dim_y;
    // this->cells.reserve(num_cells);
    for (int i=0; i < num_cells; ++i) {
        this->cells.push_back(Cell());
    }
}

void World::set_carrying_capacity(int carrying_capacity) {
    // std::for_each(this->cells.begin(), this->cells.end(), Cell::set_carrying_capacity());
    carrying_capacity++;
}

///////////////////////////////////////////////////////////////////////////////
// TEST

int main(int argc, char * argv[]) {
    if (argc < 4) {
        std::cout << "usage: " << argv[0] <<  " <DIM-X> <DIM-Y> <CELL-CARRYING-CAPACITY>\n";
        exit(1);
    }
    int dim_x = atoi(argv[1]);
    int dim_y = atoi(argv[2]);
    //int cc = atoi(argv[3]);
    
    Species sp1;
    SpeciesContainer spp;
    spp.push_back(sp1);
   
	World world(dim_x, dim_y, spp);
	return 0;
}







