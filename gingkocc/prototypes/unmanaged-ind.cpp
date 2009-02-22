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
#include <cstring>
#include <list>
#include <iostream>
#include <string>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// SUPPORT FUNCTIONS

double random_uniform() {
    return (double)rand()/double(RAND_MAX);
}

int poisson_rv(int rate) {
    const int MAX_EXPECTATION = 64;
    if (rate > MAX_EXPECTATION) {
        double r = rate/2.0;
        return poisson_rv(r) + poisson_rv(r);
    }
    double L = exp(-1.0 * rate);
    double p = 1.0;
    double k = 0.0;    
    while (p >= L) {
        k += 1.0;
        p *= random_uniform();
    }
    return k - 1.0;
}

///////////////////////////////////////////////////////////////////////////////
// POPULATION ECOLOGY AND GENETICS

class Species;
class Population;
class Cell;
typedef std::vector<Cell> Cells;
typedef std::vector<Cell>::iterator CellIterator;
typedef std::vector<Species> SpeciesContainer;
typedef std::vector<Species>::iterator SpeciesIterator;
typedef std::vector<Species>::const_iterator SpeciesConstIterator;
typedef std::vector<Population>::iterator VecPopIterator;


const unsigned genotypeLen = 10;
#if defined(STATIC_GENOTYPE_LENGTH)
	typedef float Genotype[genotypeLen];
#else
	typedef std::vector<float> Genotype;
#endif
/// A single organism of a population of a particular species.
/// Responsible for tracking (non-neutral) genotype and neutral marker 
/// histories. 
class Individual {
    public:
        /// default constructor
        Individual() {
            // TODO!
#			if !defined(STATIC_GENOTYPE_LENGTH)
            	genotype.resize(genotypeLen);
#			endif
        }            
                   
        /// instantiates individual with given genotype
#		if defined(STATIC_GENOTYPE_LENGTH)
			Individual(const Genotype& g) {
				memcpy(this->genotype, g, sizeof(Genotype));
			}            
#		else
			Individual(const Genotype& g)
				:genotype(g)
				{}
#		endif
                
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
        Population(const Species* sp=NULL, const Cell* c=NULL)
            : species(sp),
              cell(c) {
        }
        void setCell(const Cell * c) {
        	this->cell = c;
        }
        void add_individual(const Individual& individual) {
            this->individuals.push_back(individual);
        }
        void assign(const unsigned n, const Individual & individual) {
        	this->individuals.assign(n, individual);
        }
        void reserve(const unsigned n) {
        	this->individuals.reserve(n);
        }
        void resize(const unsigned n) {
        	this->individuals.resize(n);
        }
        int size() const {
            return this->individuals.size(); 
        }
    
    private:
        const Species*    species; // the species to which this population belongs
        const Cell*       cell;    // the current location of this population
        std::vector<Individual> individuals; // members of this population        
}; // Population

/// A collection of Populations sharing the same ecologies (e.g. movement, 
/// fitness/survival functions, breeding pool)
class Species {
    public:
        Species() {}
        Species(const char* sp_label) 
            : label(sp_label) {       
        }
        ~Species() {}
        Population get_population(Cell* cell=NULL, int mean_size=0, int max_size=0) const;
        void initialize_population(Population * pop, Cell* cell=NULL, int mean_size=0, int max_size=0) const;
        
    private:
        std::string label;
//         std::list<Population*> populations;
}; // Species

/// Returns a Population object with the Population object's Species pointer
/// set to self.
/// If Cell* is given, then the Population object's Cell pointer is set.
/// If mean_size > 0, then individuals are created and added to the population,
/// with the number of individuals >= 0 and <= max_size.
/// The exact number of individuals is currently implemented as a Poisson 
/// distributed random number with mean of size. 
/// Future implementations and/or derived classes will take into account the 
/// environment of cell, with favorable environments resulting in greater 
/// numbers of individuals.
Population Species::get_population(Cell* cell, int mean_size, int max_size) const {
    Population p = Population(this, cell);
    this->initialize_population(&p, cell, mean_size, max_size);
    return p;
}

void Species::initialize_population(Population* popPtr, Cell* cell, int mean_size, int max_size) const {
	if (popPtr == 0L)
		return;
    Population & p = *popPtr;
	p.setCell(cell);
    if (mean_size > 0) {
        int n = poisson_rv(mean_size);
        while ((max_size > 0) and (n > max_size)) {
            n = poisson_rv(mean_size);
        }
#		if 0
		//@Jeet, it is much more efficient to use assign, or some "range insertion" 
		//	rather than a loop with push_back.  If you do lots of push_backs it is
		//	usually a good idea to call reserve to make sure the vector is big 
		//	enough before you start pushing back
			//p.reserve(n);
			for (; n > 0; --n) {
				p.add_individual(Individual());
			}
#		else 
			p.assign(n, Individual());
			//p.resize(n);
#		endif
    }
}

///////////////////////////////////////////////////////////////////////////////
// SPATIAL FRAMEWORK

/// The fundamental atomic spatial unit of the world.
class Cell {
    public:
    
        Cell(const SpeciesContainer* sp)
            : species(sp) {}
            
        Cell(const Cell& c)
            : species(c.species),
              populations(c.populations) {
            this->carrying_capacity = c.carrying_capacity;
        }            
    
        void set_carrying_capacity(int cc) {
            this->carrying_capacity = cc;
        }
        
        void initialize_populations();
        
        // for debugging
        int num_individuals() {
            return this->populations[0].size();
        }
    
    private:
        int                        carrying_capacity;
        const SpeciesContainer*    species;
        std::vector<Population>    populations;    

}; // Cell

// creates slots for populations, corresponding to species
#if 0
	void Cell::initialize_populations() {
		for (SpeciesConstIterator sp=this->species->begin(); 
				sp != this->species->end();
				sp++) {
			this->populations.push_back(sp->get_population(this, this->carrying_capacity/2, this->carrying_capacity)); 
		}
	}
#elif 0
	void Cell::initialize_populations() {
		for (SpeciesConstIterator sp=this->species->begin(); 
				sp != this->species->end();
				sp++) {
	        this->populations.push_back(Population());
    	    Population & population = *this->populations.rbegin();
        	sp->initialize_population(&population, this, this->carrying_capacity/2, this->carrying_capacity);
    	}
    }
#else
	void Cell::initialize_populations() {
		this->populations.resize(this->species->size());
		SpeciesConstIterator spIt = this->species->begin();
		VecPopIterator pIt = this->populations.begin();
		for (; spIt != this->species->end(); ++pIt, ++spIt)
				spIt->initialize_population(&(*pIt), this, this->carrying_capacity/2, this->carrying_capacity);
	}
#endif

/// The world.
class World {

    public:
        World();
        World(int dim_x, int dim_y, const SpeciesContainer& spp);
        void generate_landscape(int dim_x, int dim_y);
        void set_cell_carrying_capacity(int carrying_capacity);
        void add_species(const Species& sp);
        void initialize_biota();
        
        // for debugging
        void dump(std::ostream& out) {
            int count=1;
            for (CellIterator i=this->cells.begin(); 
                i != this->cells.end(); ++i, ++count) {
                    out << count << ": " << i->num_individuals() << "\n";
            }            
        }
        
    private:
        int                 dim_x;
        int                 dim_y;
        std::vector<Cell>   cells;
        SpeciesContainer    species;   
        
}; // World

/// default constructor
World::World() {}

/// constructor: calls landscape initializer
World::World(int dim_x, int dim_y, const SpeciesContainer& spp)
    : species(spp) {   
    this->generate_landscape(dim_x, dim_y);
}

/// generates the spatial framework
void World::generate_landscape(int dim_x, int dim_y) {
    this->dim_x = dim_x;
    this->dim_y = dim_y;
    int num_cells = dim_x * dim_y;
    // this->cells.reserve(num_cells);
    for (int i=0; i < num_cells; ++i) {                    
        this->cells.push_back(Cell(&this->species)); // Cell objects created here, ownership = World.cells
    }
}

/// sets the carrying capacity for all cells
void World::set_cell_carrying_capacity(int carrying_capacity) {
    // std::for_each(this->cells.begin(), this->cells.end(), Cell::set_carrying_capacity());
    for (CellIterator i=this->cells.begin(); i != this->cells.end(); ++i) {
        i->set_carrying_capacity(carrying_capacity);
    }
}

/// adds a species to the World
void World::add_species(const Species& sp) {
    this->species.push_back(sp);
}

/// initializes biota over landscape
/// must be called after landscape has been generated and species
/// list populated
void World::initialize_biota() {
    for (CellIterator i=this->cells.begin(); i != this->cells.end(); ++i) {
        i->initialize_populations();
    }
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
    int cc = atoi(argv[3]);
    
    // build world
    World world;   	
	world.generate_landscape(dim_x, dim_y);
	world.set_cell_carrying_capacity(cc);
	world.add_species(Species("snail"));
	world.initialize_biota();
	world.dump(std::cout);
	return 0;
}


