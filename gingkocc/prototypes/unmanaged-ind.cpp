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
// SUPPORT FUNCTIONS: WILL BE PLACED IN CLASS RandomVariate

double uniform_variate() {
    return double(rand())/double(RAND_MAX);
}

// from Knuth, The Art of Computer Programming, Sec 3.4.1, Algorithm P
double standard_normal_variate() {

    // since this method generates two variates at a time,
    // we store the second one to be returned on the next call
    static double stored_variate = 0;
    static bool return_stored=false;    

    if (return_stored) {
        return_stored = false;
        return stored_variate;
    }

    double u1;
    double u2;
    double v1;
    double v2;
    double s = 1; 
    
    while (s >= 1.0) {
        u1 = uniform_variate();
        u2 = uniform_variate();
        v1 = 2.0 * u1 - 1.0;
        v2 = 2.0 * u2 - 1.0;
        s = pow(v1, 2) + pow(v2, 2);
    }
    
    double polar = sqrt( (-2 * log(s)) / s );
    double x1 = v1 * polar;
    stored_variate = v2 * polar;
    
    return x1;
}

double normal_variate(double mean, double sd) {
    return standard_normal_variate() * sd + mean;
}

int poisson_variate(int rate) {
    const int MAX_EXPECTATION = 64;
    if (rate > MAX_EXPECTATION) {
        double r = rate/2.0;
        return poisson_variate(r) + poisson_variate(r);
    }
    double L = exp(-1.0 * rate);
    double p = 1.0;
    double k = 0.0;    
    while (p >= L) {
        k += 1.0;
        p *= uniform_variate();
    }
    return k - 1.0;
}

template <typename T>
T random_sample(const std::vector<T>& v) {
    return v[rand() % v.size()];
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

#if defined(STATIC_GENOTYPE_LENGTH)
	const unsigned SIZE_OF_INDIVIDUAL = sizeof(Individual);
#else
	const unsigned SIZE_OF_INDIVIDUAL = sizeof(Individual) +  genotypeLen*sizeof(float) ; // we have to take into account the size of the vector's internal storage
#endif

/// A single population of a particular species.
/// Responsible for managing collections of individuals, and relating them to 
/// their species.
class Population : public std::vector<Individual> {
    public:
        Population(const Species* sp=NULL, const Cell* c=NULL)
            : species(sp),
              cell(c) {
        }
        void setCell(const Cell * c) {
        	this->cell = c;
        }

    private:
        const Species*    species; // the species to which this population belongs
        const Cell*       cell;    // the current location of this population
     
}; // Population

/// A collection of Populations sharing the same ecologies (e.g. movement, 
/// fitness/survival functions, breeding pool)
class Species {
    public:
        Species() {}
        Species(const char* sp_label) 
            : label(sp_label) {       
        }
        virtual ~Species() {}
        virtual Population get_population(Cell* cell=NULL, int mean_size=0, int max_size=0) const;
        virtual Population& reproduce(Population& cur_gen) const;
        virtual void initialize_population(Population* popPtr, Cell* cell, int mean_size, int max_size) const;
        
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
        int n = poisson_variate(mean_size);
        while ((max_size > 0) and (n > max_size)) {
            n = poisson_variate(mean_size);
        }
        p.assign(n, Individual());
    }
}

/// Returns next generation.
/// Derived classes should override this to implement different reproduction
/// models. 
/// Current model: next gen population size is a Poisson distributed random
/// variate with mean equal to current population size. 
//  Each offspring in the next generation randomly selects two parents 
//  from the current generation.
/// Of course, given no pop gen component right now, the latter does not hold.
Population& Species::reproduce(Population& cur_gen) const {
    static Population next_gen;
    // we assume that each individual produces a Poisson distributed number
    // of offspring with mean of 2; as there are cur_gen.size() individuals
    // the total number of offspring is cur_gen.size() * Poisson(2); the sum
    // of n Poisson variables with mean M = Poisson(n*M).
    // seems to only work when num_gens < 7 ... (at least, with no carrying
    // capacity enforced.
//     unsigned int next_gen_size = poisson_variate(cur_gen.size() * 2);
//     next_gen.assign(next_gen_size, Individual());

    // tweak: next gen population size is normally distributed with mean =
    // current pop size, and sd = 10% of current pop size
    unsigned int next_gen_size = normal_variate(cur_gen.size(), cur_gen.size()/10);
    next_gen.assign(next_gen_size, Individual());

    cur_gen = next_gen; // copy vals
    return cur_gen; 
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
        void reproduce_populations();
        
        // for debugging
        int num_individuals() {
            return this->populations[0].size();
        }
        int ind_capacity() {
            return this->populations[0].capacity();
        }
    
    private:
        int                        carrying_capacity;
        const SpeciesContainer*    species;
        std::vector<Population>    populations;    

}; // Cell

/// creates slots for populations, corresponding to species
void Cell::initialize_populations() {
    this->populations.resize(this->species->size());
    SpeciesConstIterator spIt = this->species->begin();
    VecPopIterator pIt = this->populations.begin();
    for (; spIt != this->species->end(); ++pIt, ++spIt)
        spIt->initialize_population(&(*pIt), this, this->carrying_capacity/2, this->carrying_capacity);
}

/// gets the next generation for each population
void Cell::reproduce_populations() {
    this->populations.resize(this->species->size());
    SpeciesConstIterator spIt = this->species->begin();
    VecPopIterator pIt = this->populations.begin();
    for (; spIt != this->species->end(); ++pIt, ++spIt)
        spIt->reproduce(*pIt);
}
	
/// The world.
class World {

    public:
        World();
        World(int dim_x, int dim_y, const SpeciesContainer& spp);
        void generate_landscape(int dim_x, int dim_y);
        void set_cell_carrying_capacity(int carrying_capacity);
        void add_species(const Species& sp);
        void initialize_biota();
        void cycle();
        
        // for debugging
        void dump(std::ostream& out) {
            int count=1;
            long indCount = 0;
            long indCapacityCount = 0;
            for (CellIterator i=this->cells.begin(); 
                i != this->cells.end(); ++i, ++count) {
                	const unsigned ni = i->num_individuals();
                    out << count << ": " << ni << "\n";
                    indCount += ni;
                    indCapacityCount += i->ind_capacity();
            }
            out << "Total individuals = " << indCount << '\n';
            out << "Total individual capacity = " << indCapacityCount << '\n';
            out << "Size of individual in bytes = " << SIZE_OF_INDIVIDUAL << std::endl;
            out << "Min. possible memory used for individuals = " << indCount*SIZE_OF_INDIVIDUAL << std::endl;
            out << "Memory used for individuals = " << indCapacityCount*SIZE_OF_INDIVIDUAL << std::endl;
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
    this->cells.reserve(num_cells);
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

/// runs a single iteration of a lifecycle 
void World::cycle() {
    // survival
    // competition
    // reproduction
    for (CellIterator cell_iter=this->cells.begin(); 
            cell_iter != this->cells.end(); 
            ++cell_iter) {
        cell_iter->reproduce_populations();
    }    
    // migration
}

///////////////////////////////////////////////////////////////////////////////
// TEST

int main(int argc, char * argv[]) {
    if (argc < 5) {
        std::cout << "usage: " << argv[0] <<  " <DIM-X> <DIM-Y> <CELL-CARRYING-CAPACITY> <NUM-GENS>\n";
        exit(1);
    }
    int dim_x = atoi(argv[1]);
    int dim_y = atoi(argv[2]);
    int cc = atoi(argv[3]);
    int num_gens = atoi(argv[4]);
    
    // build world
    World world;   	
	world.generate_landscape(dim_x, dim_y);
	world.set_cell_carrying_capacity(cc);
	world.add_species(Species("snail"));
	world.initialize_biota();
	for (int i=1; i <= num_gens; ++i) {
	    world.cycle();
	}
	world.dump(std::cout);
	return 0;
}


