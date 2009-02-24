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
#include <ctime>
#include <cassert>

/************************* SUPPORT CLASSES AND METHODS ***********************/
 
///////////////////////////////////////////////////////////////////////////////
//! Wraps random number generator seed etc.
class RandomNumberGenerator {

    public:
        RandomNumberGenerator();
        RandomNumberGenerator(unsigned int seed);
        void set_seed(unsigned int seed);

        double uniform();   // [0, 1)
        double randint(int a, int b);
        double standard_normal();
        double normal(double mean, double sd);
        unsigned int poisson(int rate);
        
        template <typename T>
        inline T& choice(std::vector<T>& v) {
            return v[rand() % v.size()];
        }
        
        template <typename T>        
        inline T& choice(T& a, T&b) {
            if (this->uniform() < 0.5) {
                return a;
            } else {
                return b;
            }
        }        
        
    private:
        unsigned int seed;

};

//! seeds using time
RandomNumberGenerator::RandomNumberGenerator() {
    this->set_seed(time(0));
}

//! seeds using given seed
RandomNumberGenerator::RandomNumberGenerator(unsigned int seed) {
    this->set_seed(seed);
}

//! seeds using given seed
void RandomNumberGenerator::set_seed(unsigned int seed) {
    this->seed = seed;
    srand(seed);
}

//! returns a uniform random real between 0 and 1
double RandomNumberGenerator::uniform() {
    return static_cast<double>(rand())/static_cast<double>(RAND_MAX);
}

//! returns a uniform random integer between >= a and <= b
double RandomNumberGenerator::randint(int a, int b) {
    return (rand() % (b-a+1)) + a;
}

//! Gaussian distribution with mean=1 and sd=0
//! from Knuth, The Art of Computer Programming, Sec 3.4.1, Algorithm P
double RandomNumberGenerator::standard_normal() {

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
        u1 = this->uniform();
        u2 = this->uniform();
        v1 = 2.0 * u1 - 1.0;
        v2 = 2.0 * u2 - 1.0;
        s = pow(v1, 2) + pow(v2, 2);
    }
    
    double polar = sqrt( (-2 * log(s)) / s );
    double x1 = v1 * polar;
    stored_variate = v2 * polar;
    
    return x1;
}

//! Gaussian with given mean and sd
double RandomNumberGenerator::normal(double mean, double sd) {
    return this->standard_normal() * sd + mean;
}

//! Poisson r.v. with given rate
unsigned int RandomNumberGenerator::poisson(int rate) {
    const int MAX_EXPECTATION = 64;
    if (rate > MAX_EXPECTATION) {
        double r = rate/2.0;
        return this->poisson(r) + this->poisson(r);
    }
    double L = exp(-1.0 * rate);
    double p = 1.0;
    double k = 0.0;    
    while (p >= L) {
        k += 1.0;
        p *= this->uniform();
    }
    return k - 1.0;
}


/*********************** POPULATION ECOLOGY AND GENETICS *********************/
/***********************        (DECLARATION)           **********************/
 
class Species;
class Population;
class Landscape;
class Cell;
class World;

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

///////////////////////////////////////////////////////////////////////////////
//! A single organism of a population of a particular species.
//! Responsible for tracking (non-neutral) genotype and neutral marker 
//! histories. 
class Individual {
    public:
    
        enum Sex {
            Male,
            Female
        };
        
        static Individual::Sex random_sex(RandomNumberGenerator& rng, 
                float female_threshold=0.5) {
            if (rng.uniform() < female_threshold) {
                return Individual::Male;
            } else {
                return Individual::Female;
            }
        }
        
        // construct a new ind from two existing ones1
        Individual(const Individual&, const Individual& female);
    
        // construct an individual de novo
        Individual(Population& population) {
            this->population = &population;
            this->sex = Individual::random_sex(this->get_rng());
            
#if !defined(STATIC_GENOTYPE_LENGTH)
            this->genotype.assign(genotypeLen, 0.0);
#else
            for (float* g=this->genotype; g < STATIC_GENOTYPE_LEN; ++g) {
                *g = 0.0;
            }
#endif
        }
        
        // clones an individual
        Individual(const Individual& ind)
            : genotype(ind.genotype) {
            this->population = ind.population;
            this->sex = ind.sex;
        }
                   
        void set_population(Population& pop) {
            this->population = &pop;
        }
        
        bool is_male() const {
            return this->sex == Individual::Male;
        }
        
        bool is_female() const {
            return this->sex == Individual::Female;
        }                
        
        RandomNumberGenerator& get_rng();

// #		if defined(STATIC_GENOTYPE_LENGTH)
// 			Individual(const Genotype& g) {
// 				memcpy(this->genotype, g, sizeof(Genotype));
// 			}            
// #		else
// 			Individual(const Genotype& g)
// 				:genotype(g)
// 				{}
// #		endif
                
    private:
        Population*     population; // host population
        Genotype        genotype;   // non-neutral genotype: maps to fitness phenotype
        Individual::Sex sex;
        

}; // Individual

#if defined(STATIC_GENOTYPE_LENGTH)
	const unsigned SIZE_OF_INDIVIDUAL = sizeof(Individual);
#else
	const unsigned SIZE_OF_INDIVIDUAL = sizeof(Individual) +  genotypeLen*sizeof(float) ; // we have to take into account the size of the vector's internal storage
#endif

///////////////////////////////////////////////////////////////////////////////
//! A single population of a particular species.
//! Responsible for managing collections of individuals, and relating them to 
//! their species.
class Population {
    public:
    
        // special case: for the static 
        Population() {        
        }
    
        // lifecycle
        Population(const Species& sp, const Cell& c) {
            this->species = &sp;
            this->cell = &c;
        }       
        
        // accessors
        void set_cell(const Cell& cell) {
        	this->cell = &cell;
        }
        const Cell& get_cell() {
        	return *(this->cell);
        }
        void set_species(const Species& sp) {
        	this->species = &sp;
        }
        const Species& get_species() {
        	return *(this->species);
        }        
        World& get_world();
        
        void assign(unsigned int n) {
            // need to set reference to this population to each individual,
            Individual ind(*this);
            this->individuals.assign(n, ind);                        
        }        
        void assign(unsigned int n, const Individual& ind) {
            // need to set reference to this population to each individual,
            Individual ind_copy(ind);
            ind_copy.set_population(*this);
            this->individuals.assign(n, ind);                        
        }
        void add(const Individual& ind) {
            this->individuals.push_back(ind);
        }
        void clear() {
            this->individuals.clear();
        }
        unsigned int capacity() {
            return this->individuals.capacity();
        }        
        unsigned int size() {
            return this->individuals.size();
        }
        void partition_by_gender(std::vector<Individual*>& males,
            std::vector<Individual*>& females) {
            males.clear();
            females.clear();
//             unsigned int est_size = static_cast<unsigned int>(this->individuals.size()/2);
//             males.reserve(est_size);
//             females.reserve(est_size);
            for (std::vector<Individual>::iterator ind = this->individuals.begin();
                    ind != this->individuals.end();
                    ++ind) {
                if (ind->is_male()) {
                    males.push_back(&(*ind));
                } else {
                    females.push_back(&(*ind));
                }
            }
        }
                                     
    private:
        const Species*          species; // the species to which this population belongs
        const Cell*             cell;    // the current location of this population
        std::vector<Individual> individuals; // the individuals of this population
     
}; // Population


///////////////////////////////////////////////////////////////////////////////
//! A collection of Populations sharing the same ecologies (e.g. movement, 
//! fitness/survival functions, breeding pool)
class Species {
    public:
        // lifecycle
        Species();
        Species(const char* sp_label);
        virtual ~Species() {}
        
        // accessors
        World& get_world() const {
            return *(this->world);    
        }        
        void set_world(World& world) {
            this->world = &world;    
        }        
        void set_index(int index) {
            this->index = index;    
        }
        int get_index() const {
            return this->index;    
        }            
                       
        // operations
        virtual void reproduce() const; // reproduce all 
        virtual Population& breed(std::vector<Individual*>& male_ptrs,
                                  std::vector<Individual*>& female_ptrs,
                                  Population& offspring) const;

    private:
        std::string     label;
        int             index;
        World*          world;
        static std::vector<Individual*> male_ptrs;
        static std::vector<Individual*> female_ptrs;
        static Population               offspring;          

}; // Species

std::vector<Individual*> Species::male_ptrs;
std::vector<Individual*> Species::female_ptrs;
Population               Species::offspring;

/********************* SPATIAL AND ENVIRONMENTAL FRAMWORK ********************/
/***********************        (DECLARATION)           **********************/

///////////////////////////////////////////////////////////////////////////////
//! The fundamental atomic spatial unit of the world.
class Cell {
    public:
    
        // lifecycle
        Cell(Landscape& landscape);
        Cell(const Cell& cell);
        
        // operators
        const Cell& operator=(const Cell& cell);

        // accessors
        void set_carrying_capacity(int cc) {
            this->carrying_capacity = cc;
        }
        void set_landscape(Landscape& landscape);
        void repopulate(const Species& species, const Population& population);
        
        // operations
        void initialize_biota();
        void seed_population(Species& sp, unsigned int size);
        void reproduce_populations();
        
        // for debugging
        int num_individuals() {
            return this->populations.at(0).size();            
        }
        int ind_capacity() {
            return this->populations.at(0).size();
        }
    
    private:
        Landscape*                 landscape;
        World*                     world;
        int                        carrying_capacity;
        std::vector<Population>    populations;    

}; // Cell

///////////////////////////////////////////////////////////////////////////////	
//! The landscape.
class Landscape {

    public:
        Landscape();
        Landscape(World& world, int dim_x, int dim_y);
        Landscape(const Landscape& landscape);
        const Landscape& operator=(const Landscape& landscape);
        
        Cells& get_cells() {
            return this->cells;
        }
        
        World& get_world() {
            return *(this->world);
        }               
        
        void set_cell_carrying_capacity(int carrying_capacity);
        
        void generate(World& world, int dim_x, int dim_y);
        
    private:
        int      dim_x;
        int      dim_y;            
        Cells    cells;
        World*   world;

};

///////////////////////////////////////////////////////////////////////////////	
//! The world.
class World {

    public:

        // lifecycle
        World();
        World(unsigned int seed);    

        // accessors
        SpeciesContainer& get_species() { 
            return this->species;
        }
        
        RandomNumberGenerator& get_rng() { 
            return this->rng;
        }
        
        
        Landscape& get_landscape() {
            return this->landscape;
        }        
        
        Cells& get_cells() {
            return this->landscape.get_cells();
        }
        
        // set up
        void generate_landscape(int dim_x, int dim_y);
        void set_cell_carrying_capacity(int carrying_capacity);
        void add_species(const Species& sp);
        void initialize_biota();
        
        // run
        void cycle();
        
        // debug
        void dump(std::ostream& out);
        
    private:
        Cells*                  cells;
        Landscape               landscape;
        SpeciesContainer        species;
        RandomNumberGenerator   rng;
        
}; // World


/*********************** POPULATION ECOLOGY AND GENETICS *********************/
/***********************        (IMPLEMENTATION         **********************/

///////////////////////////////////////////////////////////////////////////////
// Individual (IMPLEMENTATION)

RandomNumberGenerator& Individual::get_rng() {
    return this->population->get_species().get_world().get_rng();
}

Individual::Individual(const Individual& male, const Individual& female) {
    this->population = female.population;
    this->population->get_world();
    this->sex = Individual::random_sex(this->get_rng());

    assert(male.genotype.size() == female.genotype.size());
    this->genotype.reserve(female.genotype.size());
    Genotype::const_iterator male_g = male.genotype.begin();
    Genotype::const_iterator female_g = female.genotype.begin();
    for ( ; female_g != female.genotype.end(); ++male_g, ++female_g) {
        float val = this->get_rng().choice(*male_g, *female_g);
        
        // here is where we mutate genotype
        // pretty simple for now
        if (this->get_rng().uniform() < 0.1) {
            val += (get_rng().uniform() - 0.5); // i.e., -0.5 to +0.5
        } 
        
        this->genotype.push_back(val);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Population (IMPLEMENTATION)

World& Population::get_world() {
    return this->species->get_world();
}  

///////////////////////////////////////////////////////////////////////////////
// Species (IMPLEMENTATION)

//! default constructor: assigns dummy values
Species::Species() {
    this->label = "Sp";
    this->index = -1;
}

//! apart from setting label, ensures index is unassigned
Species::Species(const char* sp_label) 
    : label(sp_label) {
    this->index = -1;
}

//! Derived classes should override this to implement different reproduction
//! models. 
void Species::reproduce() const {          
    for (CellIterator cell = this->world->get_cells().begin();
            cell != this->world->get_cells().end();
            ++cell) {        
        // assuming these delete without actually resizing capacity
        this->male_ptrs.clear();
        this->female_ptrs.clear();
        this->offspring.clear();
        this->offspring.set_cell(*cell);
        this->offspring.set_species(*this);
        cell->repopulate(*this, this->breed(male_ptrs, female_ptrs, offspring));
    }
}    

//! Derived classes should override this to implement different reproduction
//! models. 
Population& Species::breed(std::vector<Individual*>& male_ptrs,
                           std::vector<Individual*>& female_ptrs,
                           Population& offspring) const {
    if (female_ptrs.size() == 0 or male_ptrs.size() == 0) {
        return offspring;
    }
    for (std::vector<Individual*>::iterator female = female_ptrs.begin();
            female != female_ptrs.end();
            ++female) {
        // sampling with replacement            
        Individual* male = this->world->get_rng().choice(male_ptrs);
        offspring.add(Individual(*male, **female));
    }
    return offspring;
}

/********************* SPATIAL AND ENVIRONMENTAL FRAMWORK ********************/
/***********************        (IMPLEMENTATION)        **********************/

///////////////////////////////////////////////////////////////////////////////
// Cell (IMPLEMENTATION)

//! constructor: needs reference to World
Cell::Cell(Landscape& landscape){
    this->set_landscape(landscape);
    this->carrying_capacity = 0;
}

//! copy constructor
Cell::Cell(const Cell& cell) {
    this->set_landscape(*cell.landscape);
    this->populations = cell.populations;
}

//! assignment
const Cell& Cell::operator=(const Cell& cell) {   
    this->set_landscape(*cell.landscape);
    this->carrying_capacity = cell.carrying_capacity;
    this->populations = cell.populations;
    return *this;
}

//! binds to landscape
void Cell::set_landscape(Landscape& landscape) {
    this->landscape = &landscape;
    this->world = &this->landscape->get_world();        
}

//! creates slots for populations, corresponding to #'s of species
void Cell::initialize_biota() {
    this->populations.reserve(this->world->get_species().size());
    for (SpeciesConstIterator spIt = this->world->get_species().begin();
            spIt != this->world->get_species().end();
            ++spIt) {
        this->populations.push_back(Population(*spIt, *this));
    }        
}

//! Adds new individuals to the population of the specified species in this
//! cell.
void Cell::seed_population(Species& sp, unsigned int size) {
    assert(sp.get_index() >= 0);
    assert(sp.get_index() <= this->populations.size());
    this->populations.at(sp.get_index()).assign(size);    
}

//! redefines a population (typically with the next gener
void Cell::repopulate(const Species& sp, const Population& population) {
    this->populations.at(sp.get_index()) = population;      
}

///////////////////////////////////////////////////////////////////////////////
// Landscape (IMPLEMENTATION)

Landscape::Landscape() {

}

Landscape::Landscape(World& world, int dim_x, int dim_y) {
    this->generate(world, dim_x, dim_y);
}

void Landscape::generate(World& world, int dim_x, int dim_y) {
    this->world = &world;
    this->dim_x = dim_x;
    this->dim_y = dim_y;
    int num_cells = dim_x * dim_y;
    this->cells.assign(num_cells, Cell(*this)); // Cell objects created here
}

//! sets the carrying capacity for all cells
void Landscape::set_cell_carrying_capacity(int carrying_capacity) {
    // std::for_each(this->cells.begin(), this->cells.end(), Cell::set_carrying_capacity());
    for (CellIterator i=this->cells.begin(); i != this->cells.end(); ++i) {
        i->set_carrying_capacity(carrying_capacity);
    }
}
	
///////////////////////////////////////////////////////////////////////////////
// World (IMPLEMENTATION)
	
//! default constructor
World::World()
    : rng(time(0)) {
}

//! constructor: calls
World::World(unsigned int seed) 
    : rng(seed) {  
}           

//! generates a new landscape
void World::generate_landscape(int dim_x, int dim_y) {
    this->landscape.generate(*this, dim_x, dim_y);
    this->cells = &landscape.get_cells();
}

//! sets the carrying capacity for each cell
void World::set_cell_carrying_capacity(int carrying_capacity) {
    this->landscape.set_cell_carrying_capacity(carrying_capacity);
}

//! adds a species to the World
void World::add_species(const Species& sp) {
    this->species.push_back(sp);
    Species& world_sp = this->species.back();
    world_sp.set_world(*this);
    world_sp.set_index(this->species.size() - 1);
}

//! initializes biota over landscape
//! must be called after landscape has been generated and species
//! list populated
void World::initialize_biota() {
    for (CellIterator i=this->get_cells().begin(); i != this->get_cells().end(); ++i) {
        i->initialize_biota();
    }
}

//! runs a single iteration of a lifecycle 
void World::cycle() {
    // survival
    // competition
    // reproduction
    for (SpeciesIterator sp_iter=this->species.begin(); 
            sp_iter != this->species.end(); 
            ++sp_iter) {
        sp_iter->reproduce();
    }    
    // migration
}

// for debugging
void World::dump(std::ostream& out) {
    int count=1;
    long indCount = 0;
    long indCapacityCount = 0;
    for (CellIterator i=this->get_cells().begin(); 
            i != this->get_cells().end(); 
            ++i, ++count) {
        const unsigned ni =  i->num_individuals();
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
	
	for (CellIterator cell = world.get_cells().begin();
	        cell != world.get_cells().end();
	        ++cell) {
	    for (SpeciesIterator sp = world.get_species().begin();
	            sp != world.get_species().end();
	            ++sp) {
	            cell->seed_population(*sp, cc);
        }	            
    }	
	
	for (int i=1; i <= num_gens; ++i) {
	    world.cycle();
	}
	
	world.dump(std::cout);
	return 0;
}


