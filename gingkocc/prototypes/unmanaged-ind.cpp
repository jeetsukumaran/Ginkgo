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

#if defined(NDEBUG)
    #define DEBUG_BLOCK(y) y;
#else
    #define DEBUG_BLOCK(y)
#endif    


/************************* SUPPORT CLASSES AND METHODS ***********************/

 
///////////////////////////////////////////////////////////////////////////////
//! Wraps random number generator seed etc.
class RandomNumberGenerator {

    public:
        RandomNumberGenerator();
        RandomNumberGenerator(unsigned int seed);
        void set_seed(unsigned int seed);

        float uniform();   // [0, 1)
        float randint(int a, int b);
        float standard_normal();
        float normal(float mean, float sd);
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
float RandomNumberGenerator::uniform() {
    return static_cast<float>(rand())/static_cast<float>(RAND_MAX);
}

//! returns a uniform random integer between >= a and <= b
float RandomNumberGenerator::randint(int a, int b) {
    return (rand() % (b-a+1)) + a;
}

//! Gaussian distribution with mean=1 and sd=0
//! from Knuth, The Art of Computer Programming, Sec 3.4.1, Algorithm P
float RandomNumberGenerator::standard_normal() {

    // since this method generates two variates at a time,
    // we store the second one to be returned on the next call
    static float stored_variate = 0;
    static bool return_stored=false;    

    if (return_stored) {
        return_stored = false;
        return stored_variate;
    }

    float u1;
    float u2;
    float v1;
    float v2;
    float s = 1; 
    
    while (s >= 1.0) {
        u1 = this->uniform();
        u2 = this->uniform();
        v1 = 2.0 * u1 - 1.0;
        v2 = 2.0 * u2 - 1.0;
        s = pow(v1, 2) + pow(v2, 2);
    }
    
    float polar = sqrt( (-2 * log(s)) / s );
    float x1 = v1 * polar;
    stored_variate = v2 * polar;
    
    return x1;
}

//! Gaussian with given mean and sd
float RandomNumberGenerator::normal(float mean, float sd) {
    return this->standard_normal() * sd + mean;
}

//! Poisson r.v. with given rate
unsigned int RandomNumberGenerator::poisson(int rate) {
    const int MAX_EXPECTATION = 64;
    if (rate > MAX_EXPECTATION) {
        float r = rate/2.0;
        return this->poisson(r) + this->poisson(r);
    }
    float L = exp(-1.0 * rate);
    float p = 1.0;
    float k = 0.0;    
    while (p >= L) {
        k += 1.0;
        p *= this->uniform();
    }
    return k - 1.0;
}


/*********************** POPULATION ECOLOGY AND GENETICS *********************/
/***********************        (DECLARATION)           **********************/

class Individual;
class Species;
class Population;
class Landscape;
class Cell;
class World;

typedef int EnvironmentFactor;
typedef std::vector<EnvironmentFactor> EnvironmentFactors;
typedef std::vector<Cell> Cells;
typedef std::vector<Species> SpeciesContainer;
typedef std::vector<Individual> Individuals;

const unsigned GENOTYPE_LEN = 10;

///////////////////////////////////////////////////////////////////////////////
//! A single organism of a population of a particular species.
//! Responsible for tracking (non-neutral) genotype and neutral marker 
//! histories. 
class Individual {
    public:
    
        typedef std::vector<float> Genotype;
    
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
        
        Individual(const Individual&, const Individual& female); // offspring
        Individual(const Population& population);                // de novo
        Individual(const Individual& ind);                       // copy
        const Individual& operator=(const Individual& ind);      // assignment
                   
        void set_population(const Population& pop);
        
        Genotype& get_genotype() {
            return this->genotype;
        }
        
        void set_fitness(float fitness) {
            this->fitness = fitness;
        }
        
        void clear_fitness() {
            this->fitness = -1;
        }                
        
        float get_fitness() const {
            return this->fitness;
        }
        
        bool is_male() const {
            return this->sex == Individual::Male;
        }
        
        bool is_female() const {
            return this->sex == Individual::Female;
        }                

    private:
        const Population*     population; // host population
        RandomNumberGenerator *rng; 
        int                   num_environment_factors;
        Genotype              genotype;   // non-neutral genotype: maps to fitness phenotype
        Individual::Sex       sex;        // male or female
        int                   movement_reserve;   // movement "currency" available (reset every round)
        float                 fitness;    // individual's fitness (calc. and cached every round)
        

}; // Individual

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
        const Cell& get_cell() const {
        	return *(this->cell);
        }
        void set_species(const Species& sp) {
        	this->species = &sp;
        }
        const Species& get_species() const {
        	return *(this->species);
        }        
        
        World& get_world();
        
        Individuals& get_individuals() {
            return this->individuals;
        }
        
        void assign(unsigned int n) {
            // need to set reference to this population to each individual
            this->individuals.clear();
            this->individuals.reserve(n);
            for (unsigned int i = 0; i < n; ++i) {            
                Individual ind(*this);
                this->individuals.push_back(ind);  
            }                
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
        void reserve(unsigned int n) {
            this->individuals.reserve(n);
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
        void partition_by_sex(std::vector<Individual*>& males,
            std::vector<Individual*>& females) {
            males.clear();
            females.clear();
//             unsigned int est_size = static_cast<unsigned int>(this->individuals.size()/2);
//             males.reserve(est_size);
//             females.reserve(est_size);
            for (Individuals::iterator ind = this->individuals.begin();
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
        Individuals individuals; // the individuals of this population
     
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
       
        void set_index(int index) {
            this->index = index;    
        }
        int get_index() const {
            return this->index;    
        }            
                       
        World& get_world() const;
        void set_world(World& world);
                       
        // life history framework
        void environmental_selection() const;
        void population_reproduction() const;     
        void population_migration() const;
        
        // fitness/survival/competition
        virtual float calc_fitness(Individual& individual, EnvironmentFactors& env) const;
        
        // species-specific models
        virtual Population& spawn_offspring(std::vector<Individual*>& male_ptrs,
                                            std::vector<Individual*>& female_ptrs,
                                            Population& offspring) const;

    protected:        
        static std::vector<Individual*> male_ptrs;
        static std::vector<Individual*> female_ptrs;
        static Population               offspring;          

    private:
        std::string                     label;                  // arbitrary identifier
        int                             index;                  // "slot" in cell's pop vector
        int                             movement_rate;          // modifier to the movement surface
        std::vector<float>              selection_strengths;    // weighted_distance = distance / (sel. strength)
        World*                          world;                  // pointer to world
        RandomNumberGenerator*          rng;                    // pointer to rng
        
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
        std::vector<Population>& get_populations() {
            return this->populations;
        }
        
        EnvironmentFactors& get_environment() {
            return this->environment;
        }
        
        void repopulate(const Species& species, const Population& population);
        void density_dependent_selection();
        
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
        int                         carrying_capacity;  // max # ind
        std::vector<Population>     populations;        // list of local pops
        EnvironmentFactors          environment;        // environmental factors
        Landscape*                  landscape;          // host landscape
        World*                      world;              // host world

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
        
        int get_num_environment_factors() {
            return this->num_environment_factors;
        }
        
        unsigned long get_current_generation() {
            return this->current_generation;
        }
        
        void density_dependent_selection();
        
        void run_cycles(unsigned long num_generations);
                        
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
        int                     num_environment_factors;
        unsigned long           current_generation;                
        
}; // World


/*********************** POPULATION ECOLOGY AND GENETICS *********************/
/***********************        (IMPLEMENTATION         **********************/

///////////////////////////////////////////////////////////////////////////////
// Individual (IMPLEMENTATION)
    
Individual::Individual(const Population& population) {
    this->set_population(population);
    this->sex = Individual::random_sex(*this->rng);
    this->genotype.assign(GENOTYPE_LEN, 0.0);
    this->movement_reserve = 0;
    this->fitness = -1;
}

Individual::Individual(const Individual& ind)
    : genotype(ind.genotype) {
    *this = ind; // genotype initialized twice?
}

Individual::Individual(const Individual& female, const Individual& male) {
    this->set_population(*female.population);
    
    // sex assignment
    this->sex = Individual::random_sex(*this->rng);
    
    // genotype inheritance
    assert(male.genotype.size() == female.genotype.size());
    this->genotype.reserve(female.genotype.size());
    Genotype::const_iterator male_g = male.genotype.begin();
    Genotype::const_iterator female_g = female.genotype.begin();
    for ( ; female_g != female.genotype.end(); ++male_g, ++female_g) {
        float val = this->rng->choice(*male_g, *female_g);        
        // here is where we mutate genotype: pretty simple for now
        // and some magic numbers
        if (this->rng->uniform() < 0.1) {
            val += (this->rng->uniform() - 0.5); // i.e., -0.5 to +0.5
        }         
        this->genotype.push_back(val);
    }
    this->movement_reserve = 0;
    this->fitness = -1;
}

const Individual& Individual::operator=(const Individual& ind) {
    this->population = ind.population;
    this->genotype = ind.genotype;
    this->sex = ind.sex;
    this->movement_reserve = ind.movement_reserve;
    this->fitness = ind.fitness;
    return *this;
}

void Individual::set_population(const Population& pop) {
    this->population = &pop;
    this->rng = &this->population->get_species().get_world().get_rng(); 
    this->num_environment_factors = this->population->get_species().get_world().get_num_environment_factors();
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

World& Species::get_world() const {
    return *(this->world);    
}        

void Species::set_world(World& world) {
    this->world = &world;
    this->rng = &world.get_rng();
} 

void Species::environmental_selection() const {

    for (Cells::iterator cell = this->world->get_cells().begin();
            cell != this->world->get_cells().end();
            ++cell) {
        Population& pop = cell->get_populations().at(this->index);
        EnvironmentFactors& env = cell->get_environment();
        Individuals& individuals = pop.get_individuals();
        Population survivors = Population(*this, *cell);
        survivors.reserve(individuals.size());
        
//##DEBUG##
DEBUG_BLOCK( std::cout << individuals.size() << " "; )
        
        
        for (Individuals::iterator i = individuals.begin();
                i != individuals.end();
                ++i) {                     
            float fitness = this->calc_fitness(*i, env);   
            if (this->rng->uniform() <= fitness) {
                survivors.add(*i);
            }         
        }                
        cell->repopulate(*this, survivors);           
    }    
    
//##DEBUG##
DEBUG_BLOCK( std::cout << std::endl; )    
    
} 

void Species::population_reproduction() const {          
    for (Cells::iterator cell = this->world->get_cells().begin();
            cell != this->world->get_cells().end();
            ++cell) {        
        // assuming these delete without actually resizing capacity
        this->male_ptrs.clear();
        this->female_ptrs.clear();
        cell->get_populations().at(this->index).partition_by_sex(this->male_ptrs, 
                this->female_ptrs);
        this->offspring.clear();
        this->offspring.set_cell(*cell);
        this->offspring.set_species(*this);
        cell->repopulate(*this, this->spawn_offspring(female_ptrs, male_ptrs, offspring));
    }
}    

void Species::population_migration() const {
} 


float Species::calc_fitness(Individual& individual, EnvironmentFactors& env) const {
    float fitness = 0;
    Individual::Genotype& gen = individual.get_genotype();
    if (env.at(0) < 0 or env.at(0) > 10) {
        fitness = 0;
    } else {
        fitness = 5; // magic number!
        EnvironmentFactors::const_iterator eiter = env.begin();
        Individual::Genotype::const_iterator giter = gen.begin();
        for ( ; eiter != env.end(); ++eiter, ++giter) {
            // TODO: 'Species' to have weighting coefficients or powers 
            // that will be added here
            fitness += *eiter * *giter;
        }
    }
    individual.set_fitness(fitness);
    return fitness;
}

//! Derived classes should override this to implement different mating
//! and reproduction models. 
Population& Species::spawn_offspring(std::vector<Individual*>& female_ptrs,
                           std::vector<Individual*>& male_ptrs,                           
                           Population& offspring) const {
    if (female_ptrs.size() == 0 or male_ptrs.size() == 0) {
        return offspring;
    }
    for (std::vector<Individual*>::iterator female = female_ptrs.begin();
            female != female_ptrs.end();
            ++female) {
        // sampling with replacement            
        Individual* male = this->world->get_rng().choice(male_ptrs);
        for (int i=0; i < 2; ++i) {
            // change to (*female)->add_offspring())
            // then migrate
            // then recreate new populations from females in-situ
            offspring.add(Individual(*(*female), *male));
        }            
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
    for (SpeciesContainer::const_iterator spIt = this->world->get_species().begin();
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

//! cross-species competition within each cell
void Cell::density_dependent_selection() { 
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
    for (Cells::iterator i=this->cells.begin(); i != this->cells.end(); ++i) {
        i->set_carrying_capacity(carrying_capacity);
    }
}
	
///////////////////////////////////////////////////////////////////////////////
// World (IMPLEMENTATION)
	
//! default constructor
World::World()
    : rng(time(0)) {
    this->current_generation = 0;
}

//! constructor: calls
World::World(unsigned int seed) 
    : rng(seed) { 
    this->current_generation = 0;    
}           

//! generates a new landscape
void World::generate_landscape(int dim_x, int dim_y) {
    this->landscape.generate(*this, dim_x, dim_y);
    this->cells = &landscape.get_cells();
    	
	// actual implementation will load from file
	this->num_environment_factors = 4;
    for (Cells::iterator cell = this->get_cells().begin();
        cell != this->get_cells().end();
        ++cell) {
        cell->get_environment().assign(this->num_environment_factors, 1);
    }      
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
    for (Cells::iterator cell=this->get_cells().begin(); 
            cell != this->get_cells().end(); 
            ++cell) {
        cell->initialize_biota();
    }
    this->current_generation = 1;
}

//! runs a single iteration of a lifecycle 
void World::cycle() {

//##DEBUG##
DEBUG_BLOCK( std::cout << "\n*** GENERATION " << this->get_current_generation() << " ***\n\n"; )

    // survival
    for (SpeciesContainer::iterator sp_iter=this->species.begin(); 
            sp_iter != this->species.end(); 
            ++sp_iter) {
        sp_iter->environmental_selection();
    }          
    
    // competition
//     for (SpeciesContainer::iterator sp_iter=this->species.begin(); 
//             sp_iter != this->species.end(); 
//             ++sp_iter) {
//         sp_iter->density_dependent_selection();
//     }       
//     
    // reproduction
    for (SpeciesContainer::iterator sp_iter=this->species.begin(); 
            sp_iter != this->species.end(); 
            ++sp_iter) {
        sp_iter->population_reproduction();
    }    
    
    // migration
    for (SpeciesContainer::iterator sp_iter=this->species.begin(); 
            sp_iter != this->species.end(); 
            ++sp_iter) {
        sp_iter->population_migration();
    }                  
}

void World::density_dependent_selection() {

}

// run multiple cycles
void World::run_cycles(unsigned long num_generations) {
    for ( ; this->current_generation < num_generations; 
            ++this->current_generation) {
        this->cycle();
    }        
}


// for debugging
void World::dump(std::ostream& out) {
    int count=1;
    long indCount = 0;
    long indCapacityCount = 0;
    for (Cells::iterator i=this->get_cells().begin(); 
            i != this->get_cells().end(); 
            ++i, ++count) {
        const unsigned ni =  i->num_individuals();
//         out << ni << " ";
        indCount += ni;
        indCapacityCount += i->ind_capacity();
    }
    out << "\nTotal individuals = " << indCount << '\n';
    out << "Total individual capacity = " << indCapacityCount << '\n';
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
	
	
	
	for (Cells::iterator cell = world.get_cells().begin();
	        cell != world.get_cells().end();
	        ++cell) {
	    for (SpeciesContainer::iterator sp = world.get_species().begin();
	            sp != world.get_species().end();
	            ++sp) {
	            cell->seed_population(*sp, cc);
        }	            
    }	
	
	world.run_cycles(num_gens);
	
	world.dump(std::cout);
	return 0;
}


