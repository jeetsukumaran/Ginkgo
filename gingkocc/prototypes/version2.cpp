///////////////////////////////////////////////////////////////////////////////
//
// version2.cpp
//
// Performance evaluation of population evolution under non-managed allocation
// of organisms, implementing:
//      - conflated fecundity and (density-dependent) survival
//      - organisms managed as multi-species vectors in each cell
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
#include <set>
#include <cstring>
#include <list>
#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
#include <cassert>
#include <iomanip>

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
        RandomNumberGenerator(unsigned long seed);
        void set_seed(unsigned long seed);

        float random();   // [0, 1)
        float randint(int a, int b);
        float standard_normal();
        float normal(float mean, float sd);
        unsigned int poisson(int rate);
        
        template <typename T>
        inline typename T::value_type& choice(T& collection) {
            return collection[rand() % collection.size()];
        }
        
        template <typename T>        
        inline T& choice(T& a, T& b) {
            if (this->random() < 0.5) {
                return a;
            } else {
                return b;
            }
        }        
        
    private:
        unsigned long _seed;

};

//! seeds using time
RandomNumberGenerator::RandomNumberGenerator() {
    this->set_seed(time(0));
}

//! seeds using given seed
RandomNumberGenerator::RandomNumberGenerator(unsigned long seed) {
    this->set_seed(seed);
}

//! seeds using given seed
void RandomNumberGenerator::set_seed(unsigned long seed) {
    this->_seed = seed;
    srand(this->_seed);
}

//! returns a uniform random real between 0 and 1
float RandomNumberGenerator::random() {
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
        u1 = this->random();
        u2 = this->random();
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
        p *= this->random();
    }
    return k - 1.0;
}

/******************************************************************************
 * FORWARD DECLARATIONS AND GLOBAL TYPEDEFS
 *****************************************************************************/

class Organism;
class Species;
class Community;
typedef int GenotypeFactor;
typedef std::vector<GenotypeFactor> GenotypeFactors;
typedef std::vector<Species*> SpeciesPool;
typedef std::vector<Organism> Organisms;
typedef std::vector<Community> Communities;

class Landscape;
class Cell;
class World;
typedef long EnvironmentalFactor;
typedef std::vector<EnvironmentalFactor> EnvironmentalFactors;
typedef std::vector<Cell> Cells;


/******************************************************************************
 * DECLARATIONS
 *****************************************************************************/
 
/* POPULATION ECOLOGY AND GENETICS *******************************************/

///////////////////////////////////////////////////////////////////////////////
//! A single organism of a population of a particular species.
//! Responsible for tracking (non-neutral) genotype and neutral marker 
//! histories. Very lightweight, with most functionality delegated to other
//! classes.
class Organism {
    public:
    
        // gender
        enum Sex {
            Male,
            Female
        };

        // lifecycle and assignment
        
        Organism(int species_index,
                           const GenotypeFactors& new_genotype,
                           Organism::Sex new_sex) 
            : _species_index(species_index),
              _genotype(new_genotype),
              _sex(new_sex),
              _fitness() {
        }
        
        //! Copy constructor.
        Organism(const Organism& ind) {
            *this = ind;
        }
        
        //! Assignment.
        const Organism& operator=(const Organism& ind) {
            this->_species_index = ind._species_index;
            this->_genotype = ind._genotype;
            this->_sex = ind._sex;
            this->_fitness = ind._fitness;
            return *this;
        }
                   
        // genotype       
        GenotypeFactors& genotype() {
            return this->_genotype;
        }
        
        // fitness & survival
        float &fitness() {
            return this->_fitness;
        }
        bool &killed() {
            return this->_killed;
        }        
        
        // meta-info
        int get_species_index() const {
            return this->_species_index;
        }
        
        bool is_male() const {
            return this->_sex == Organism::Male;
        }
        
        bool is_female() const {
            return this->_sex == Organism::Female;
        }                

    private:
        int              _species_index;    // species
        GenotypeFactors  _genotype;         // non-neutral genotype: maps to fitness phenotype
        Organism::Sex    _sex;              // male or female
        float            _fitness;          // cache this organism's fitness
        bool             _killed;           // flag an organism as dead: allowing for use of std::remove_if() and std::resize() or v.erase()
        
}; // Organism

///////////////////////////////////////////////////////////////////////////////
//! A collection of processes and properties that determine the ecologies of
//! organisms.
class Species {

    public:
    
        // --- lifecycle and assignment ---        
        Species(int index,
                const char* label, 
                int num_fitness_factors,
                RandomNumberGenerator& rng);
        Species(const Species& species);                
        ~Species() {}        
                
        // --- access and mutation ---
        int get_index() const {
            return this->_index;
        }
        void set_num_fitness_factors(int i) {
            this->_num_fitness_factors = i;
        }
        int get_num_fitness_factors() const {
            return this->_num_fitness_factors;
        }
        void set_index(int i) {
            this->_index = i;
        }        
        int get_movement_rate() const {
            return this->_movement_rate;
        }
        void set_movement_rate(int i) {
            this->_movement_rate = i;
        }        
        std::vector<float>& selection_strengths() {
            return this->_selection_strengths;
        }
        std::vector<unsigned>& movement_surface() {
            return this->_movement_surface;
        }
        float get_mutation_rate() const {
            return this->_movement_rate;
        }
        void set_mutation_rate(float i) {
            this->_movement_rate = i;
        }
        GenotypeFactor get_max_mutation_size() const {
            return this->_max_mutation_size;
        }
        void set_max_mutation_size(GenotypeFactor i) {
            this->_max_mutation_size = i;
        }
        int get_mean_reproductive_rate() const {
            return this->_mean_reproductive_rate;
        }
        void set_mean_reproductive_rate(int i) {
            this->_mean_reproductive_rate = i;
        }
        int get_reproductive_rate_mutation_size() const {
            return this->_reproductive_rate_mutation_size;
        }
        void set_reproductive_rate_mutation_size(int i) {
            this->_reproductive_rate_mutation_size = i;
        }
        GenotypeFactors& default_genotype() {
            return this->_default_genotype;
        }
        
        // --- organism generation and reproduction ---
        Organism::Sex get_random_sex(float female_threshold=0.5) const {
            if (this->_rng.random() < female_threshold) {
                return Organism::Male;
            } else {
                return Organism::Female;
            }
        }  
        
        Organism new_organism() const {
            return Organism(this->_index, this->_default_genotype, this->get_random_sex());
        }
                                
    private:
        const Species& operator=(const Species& species);

        int                         _index;                     // "slot" in cell's pop vector    
        std::string                 _label;                     // arbitrary identifier
        int                         _num_fitness_factors;       // so genotypes of appropriate length can be composed                
        int                         _movement_rate;             // modifier to the global movement surface to derive the species-specific movement surface
        std::vector<float>          _selection_strengths;       // weighted_distance = distance / (sel. strength)
        std::vector<unsigned>       _movement_surface;          // the movement surface: the "cost" to enter into every cell on the landscape
        float                       _mutation_rate;             // rate of mutations
        GenotypeFactor              _max_mutation_size;         // window "size" of mutations
        int                         _mean_reproductive_rate;    // "base" reproductive rate
        int                         _reproductive_rate_mutation_size;  // if reprod. rate evolves, size of step
        GenotypeFactors             _default_genotype;          // genotype of individuals generated de novo        
        RandomNumberGenerator&      _rng;                       // rng to use

}; // Species

/* LANDSCAPE AND SPATIAL RELATIONS *******************************************/

///////////////////////////////////////////////////////////////////////////////
//! The fundamental atomic spatial unit of the world.
class Cell {
    public:
    
        // --- lifecycle and assignment ---
        Cell(long index,
             long x, 
             long y, 
             int num_environmental_factors,
             Landscape& landscape, 
             const SpeciesPool& species, 
             RandomNumberGenerator& rng);
        ~Cell();
        
        // --- geospatial ---
        long get_index() const {
            return this->_index;
        }
        long get_x() const {
            return this->_x;
        }
        long get_y() const {
            return this->_y;
        }
        
        // --- abiotic ---
        long get_carrying_capacity(long cc) const {
            return this->_carrying_capacity;
        }        
        void set_carrying_capacity(long cc) {
            this->_carrying_capacity = cc;
        }        
        int add_environment_factor(EnvironmentalFactor e) {
            this->_environment.push_back(e);
            return this->_environment.size();
        }
        void set_environment_factor(int idx, EnvironmentalFactor e) {
            assert(idx < this->_enviroment.size());
            this->_environment[idx] = e;
        }
        EnvironmentalFactor get_environment_factor(int idx) const {
            assert(idx < this->_enviroment.size());
            return this->_environment[idx];
        }
        int get_num_environmental_factors() const {
            return this->_environment.size();
        }
        int get_movement_cost() const {
            return this->_movement_cost;
        }
        int set_movement_cost(int cost) {
            this->_movement_cost = cost;
        }              
        
        // --- biotic ---
        long num_organisms() const {
            return this->_organisms.size();
        }
        void add_new_organisms(int species_index, long num) {
            this->_organisms.reserve(this->_organisms.size() + num);
            for ( ; num > 0; --num) {
                this->_organisms.push_back(this->_species.at(species_index)->new_organism());
            }
        }
        
    private:
        // disable copying/assignment
        const Cell& operator=(const Cell& cell);
        Cell(const Cell& cell);        
        
    private:        
        long                         _carrying_capacity;    // max # ind
        long                         _index;                // cell index
        long                         _x;                    // x-coordinate
        long                         _y;                    // y-coordinate
        EnvironmentalFactors        _environment;           // environmental factors
        int                         _movement_cost;         // the base movement penalty for entering this cell
        Organisms                   _organisms;             // the individual organisms of this biota
        
        Landscape&                  _landscape;             // host landscape
        const SpeciesPool&          _species;               // species pool
        RandomNumberGenerator&      _rng;                   // random number generator

}; // Cell

///////////////////////////////////////////////////////////////////////////////	
//! The landscape.
class Landscape {

    public:
    
        // --- lifecycle and assignment ---        
        Landscape(const SpeciesPool& species, RandomNumberGenerator& rng);
        ~Landscape();
        
        // --- initialization and set up ---
        void generate(long size_x, long size_y, int num_environmental_factors); 
 
        // --- landscape access, control and mutation ---                                      
        void set_cell_carrying_capacity(long carrying_capacity);
        
        // --- cell access and spatial mapping ---
        Cell& operator()(long x, long y) {
            return *this->_cells[this->xy_to_index(x, y)];
        }
        Cell& operator[](long index) {
            return *this->_cells[index];
        }            
        Cell& at(long x, long y) {
            return *this->_cells.at(this->xy_to_index(x, y));
        }
        Cell& at(long index) {
            return *this->_cells.at(index);
        }
        long index_to_x(long index) const {
            return index % this->_size_x;          
        }
        long index_to_y(long index) const {
            return static_cast<long>(index / this->_size_x);            
        }
        long xy_to_index(long x, long y) const {
            return (y * this->_size_x) + x;
        }
        long size() const {
            return this->_size_x * this->_size_y;
        }
        
        // --- debugging ---
        void dump() {
            for (long y = 0; y < this->_size_y; ++y) {
                for (long x = 0; x < this->_size_x; ++x) {
                    std::cout <<  this->operator()(x,y).num_organisms() << " ";
                }
                std::cout << std::endl;
            }
        }
                
    private:
        long                        _size_x;                // size of the landscape in the x dimension
        long                        _size_y;                // size of the landscape in the y dimension        
        long                        _size;                  // == x * y, cached here
        std::vector<Cell*>          _cells;                 // cells of the landscape
        
        const SpeciesPool&          _species;               // species pool
        RandomNumberGenerator&      _rng;                   // random number generator        
};

///////////////////////////////////////////////////////////////////////////////	
//! The world.
class World {

    public:
    
        // --- lifecycle --

        World();
        World(unsigned long seed);
        ~World();
        
        // --- access and mutation ---

        RandomNumberGenerator& rng() {
            return this->_rng;
        }  
        Landscape& landscape() {
            return this->_landscape;
        }
        int get_num_fitness_factors() const {
            return this->_num_fitness_factors;
        }
        
        // --- initialization and set up ---
        
        // Should be set up in order: 
        // 1. landscape (with must be generated so cells exist
        // 2. number of fitness factors (= number of environmental factors
        //    = number of genotype factors) must be set
        // 3. movement costs must be defined for each cell 
        //    so that when species are added the movement surface can be 
        //    calculated.
        // 4. carrying capacity should be set,
        // In actual implementation, will be handled by a "scheduler" object.
        
        void generate_landscape(long size_x, long size_y, int num_environmental_factors);
        void set_environmental_factors(); 
        void set_movement_costs();                
        void set_cell_carrying_capacity(long carrying_capacity);
                        
        // to kick start
        Species& new_species(const char *label);        
        void seed_population(long x, long y, int species_index, long size);
        void seed_population(long cell_index, int species_index, long size);
        
        // --- simulation cycles ---
        void cycle();
        
    private:
        SpeciesPool                         _species_pool;
        RandomNumberGenerator               _rng;
        Landscape                           _landscape;        
        int                                 _num_fitness_factors;
        unsigned long                       _current_generation;                
        
}; // World

/******************************************************************************
 * DEFINITIONS
 *****************************************************************************/

///////////////////////////////////////////////////////////////////////////////	
// Species

// --- lifecycle and assignment ---

Species::Species(int index,
                 const char* label, 
                 int num_fitness_factors,
                 RandomNumberGenerator& rng) 
    : _index(index),
      _label(label),
      _num_fitness_factors(num_fitness_factors),
      _rng(rng) {
    this->_index = index;
    this->_movement_rate = 1;
    this->_mutation_rate = 0.1;
    this->_max_mutation_size = 1;
    this->_mean_reproductive_rate = 6;
    this->_reproductive_rate_mutation_size = 1;    
}

Species::Species(const Species& species)
    : _index(species._index),
      _label(species._label),
      _num_fitness_factors(species._num_fitness_factors),
      _rng(species._rng) {
    *this = species;
}

const Species& Species::operator=(const Species& species) {
    this->_label = species._label;
    this->_index = species._index;
    this->_movement_rate = species._movement_rate;
    this->_mutation_rate = species._mutation_rate;
    this->_max_mutation_size = species._max_mutation_size;
    this->_mean_reproductive_rate = species._mean_reproductive_rate;
    this->_reproductive_rate_mutation_size = species._reproductive_rate_mutation_size;    
    this->_default_genotype = species._default_genotype;
    return *this;
}

///////////////////////////////////////////////////////////////////////////////	
// Cell

// --- lifecycle and assignment ---

Cell::Cell(long index, 
           long x, 
           long y, 
           int num_environmental_factors,
           Landscape& landscape, 
           const SpeciesPool& species, 
           RandomNumberGenerator& rng)     
    : _index(index),
      _x(x),
      _y(y),
      _landscape(landscape),
      _species(species),
      _rng(rng) {      
    this->_carrying_capacity = 0;
    this->_movement_cost = 1;
    this->_environment.assign(num_environmental_factors, 0.0);
}

///////////////////////////////////////////////////////////////////////////////	
// Landscape

// --- lifecycle and assignment --- 

Landscape::Landscape(const SpeciesPool& species, RandomNumberGenerator& rng)
    : _species(species),
      _rng(rng) {
    this->_size_x = 0;
    this->_size_y = 0;
    this->_size = 0 * 0;
}

// clean up cells
Landscape::~Landscape() {
    for (std::vector<Cell*>::iterator cell = this->_cells.begin();
            cell != this->_cells.end();
            ++cell) {
        delete *cell; 
    }  
}

// --- initialization and set up ---

void Landscape::generate(long size_x, long size_y, int num_environmental_factors) {
    this->_size_x = size_x;
    this->_size_y = size_y;
    this->_size = size_x * size_y;
    this->_cells.reserve(this->_size);
    for (long x = 0, index = 0; x < size_x; ++x) {
        for (long y = 0; y < size_y; ++y, ++index) {
            Cell* cell = new Cell(index, x, y, num_environmental_factors, *this, this->_species, this->_rng);
            this->_cells.push_back(cell);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////	
// World

// --- lifecycle and assignment --- 

//! default constructor
World::World()
    : _species_pool(),
      _rng(time(0)),
      _landscape(_species_pool, _rng) {
    this->_current_generation = 0;
}

//! constructor: calls
World::World(unsigned long seed) 
    : _species_pool(),
      _rng(seed),
      _landscape(_species_pool, _rng) {
    this->_current_generation = 0;    
}    

//! clean up species pool
World::~World() {
    for (std::vector<Species*>::iterator sp = this->_species_pool.begin();
            sp != this->_species_pool.end();
            ++sp) {
        delete *sp;            
    }            
}

// --- initialization and set up ---

//! Creates a new landscape.
void World::generate_landscape(long size_x, long size_y, int num_environmental_factors) {
    this->_num_fitness_factors = num_environmental_factors;
    this->_landscape.generate(size_x, size_y, num_environmental_factors);
}

//! Actual implementation will load from file(s).
void World::set_environmental_factors() {

}

//! Sets the (uniform) carrying capacity for all the cells on the landscape.
void World::set_cell_carrying_capacity(long carrying_capacity) {
    for (long i = 0; i < this->_landscape.size(); ++i) {
        this->_landscape[i].set_carrying_capacity(carrying_capacity);
    }        
}

//! Adds a new species definition to this world.
Species& World::new_species(const char* label) {
    Species* sp = new Species(this->_species_pool.size(),
                              label, 
                              this->_num_fitness_factors, 
                              this->_rng);
    this->_species_pool.push_back(sp);
    return *sp;
}

//! Populates the cell at (x,y) with organisms of the given species.
void World::seed_population(long x, long y, int species_index, long size) {
    this->_landscape.at(x, y).add_new_organisms(species_index, size);
}

//! Populates the cell cell_index with organisms of the given species.
void World::seed_population(long cell_index, int species_index, long size) {
    this->_landscape.at(cell_index).add_new_organisms(species_index, size);
}

/******************************************************************************
 * MAIN
 *****************************************************************************/

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cout << "usage: " << argv[0] <<  " <DIM-X> <DIM-Y> <CELL-CARRYING-CAPACITY> <NUM-CELLS-TO-POPULATE> <NUM-GENS>\n";
        exit(1);
    }

    int size_x = atoi(argv[1]);
    int size_y = atoi(argv[2]);
    int cc = atoi(argv[3]);
    int num_cells_init = atoi(argv[4]);
    int num_gens = atoi(argv[5]);
    int num_env_factors = 4;
    
    World   world;

//##DEBUG##
DEBUG_BLOCK( std::cout << "(generating landscape)\n"; )
    
	world.generate_landscape(size_x, size_y, num_env_factors);
	
//##DEBUG##
DEBUG_BLOCK( std::cout << "(setting carrying capacity)\n"; )

	world.set_cell_carrying_capacity(cc);

//##DEBUG##
DEBUG_BLOCK( std::cout << "(adding species)\n"; )	
	
	Species& sp1 = world.new_species("gecko");
	sp1.selection_strengths().assign(num_env_factors, 1);
	sp1.default_genotype().assign(num_env_factors, world.rng().randint(-2, 2));
	
//##DEBUG##
DEBUG_BLOCK( std::cout << "(seeding populations)\n"; )
    
    long max_index = (size_x * size_y)-1;
    long cell_index = 0;
    for (std::set<long> seeded; num_cells_init > 0; --num_cells_init) {
        do {
            cell_index = world.rng().randint(0, max_index);
        } while ((seeded.find(cell_index) != seeded.end()) and seeded.size() < max_index+1);
        if (seeded.size() >= max_index+1) {
            break;
        }
        seeded.insert(cell_index);
        world.seed_population(cell_index,
                              sp1.get_index(),
                              cc);
    }
    
    world.landscape().dump();

}





