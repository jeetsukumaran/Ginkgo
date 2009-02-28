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
typedef std::vector<Organism> Organisms;
typedef std::vector<Species> SpeciesCollection;
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
        Species(const char* label = "Sp");
        Species(const Species& species);
        virtual ~Species() {}        
        const Species& operator=(const Species& species);
        
        // --- setup and initialization ---
        void initialize(const char* label);
        void set_world(World& world);

    private:
        std::string                 _label;                     // arbitrary identifier
        int                         _index;                     // "slot" in cell's pop vector
        int                         _movement_rate;             // modifier to the global movement surface to derive the species-specific movement surface
        std::vector<float>          _selection_strengths;       // weighted_distance = distance / (sel. strength)
        std::vector<unsigned>       _movement_surface;          // the movement surface: the "cost" to enter into every cell on the landscape
        float                       _mutation_rate;             // rate of mutations
        float                       _max_mutation_size;         // window "size" of mutations
        int                         _mean_reproductive_rate;    // "base" reproductive rate
        int                         _reproductive_rate_mutation_size;  // if reprod. rate evolves, size of step
        GenotypeFactors             _default_genotype;          // genotype of individuals generated de novo
        
        World*                      _world;                     // pointer to world
        Landscape*                  _landscape;                 // pointer to landscape
        RandomNumberGenerator*      _rng;                       // pointer to rng
        int                         _num_environmental_factors; // so genotypes of appropriate length can be composed        

}; // Species

/* LANDSCAPE AND SPATIAL RELATIONS *******************************************/

///////////////////////////////////////////////////////////////////////////////
//! The fundamental atomic spatial unit of the world.
class Cell {
    public:
    
        // --- lifecycle and assignment ---
        Cell(Landscape& landscape, long index, long x, long y);
        Cell(const Cell& cell);
        ~Cell();
        const Cell& operator=(const Cell& cell);
        
        // --- accessors and mutators
        void set_landscape_ref(Landscape& landscape);

    private:
        long                         _carrying_capacity;    // max # ind
        long                         _index;                // cell index
        long                         _x;                    // x-coordinate
        long                         _y;                    // y-coordinate
        EnvironmentalFactors        _environment;           // environmental factors
        Organisms                   _organisms;             // the individual organisms of this biota
        
        const SpeciesCollection*    _species_pool;          // all the species in the landscape        
        Landscape*                  _landscape;             // host landscape
        World*                      _world;                 // host world

}; // Cell

///////////////////////////////////////////////////////////////////////////////	
//! The landscape.
class Landscape {

    public:
    
        // --- lifecycle and assignment ---        
        Landscape();
        Landscape(World& world, long dim_x, long dim_y);
        ~Landscape();
        
        // --- initialization and set up ---
        void generate(World& world, long size_x, long size_y); 
        
        // --- accessor and mutators ---        
        Cells& cells() {
            return this->_cells;
        }
        
        SpeciesCollection& species_pool();
                                
        // --- landscape access, control and mutation                                       
        void set_cell_carrying_capacity(long carrying_capacity);        
        
    private:
        long            _size_x;   // size of the landscape in the x dimension
        long            _size_y;   // size of the landscape in the y dimension        
        Cells           _cells;     // cells of the landscape
        World*          _world;     // pointer to the host world
};

///////////////////////////////////////////////////////////////////////////////	
//! The world.
class World {

    public:
    
        // --- lifecycle --
        World();
        World(unsigned long seed);
        
        // --- access and mutation ---
        RandomNumberGenerator& rng() {
            return this->_rng;
        }
        SpeciesCollection& species_pool() {
            return this->_species_pool;
        }     
        Landscape& landscape() {
            return this->_landscape;
        }
        int get_num_environmental_factors() const {
            return this->_num_environmental_factors;
        }
        
        // --- initialization and set up ---
        void generate_landscape(long size_x, long size_y);
        void set_cell_carrying_capacity(long carrying_capacity);
        void add_species(const Species& species);
        void seed_population(long x, long y, int species_index, long size);
        void seed_population(long cell_index, int species_index, long size);
        
        // --- simulation cycles ---
        void cycle();
        
    private:
        Landscape               _landscape;
        SpeciesCollection       _species_pool;
        RandomNumberGenerator   _rng;
        int                     _num_environmental_factors;
        unsigned long           _current_generation;                
        
}; // World

/******************************************************************************
 * DEFINITIONS
 *****************************************************************************/

///////////////////////////////////////////////////////////////////////////////	
// Species

// --- lifecycle and assignment ---

Species::Species(const char* label) {
    this->initialize(label);
}

Species::Species(const Species& species) {
    *this = species;
}

const Species& Species::operator=(const Species& species) {
    this->set_world(*species._world);
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

// --- setup and initialization ---

void Species::initialize(const char* label) {
    this->_label = label;
    this->_index = -1;
    this->_movement_rate = 1;
    this->_mutation_rate = 0.1;
    this->_max_mutation_size = 1;
    this->_mean_reproductive_rate = 6;
    this->_reproductive_rate_mutation_size = 1;
}

void Species::set_world(World& world) {
    if (&world == NULL) {
        this->_world = NULL;
        this->_landscape = NULL;
        this->_rng = NULL;
        this->_num_environmental_factors = NULL;
    } else {
        this->_world = &world;
        this->_landscape = &world.landscape();
        this->_rng = &world.rng();
        this->_num_environmental_factors = world.get_num_environmental_factors();                
    }
}

///////////////////////////////////////////////////////////////////////////////	
// Cell

// --- lifecycle and assignment ---

Cell::Cell(Landscape& landscape, long index, long x, long y)
    : _index(index),
      _x(x),
      _y(y) {
    this->set_landscape_ref(landscape);      
}

Cell::Cell(const Cell& cell) {
    *this = cell;
}

Cell::~Cell() {
}

const Cell& Cell::operator=(const Cell& cell) {
    this->set_landscape_ref(*cell._landscape);
    this->_index = cell._index;
    this->_x = cell._x;
    this->_y = cell._y;
    return *this;
}

// --- accessors and mutators ---

void Cell::set_landscape_ref(Landscape& landscape) {
    if (&landscape != NULL) {
        this->_landscape = &landscape;
        this->_species_pool = &landscape.species_pool();
    } else {
        this->_landscape = NULL;
        this->_species_pool = NULL;
    }
}

///////////////////////////////////////////////////////////////////////////////	
// Landscape

// --- lifecycle and assignment --- 

Landscape::Landscape() {

}

Landscape::Landscape(World& world, long size_x, long size_y) {
    this->generate(world, size_x, size_y);
}

Landscape::~Landscape() {
}

// --- accessor and mutators ---        

SpeciesCollection& Landscape::species_pool() {
    assert(this->_world != NULL);
    return this->_world->species_pool();
}

// --- initialization and set up ---

void Landscape::generate(World& world, long size_x, long size_y) {
    this->_world = &world;
    this->_size_x = size_x;
    this->_size_y = size_y;
    long num_cells = size_x * size_y;
    for (long x = 0, index = 0; x < size_x; ++x) {
        for (long y = 0; y < size_y; ++y, ++index) {
            this->_cells.push_back(Cell(*this, index, x, y));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////	
// World

// --- lifecycle and assignment --- 

//! default constructor
World::World()
    : _rng(time(0)) {
    this->_current_generation = 0;
}

//! constructor: calls
World::World(unsigned long seed) 
    : _rng(seed) { 
    this->_current_generation = 0;    
}    

// --- initialization and set up ---


//! Creates a new landscape.
void World::generate_landscape(long size_x, long size_y) {

}

//! Sets the (uniform) carrying capacity for all the cells on the landscape.
void World::set_cell_carrying_capacity(long carrying_capacity) {

}

//! Adds a new species definition to this world.
void World::add_species(const Species& species) {

}

//! Populates the cell at (x,y) with organisms of the given species.
void World::seed_population(long x, long y, int species_index, long size) {

}

//! Populates the cell cell_index with organisms of the given species.
void World::seed_population(long cell_index, int species_index, long size) {

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
    
    World   world;

//##DEBUG##
DEBUG_BLOCK( std::cout << "(generating landscape)\n"; )
    
	world.generate_landscape(size_x, size_y);
	
//##DEBUG##
DEBUG_BLOCK( std::cout << "(setting carrying capacity)\n"; )

	world.set_cell_carrying_capacity(cc);

//##DEBUG##
DEBUG_BLOCK( std::cout << "(adding species)\n"; )	
	
	world.add_species(Species("snail"));
	
//##DEBUG##
DEBUG_BLOCK( std::cout << "(seeding populations)\n"; )

    for (; num_cells_init > 0; --num_cells_init) {
        world.seed_population(world.rng().randint(0, size_x-1),
                              world.rng().randint(0, size_y-1),
                              cc);
    }
    

}





