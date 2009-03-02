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
// You should have received a copy of the GNU General Public License aCellIndexType
// with this program. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#include <algorithm>
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
#include <iterator>
#include <utility>

#define DEBUG_BLOCK(y) y;  


/************************* SUPPORT CLASSES AND METHODS ***********************/

///////////////////////////////////////////////////////////////////////////////
//! Wraps random number generator seed etc.
class RandomNumberGenerator {

    public:
        RandomNumberGenerator();
        RandomNumberGenerator(unsigned long seed);
        void set_seed(unsigned long seed);

        float random();   // [0, 1)
        long randint(int a, int b);
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
long RandomNumberGenerator::randint(int a, int b) {
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
typedef std::vector<Species*> SpeciesPool;
typedef std::vector<Organism> Organisms;
typedef std::vector<Community> Communities;

class Landscape;
class Cell;
class World;
typedef std::vector<Cell> Cells;
typedef long CellIndexType;

typedef int FitnessFactorType;
typedef std::vector<FitnessFactorType> FitnessFactors;


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
        
        Organism(unsigned species_index, const FitnessFactors& new_genotype, Organism::Sex new_sex) 
            : _species_index(species_index),
              _genotype(new_genotype),
              _sex(new_sex),
              _fitness() {
            this->_expired = false;              
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
            this->_expired = ind._expired;
            return *this;
        }
                   
        // genotype       
        const FitnessFactors& genotype() const {
            return this->_genotype;
        }
        
        // fitness & survival
        float &fitness() {
            return this->_fitness;
        }
        bool is_expired() const {
            return this->_expired;
        }
        void set_expired(bool val) {
            this->_expired = val;
        }          
        
        // meta-info
        unsigned species_index() const {
            return this->_species_index;
        }
        
        bool is_male() const {
            return this->_sex == Organism::Male;
        }
        
        bool is_female() const {
            return this->_sex == Organism::Female;
        }                

    private:
        unsigned        _species_index;     // species
        FitnessFactors  _genotype;          // non-neutral genotype: maps to fitness phenotype
        Organism::Sex   _sex;               // male or female
        float           _fitness;           // cache this organism's fitness
        bool            _expired;           // flag an organism to be removed allowing for use of std::remove_if() and std::resize() or v.erase()
        
}; // Organism

///////////////////////////////////////////////////////////////////////////////
//! A collection of processes and properties that determine the ecologies of
//! organisms.
class Species {

    public:
    
        // --- lifecycle and assignment ---        
        Species(unsigned index,
                const char* label, 
                unsigned num_fitness_factors,
                RandomNumberGenerator& rng);              
        ~Species() {}        
                
        // --- access and mutation ---
        unsigned get_index() const {
            return this->_index;
        }
        void set_num_fitness_factors(unsigned i) {
            this->_num_fitness_factors = i;
        }
        unsigned get_num_fitness_factors() const {
            return this->_num_fitness_factors;
        }
        void set_index(unsigned i) {
            this->_index = i;
        }            
        float get_mutation_rate() const {
            return this->_mutation_rate;
        }
        void set_mutation_rate(float i) {
            this->_mutation_rate = i;
        }
        FitnessFactorType get_max_mutation_size() const {
            return this->_max_mutation_size;
        }
        void set_max_mutation_size(FitnessFactorType i) {
            this->_max_mutation_size = i;
        }
        unsigned get_mean_reproductive_rate() const {
            return this->_mean_reproductive_rate;
        }
        void set_mean_reproductive_rate(unsigned i) {
            this->_mean_reproductive_rate = i;
        }     
        unsigned get_reproductive_rate_mutation_size() const {
            return this->_reproductive_rate_mutation_size;
        }
        void set_reproductive_rate_mutation_size(unsigned i) {
            this->_reproductive_rate_mutation_size = i;
        }
        int get_movement_capacity() const {
            return this->_movement_capacity;
        }
        void set_movement_capacity(int i) {
            this->_movement_capacity = i;
        }   
        
        void set_movement_costs(const std::vector<int>& costs) {
            this->_movement_costs = costs;
        }
        int movement_cost(CellIndexType i) {
            assert( (i >= 0 ) and (static_cast<unsigned>(i) < this->_movement_costs.size()) );
            return this->_movement_costs[i];
        }
        
        void set_selection_strengths(const std::vector<float>& strengths) {
            this->_selection_strengths = strengths;
        }
        void set_default_genotype(const FitnessFactors& genotype) {
            this->_default_genotype = genotype;
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
        
        Organism new_organism(const Organism& female, const Organism& male) const {
            return Organism(this->_index, 
                            this->compose_offspring_genotype(female.genotype(), male.genotype()),
                            this->get_random_sex());
        }        
        
        FitnessFactors compose_offspring_genotype(const FitnessFactors& female_genotype, 
                const FitnessFactors& male_genotype) const {
            FitnessFactors offspring_genotype;
            assert(male_genotype.size() == female_genotype.size());
            offspring_genotype.reserve(female_genotype.size());
            FitnessFactors::const_iterator male_g = male_genotype.begin();
            FitnessFactors::const_iterator female_g = female_genotype.begin();
            for ( ; female_g != female_genotype.end(); ++male_g, ++female_g) {
                FitnessFactorType genotype_value = this->_rng.choice(*male_g, *female_g);  
                if (this->_rng.random() < this->_mutation_rate) {
                    genotype_value += this->_rng.randint(-this->_max_mutation_size,
                                                         this->_max_mutation_size);
                }
                offspring_genotype.push_back(genotype_value);
            }
            return offspring_genotype;
        }        
                                
    private:
        // declared as private (and undefined) to prevent copying/assignment
        Species(const Species& species);    
        const Species& operator=(const Species& species);
        
    private:        
        unsigned                    _index;                     // "slot" in cell's pop vector    
        std::string                 _label;                     // arbitrary identifier
        unsigned                    _num_fitness_factors;       // so genotypes of appropriate length can be composed                
        std::vector<float>          _selection_strengths;       // weighted_distance = distance / (sel. strength)
        float                       _mutation_rate;             // rate of mutations
        FitnessFactorType           _max_mutation_size;         // window "size" of mutations
        unsigned                    _mean_reproductive_rate;    // "base" reproductive rate
        unsigned                    _reproductive_rate_mutation_size;  // if reprod. rate evolves, size of step
        std::vector<int>            _movement_costs;            // the movement surface: the "cost" to enter into every cell on the landscape
        int                         _movement_capacity;         // number of cells per round an individual can move
        FitnessFactors              _default_genotype;          // genotype of individuals generated de novo        
        RandomNumberGenerator&      _rng;                       // rng to use

}; // Species

/* LANDSCAPE AND SPATIAL RELATIONS *******************************************/

///////////////////////////////////////////////////////////////////////////////
//! The fundamental atomic spatial unit of the world.
class Cell {
    public:
    
        // --- lifecycle and assignment ---
        Cell(CellIndexType index,
             CellIndexType x, 
             CellIndexType y, 
             unsigned num_environmental_factors,
             Landscape& landscape, 
             const SpeciesPool& species, 
             RandomNumberGenerator& rng);
        ~Cell() {};
        
        // --- geospatial ---
        CellIndexType get_index() const {
            return this->_index;
        }
        CellIndexType get_x() const {
            return this->_x;
        }
        CellIndexType get_y() const {
            return this->_y;
        }
        
        // --- abiotic ---
        CellIndexType get_carrying_capacity() const {
            return this->_carrying_capacity;
        }        
        void set_carrying_capacity(CellIndexType cc) {
            this->_carrying_capacity = cc;
        }        
        unsigned add_environment_factor(FitnessFactorType e) {
            this->_environment.push_back(e);
            return this->_environment.size();
        }
        void set_environment_factor(unsigned idx, FitnessFactorType e) {
            assert(idx < this->_environment.size());
            this->_environment[idx] = e;
        }
        FitnessFactorType get_environment_factor(unsigned idx) const {
            assert(idx < this->_environment.size());
            return this->_environment[idx];
        }
        unsigned get_num_environmental_factors() const {
            return this->_environment.size();
        }            
        
        // --- biotic ---
        CellIndexType num_organisms() const {
            return this->_organisms.size();
        }
        void generate_new_organisms(unsigned species_index, CellIndexType num) {
            this->_organisms.reserve(this->_organisms.size() + num);
            for ( ; num > 0; --num) {
                this->_organisms.push_back(this->_species.at(species_index)->new_organism());
            }
        }
        void insert_organism(const Organism organism) {
            this->_organisms.push_back(organism);
        }
                
        // --- primary biogeographical and evolutionary processes ---
        void survival();
        void competition();
        void reproduction();
        void migration();
        
        // --- supporting biogeographical and evolutionary processes ---
        void extract_breeding_groups(unsigned species_index, 
            std::vector<const Organism*>& female_ptrs,
            std::vector<const Organism*>& male_ptrs) const;
        
    private:
        // disable copying/assignment
        const Cell& operator=(const Cell& cell);
        Cell(const Cell& cell);        
        
    private:        
        CellIndexType               _carrying_capacity;     // max # ind
        CellIndexType               _index;                 // cell index
        CellIndexType               _x;                     // x-coordinate
        CellIndexType               _y;                     // y-coordinate
        FitnessFactors              _environment;           // environmental factors
        Organisms                   _organisms;             // the individual organisms of this biota
        
        Landscape&                  _landscape;             // host landscape
        const SpeciesPool&          _species;               // species pool
        RandomNumberGenerator&      _rng;                   // random number generator
        
        static std::vector<const Organism*> breeding_female_ptrs;    // scratch space for breeding
        static std::vector<const Organism*> breeding_male_ptrs;      // scratch space for breeding
        static Organisms new_offspring;                         // scratch space for next gen

}; // Cell

///////////////////////////////////////////////////////////////////////////////	
//! The landscape.
class Landscape {

    public:
    
        typedef std::pair<Organism, CellIndexType>   MigrationEvent;
        typedef std::vector<MigrationEvent>          MigrationEvents;
    
        // --- lifecycle and assignment ---        
        Landscape(const SpeciesPool& species, RandomNumberGenerator& rng);
        ~Landscape();
        
        // --- initialization and set up ---
        void generate(CellIndexType size_x, CellIndexType size_y, unsigned num_environmental_factors); 
 
        // --- landscape access, control and mutation ---                                      
        void set_cell_carrying_capacity(CellIndexType carrying_capacity);
        
        // --- cell access and spatial mapping ---
        Cell& operator()(CellIndexType x, CellIndexType y) {
            return *this->_cells[this->xy_to_index(x, y)];
        }
        Cell& operator[](CellIndexType index) {
            return *this->_cells[index];
        }            
        Cell& at(CellIndexType x, CellIndexType y) {
            return *this->_cells.at(this->xy_to_index(x, y));
        }
        Cell& at(CellIndexType index) {
            return *this->_cells.at(index);
        }
        CellIndexType index_to_x(CellIndexType index) const {
            return index % this->_size_x;          
        }
        CellIndexType index_to_y(CellIndexType index) const {
            return static_cast<CellIndexType>(index / this->_size_x);            
        }
        CellIndexType xy_to_index(CellIndexType x, CellIndexType y) const {
            return (y * this->_size_x) + x;
        }
        CellIndexType size() const {
            return this->_size;
        }
        CellIndexType size_x() const {
            return this->_size_x;
        }
        CellIndexType size_y() const {
            return this->_size_y;
        }
        CellIndexType random_neighbor(CellIndexType i) {
            static CellIndexType x = 0;
            static CellIndexType y = 0;
            
//  *** if CellIndexType is unsigned ***            
//             x = this->index_to_x(i);
//             if ( x == 0 ) {
//                 x += this->_rng.randint(0, 1);
//             } else if ( x == this->_size_x-1 ) {
//                 x += this->_rng.randint(-1, 0);
//             } else {
//                 x += this->_rng.randint(-1, 1);
//             }
//             
//             y = this->index_to_y(i);
//             if ( y == 0 ) {
//                 y += this->_rng.randint(0, 1);
//             } else if ( y == this->_size_y-1 ) {
//                 y += this->_rng.randint(-1, 0);
//             } else {
//                 y += this->_rng.randint(-1, 1);
//             }            
            
            x = this->index_to_x(i) + this->_rng.randint(-1, 1); // to reflect: % this->_size_x; 
            y = this->index_to_y(i) + this->_rng.randint(-1, 1); // to reflect: % this->_size_y; 
            if (x >= this->_size_x) {
                x = this->_size_x - 1;
            } else if (x < 0) {
                x = 0;
            }
            if (y >= this->_size_y) {
                y = this->_size_y - 1;
            } else if (y < 0) {
                y = 0;
            }
            return this->xy_to_index(x, y);
        }
        
        // --- migration and movement ---
        void clear_migrants() {
            this->_migrants.clear();
        }
        void add_migrant(const Organism& organism, CellIndexType dest_cell) {
            this->_migrants.push_back(std::make_pair(organism, dest_cell));
        }
        void process_migrants() {
            for (MigrationEvents::iterator m = this->_migrants.begin();
                    m != this->_migrants.end();
                    ++m) {
                this->_cells.at(m->second)->insert_organism(m->first);                
            }
            this->_migrants.clear();
        }
        
        // --- debugging ---
        unsigned long dump(std::ostream& output = std::cout) {
            unsigned long num = 0;
            unsigned long total = 0;
            for (CellIndexType y = 0; y < this->_size_y; ++y) {
                for (CellIndexType x = 0; x < this->_size_x; ++x) {
                    num = this->operator()(x,y).num_organisms();
                    total += num;
                    output << std::setw(4) << num << " ";
                }
                output << std::endl;
            }
            output << "---\nTotal organisms: " << total << std::endl;            
            return total;
        }        
                
    private:
        CellIndexType               _size_x;                // size of the landscape in the x dimension
        CellIndexType               _size_y;                // size of the landscape in the y dimension        
        CellIndexType               _size;                  // == x * y, cached here
        std::vector<Cell*>          _cells;                 // cells of the landscape
        
        // cache to collect migrating organisms over a round of migration
        MigrationEvents             _migrants;
        
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
        unsigned get_num_fitness_factors() const {
            return this->_num_fitness_factors;
        }
        void set_num_fitness_factors(unsigned num_fitness_factors) {
            this->_num_fitness_factors = num_fitness_factors;
        }        
        
        // --- initialization and set up ---
        void generate_landscape(CellIndexType size_x, CellIndexType size_y, unsigned num_environmental_factors);
        
        // actual implementation will load this from a file, as will
        // cell environment and species travel costs
        void set_cell_carrying_capacity(CellIndexType carrying_capacity) {
            for (CellIndexType i = 0; i < this->_landscape.size(); ++i) {
                this->_landscape[i].set_carrying_capacity(carrying_capacity);
            }        
        }
        
        // --- species configuration ---
        void set_species_movement_costs(unsigned species_index, const std::vector<int>& costs) {
            assert(species_index < this->_species_pool.size());
            assert(costs.size() == static_cast<unsigned long>(this->_landscape.size()));
            this->_species_pool[species_index]->set_movement_costs(costs);
        }
        void set_species_selection_strengths(unsigned species_index, const std::vector<float>& strengths) {
            assert(species_index < this->_species_pool.size());        
            assert(strengths.size() == this->_num_fitness_factors);
            this->_species_pool[species_index]->set_selection_strengths(strengths);
        }
        void set_species_default_genotype(unsigned species_index, const FitnessFactors& genotype) {
            assert(species_index < this->_species_pool.size());        
            assert(genotype.size() == this->_num_fitness_factors);
            this->_species_pool[species_index]->set_default_genotype(genotype);
        }        
                                
        // to kick start
        Species& new_species(const char *label);        
        void seed_population(CellIndexType x, CellIndexType y, unsigned species_index, CellIndexType size);
        void seed_population(CellIndexType cell_index, unsigned species_index, CellIndexType size);
        
        // --- simulation cycles ---
        void cycle();
        void run(int num_generations);
        
    private:
        SpeciesPool                         _species_pool;
        RandomNumberGenerator               _rng;
        Landscape                           _landscape;        
        unsigned                            _num_fitness_factors;
        CellIndexType                       _current_generation;                
        
}; // World

/******************************************************************************
 * DEFINITIONS
 *****************************************************************************/

///////////////////////////////////////////////////////////////////////////////	
// Species

// --- lifecycle and assignment ---

Species::Species(unsigned index,
                 const char* label, 
                 unsigned num_fitness_factors,
                 RandomNumberGenerator& rng) 
    : _index(index),
      _label(label),
      _num_fitness_factors(num_fitness_factors),
      _rng(rng) {
    this->_index = index;
    this->_mutation_rate = 0.1;
    this->_max_mutation_size = 1;
    this->_mean_reproductive_rate = 6;
    this->_reproductive_rate_mutation_size = 1;
    this->_selection_strengths.assign(this->_num_fitness_factors, 1);
    this->_default_genotype.assign(this->_num_fitness_factors, 0.0);
    this->_movement_capacity = 1;
}


///////////////////////////////////////////////////////////////////////////////	
// Cell

std::vector<const Organism*> Cell::breeding_female_ptrs; // scratch space for breeding
std::vector<const Organism*> Cell::breeding_male_ptrs;   // scratch space for breeding
Organisms Cell::new_offspring;             // scratch space for next gen

// --- lifecycle and assignment ---

Cell::Cell(CellIndexType index, 
           CellIndexType x, 
           CellIndexType y, 
           unsigned num_environmental_factors,
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
    this->_environment.assign(num_environmental_factors, 0.0);
}

// --- primary biogeographical and evolutionary processes ---

void Cell::survival() {

}

void Cell::competition() {

}

void Cell::reproduction() {
    for (SpeciesPool::const_iterator sp = this->_species.begin(); sp != this->_species.end(); ++sp) {
        Cell::breeding_female_ptrs.clear();
        Cell::breeding_male_ptrs.clear();
        Cell::new_offspring.clear();
        
        // species-level reproduction rate for now: later this will be at the 
        // organism level and subject to evolution
        unsigned num_offspring = (*sp)->get_mean_reproductive_rate();
        
        this->extract_breeding_groups((*sp)->get_index(), Cell::breeding_female_ptrs, Cell::breeding_male_ptrs);
        if ( (Cell::breeding_female_ptrs.size() == 0) or (Cell::breeding_male_ptrs.size() == 0)) {
            this->_organisms.clear();
        } else {
            for (std::vector<const Organism*>::iterator fptr = Cell::breeding_female_ptrs.begin();
                    fptr != Cell::breeding_female_ptrs.end();
                    ++fptr) {
                for (unsigned n = 0; n <= num_offspring; ++n) {                    
                    const Organism* male = this->_rng.choice(breeding_male_ptrs);
                    this->new_offspring.push_back((*sp)->new_organism(*(*fptr), *male));
                }                    
            }
            // this->_organisms.swap = Cell::new_offspring
            this->_organisms.swap(Cell::new_offspring);
        }         
    }
}

void Cell::migration() {

    for (Organisms::iterator og = this->_organisms.begin(); og != this->_organisms.end(); ++og) {
    
        assert(og->species_index() < this->_species.size());          
        assert(!og->is_expired());        
        
        Species& sp = *this->_species[og->species_index()];
        int movement = sp.get_movement_capacity();
        CellIndexType curr_idx = this->_index;
                        
        while (movement > 0) {
            CellIndexType dest_idx = this->_landscape.random_neighbor(curr_idx);
            movement -= sp.movement_cost(dest_idx);
            if (movement >= 0) {
                curr_idx = dest_idx;
            }
        } 
        
        if (curr_idx != this->_index) {
            this->_landscape.add_migrant(*og, curr_idx);
            og->set_expired(true);            
        }

    }
    Organisms::iterator end_unexpired = std::remove_if(this->_organisms.begin(), 
        this->_organisms.end(), 
        std::mem_fun_ref(&Organism::is_expired));
    this->_organisms.erase(end_unexpired, this->_organisms.end());
}

// --- supporting biogeographical and evolutionary processes ---

//! Extracts pointers to male and female organisms of a particular species
void Cell::extract_breeding_groups(unsigned species_index, 
                                        std::vector<const Organism*>& female_ptrs,
                                        std::vector<const Organism*>& male_ptrs) const {
    assert(species_index < this->_species.size());    
    for (Organisms::const_iterator og = this->_organisms.begin(); og != this->_organisms.end(); ++og) {
        if (og->species_index() == species_index) {
            if (og->is_female()) {
                female_ptrs.push_back(&(*og));
            } else {
                male_ptrs.push_back(&(*og));
            }
        }
    }
}                                         

///////////////////////////////////////////////////////////////////////////////	
// Landscape

// --- lifecycle and assignment --- 

Landscape::Landscape(const SpeciesPool& species, RandomNumberGenerator& rng)
    : _species(species),
      _rng(rng) {
    this->_size_x = 0;
    this->_size_y = 0;
    this->_size = 0;
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

void Landscape::generate(CellIndexType size_x, CellIndexType size_y, unsigned num_environmental_factors) {
    this->_size_x = size_x;
    this->_size_y = size_y;
    this->_size = size_x * size_y;
    this->_cells.reserve(this->_size);
    for (CellIndexType x = 0, index = 0; x < size_x; ++x) {
        for (CellIndexType y = 0; y < size_y; ++y, ++index) {
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
void World::generate_landscape(CellIndexType size_x, CellIndexType size_y, unsigned num_environmental_factors) {
    this->_num_fitness_factors = num_environmental_factors;
    this->_landscape.generate(size_x, size_y, num_environmental_factors);
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
void World::seed_population(CellIndexType x, CellIndexType y, unsigned species_index, CellIndexType size) {
    this->_landscape.at(x, y).generate_new_organisms(species_index, size);
}

//! Populates the cell cell_index with organisms of the given species.
void World::seed_population(CellIndexType cell_index, unsigned species_index, CellIndexType size) {
    this->_landscape.at(cell_index).generate_new_organisms(species_index, size);
}

// --- species configuration ---


// --- simulation cycles ---
void World::cycle() {
    for (CellIndexType i = this->_landscape.size()-1; i >= 0; --i) {
        this->_landscape[i].survival();
        this->_landscape[i].competition();
        this->_landscape[i].reproduction(); 
        this->_landscape[i].migration();
    }
    this->_landscape.process_migrants();
}

void World::run(int num_generations) {
    for ( ; num_generations >= 0; --num_generations, ++(this->_current_generation)) {
//##DEBUG##
DEBUG_BLOCK( std::cerr << "\n#### GENERATION " << this->_current_generation << " ####\n\n"; )
DEBUG_BLOCK( this->_landscape.dump(std::cout); )        
        this->cycle();        
    }
}

/******************************************************************************
 * MAIN
 *****************************************************************************/

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "usage: " << argv[0] <<  " <DIM-X> <DIM-Y> <CELL-CARRYING-CAPACITY> <NUM-CELLS-TO-POPULATE> <NUM-GENS>\n";
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
DEBUG_BLOCK( std::cerr << "(generating landscape)\n"; )
    
	world.generate_landscape(size_x, size_y, num_env_factors);
	
//##DEBUG##
DEBUG_BLOCK( std::cerr << "(setting carrying capacity)\n"; )

	world.set_cell_carrying_capacity(cc);

//##DEBUG##
DEBUG_BLOCK( std::cerr << "(adding species)\n"; )	
	
	Species& sp1 = world.new_species("gecko");
	
	std::vector<int> costs;
	costs.assign(size_x * size_y, 1);
	sp1.set_movement_costs(costs);
	sp1.set_movement_capacity(1);
	
	std::vector<FitnessFactorType> genotype;
	genotype.reserve(num_env_factors);
	for (int i = 0; i < num_env_factors; ++i) {
	    genotype.push_back(static_cast<FitnessFactorType>(world.rng().randint(-10, 10)));
	}
		
//##DEBUG##
DEBUG_BLOCK( std::cerr << "(seeding populations)\n"; )
    
    unsigned long max_index = (size_x * size_y)-1;
    CellIndexType cell_index = 0;
    for (std::set<CellIndexType> seeded; num_cells_init > 0; --num_cells_init) {
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
    
//##DEBUG##
DEBUG_BLOCK( std::cerr << "(running cycles)\n"; )
    world.run(num_gens);


    std::cerr << "\n#### FINAL STATUS ####\n\n"; 
    world.landscape().dump(std::cerr);
}





