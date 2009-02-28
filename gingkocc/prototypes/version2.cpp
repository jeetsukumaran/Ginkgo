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
        RandomNumberGenerator(unsigned int seed);
        void set_seed(unsigned int seed);

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
        unsigned int _seed;

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
 * FOREWARD DECLARATIONS AND GLOBAL TYPEDEFS
 *****************************************************************************/

class Organism;
class Species;
class Community;
typedef int GenotypeFactor;
typedef std::vector<int> GenotypeFactors;
typedef std::vector<Organism> Organisms;
typedef std::vector<Species> SpeciesCollection;
typedef std::vector<Community> Communities;

class Landscape;
class Cell;
class World;
typedef int EnvironmentalFactor;
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
    
        // -- lifecycle and assignment --
        Species();
        Species(const char* label);
        Species(const Species& species);
        virtual ~Species() {}        
        const Species& operator=(const Species& species);

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
        RandomNumberGenerator*      _rng;                       // pointer to rng
        int                         _num_environmental_factors; // so genotypes of appropriate length can be composed        

}; // Species

/* LANDSCAPE AND SPATIAL RELATIONS *******************************************/

///////////////////////////////////////////////////////////////////////////////
//! The fundamental atomic spatial unit of the world.
class Cell {
    public:
    
        // lifecycle and assignment
        Cell(Landscape& landscape, int index, int x, int y);
        Cell(const Cell& cell);
        const Cell& operator=(const Cell& cell);

    private:
        int                         _carrying_capacity;     // max # ind
        int                         _index;                 // cell index
        int                         _x;                     // x-coordinate
        int                         _y;                     // y-coordinate
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
        Landscape(World& world, int dim_x, int dim_y);
        Landscape(const Landscape& landscape);
        const Landscape& operator=(const Landscape& landscape);
        
        // --- accessor and mutators ---
        
        Cells& cells() {
            return this->_cells;
        }
                
                                
        // --- landscape access, control and mutation                                
        void generate(World& world, int dim_x, int dim_y);        
        void set_cell_carrying_capacity(int carrying_capacity);        
        
    private:
        int             _size_x;    // size of the landscape in the x dimension
        int             _size_y;    // size of the landscape in the y dimension        
        Cells           _cells;     // cells of the landscape
        World*          _world;     // pointer to the host world

};

///////////////////////////////////////////////////////////////////////////////	
//! The world.
class World {

    public:

        
    private:
        Cells*                  _cells;
        Landscape               _landscape;
        SpeciesCollection       _species_pool;
        RandomNumberGenerator   _rng;
        int                     _num_environmental_factors;
        unsigned long           _current_generation;                
        
}; // World

/******************************************************************************
 * DEFINITIONS
 *****************************************************************************/

/* POPULATION ECOLOGY AND GENETICS *******************************************/



/******************************************************************************
 * MAIN
 *****************************************************************************/

int main(int argc, char* argv[]) {

}





