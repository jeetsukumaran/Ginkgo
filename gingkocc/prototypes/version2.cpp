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
        
        // fitness
        float set_fitness(float fitness) {
            return (this->_fitness = fitness);
        }
        
        void clear_fitness() {
            this->_fitness = -1;
        }                
        
        float get_fitness() const {
            return this->_fitness;
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
        float            _fitness;          // organism's fitness (calc. and cached every round)   

}; // Organism

///////////////////////////////////////////////////////////////////////////////
//! A collection of all interacting organisms in a single cell on the landscape. 
//! This will include multiple species, and thus defines the ecological 
//! community of the cell (i.e., the level which competition occurs). When the 
//! individuals are partitioned into collections of distinct species, each 
//! collection defines the breeding community of this cell.
class Community {

    public:
    
        // -- lifecycle and assignment --
        Community(const Cell& cell);
        Community(const Community& community);
        ~Community() {};
        const Community& operator=(const Community& community);
        
        // -- core accessors/mutators --        
        void set_cell(const Cell& cell);

    private:
        const SpeciesCollection*    _species_pool;  // all the species in the landscape
        const Cell*                 _cell;          // the host cell of this biota
        Organisms                   _organisms;     // the individual organisms of this biota

}; // Community

/******************************************************************************
 * DEFINITIONS
 *****************************************************************************/

/* POPULATION ECOLOGY AND GENETICS *******************************************/

/* Community */

// -- lifecycle and assignment --

Community::Community(const Cell& c) {
    this->set_cell(c);
}

Community::Community(const Community& community) {
    *this = community;
}

const Community& Community::operator=(const Community& community) {
    if (community._cell != NULL) {
        this->set_cell(*community._cell);
    } else {
        this->_cell = NULL;
        this->_species_pool = NULL;    
    }
    this->_organisms = community._organisms;
    return *this;
}

// -- core accessors/mutators --

void Community::set_cell(const Cell& c) {
    this->_cell = NULL;
    this->_species_pool = NULL;
}

/******************************************************************************
 * MAIN
 *****************************************************************************/

int main(int argc, char* argv[]) {

}





