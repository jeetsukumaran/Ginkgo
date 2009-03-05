///////////////////////////////////////////////////////////////////////////////
//
// GINGKO Biogeographical Evolution Simulator.
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

#if !defined(GINGKO_BIOSYS_H)
#define GINGKO_BIOSYS_H

#include <cassert>
#include <vector>
#include <string>
#include <cmath>

#include <iostream>


#include "gingko_defs.h"
#include "random.h"

namespace gingko {

class Organism;
class Species;

typedef std::vector<Species *>  SpeciesPointerVector;
typedef std::vector<Organism>   OrganismVector;

///////////////////////////////////////////////////////////////////////////////
// Tracks an organisms pedigree.
//
class HaploidGenealogy {

	public:
	    
		HaploidGenealogy()
		: parent_(0L),
		  reference_count_(1)
		{ }
		
		void inherit(HaploidGenealogy * parent) {
			this->parent_ = parent;
			if (parent)
				parent->increment_count();
		}
		
		~HaploidGenealogy() {
			if (this->parent_)
				this->parent_->decrement_count();
			assert(this->reference_count_ == 0 || this->reference_count_ == 1);
		}
		
		void decrement_count() {
			if (this->reference_count_ == 1)
				delete this; // never do this!!
			this->reference_count_ -= 1;
		}
		
		void increment_count() {
			this->reference_count_ += 1;
		}
		
    private:		
		HaploidGenealogy *      parent_;
		unsigned                reference_count_;
		
}; 
// HaploidGenealogy
///////////////////////////////////////////////////////////////////////////////


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
        
        Organism(unsigned species_index, 
                 unsigned num_fitness_factors, 
                 const FitnessFactors& new_genotype, 
                 Organism::Sex new_sex) 
            : species_index_(species_index),
              num_fitness_factors_(num_fitness_factors),              
              sex_(new_sex),
              fitness_(-1),
              genealogy_(0L),              
              expired_(false) {
            memcpy(this->genotype_, new_genotype, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));
		}
        
        //! Copy constructor.
        Organism(const Organism& ind)
        	: genealogy_(0L) {
            *this = ind;
        }

        ~Organism() {
            if (this->genealogy_)
            	this->genealogy_->decrement_count();
        }
        
        //! Assignment.
        const Organism& operator=(const Organism& ind) {
            if (this == &ind) {
                return *this;
            }
            this->species_index_ = ind.species_index_;
            this->num_fitness_factors_ = ind.num_fitness_factors_;
            memcpy(this->genotype_, ind.genotype_, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));
            this->sex_ = ind.sex_;
            this->fitness_ = ind.fitness_;
            this->expired_ = ind.expired_;
            if (this->genealogy_)
            	this->genealogy_->decrement_count();
            this->genealogy_ = ind.genealogy_;
            if (this->genealogy_)
            	this->genealogy_->increment_count();
            return *this;
        }
        
        // for sorting
        bool operator<(const Organism& other) const {
            assert(this->fitness_ >= 0);
            assert(other.fitness_ >= 0); 
            return this->fitness_ < other.fitness_; 
        } 
                   
        // genotype       
        const FitnessFactors& genotype() const {
            return this->genotype_;
        }
        
        // fitness & survival
        float get_fitness() const {
            return this->fitness_;
        }
        void set_fitness(float fitness) {
            this->fitness_ = fitness;
        }
        bool is_expired() const {
            return this->expired_;
        }
        void set_expired(bool val) {
            this->expired_ = val;
        }          
        
        // meta-info
        unsigned species_index() const {
            return this->species_index_;
        }
        
        bool is_male() const {
            return this->sex_ == Organism::Male;
        }
        
        bool is_female() const {
            return this->sex_ == Organism::Female;
        }                
				
		void inherit_genealogy(const Organism& organism) {		    
			assert(this->genealogy_ == 0L);
			this->genealogy_ = new HaploidGenealogy();
			this->genealogy_->inherit(organism.genealogy_);
		}
		
    private:
        unsigned        species_index_;         // species
        unsigned        num_fitness_factors_;   // number of factors effecting fitness        
        FitnessFactors  genotype_;              // non-neutral genotype: maps to fitness phenotype
        Organism::Sex   sex_;                   // male or female
        float           fitness_;               // cache this organism's fitness
        HaploidGenealogy *      genealogy_;              // track the pedigree of this organism        
        bool            expired_;               // flag an organism to be removed allowing for use of std::remove_if() and std::resize() or v.erase()
        
};
// Organism
///////////////////////////////////////////////////////////////////////////////

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
            return this->index_;
        }
        void set_num_fitness_factors(unsigned i) {
            this->num_fitness_factors_ = i;
        }
        unsigned get_num_fitness_factors() const {
            return this->num_fitness_factors_;
        }
        void set_index(unsigned i) {
            this->index_ = i;
        }            
        float get_mutation_rate() const {
            return this->mutation_rate_;
        }
        void set_mutation_rate(float i) {
            this->mutation_rate_ = i;
        }
        FitnessFactorType get_max_mutation_size() const {
            return this->max_mutation_size_;
        }
        void set_max_mutation_size(FitnessFactorType i) {
            this->max_mutation_size_ = i;
        }
        unsigned get_mean_reproductive_rate() const {
            return this->mean_reproductive_rate_;
        }
        void set_mean_reproductive_rate(unsigned i) {
            this->mean_reproductive_rate_ = i;
        }     
        unsigned get_reproductive_rate_mutation_size() const {
            return this->reproductive_rate_mutation_size_;
        }
        void set_reproductive_rate_mutation_size(unsigned i) {
            this->reproductive_rate_mutation_size_ = i;
        }
        int get_movement_capacity() const {
            return this->movement_capacity_;
        }
        void set_movement_capacity(int i) {
            this->movement_capacity_ = i;
        }   
        
        void set_movement_costs(const std::vector<int>& costs) {
            this->movement_costs_ = costs;
        }
        int movement_cost(CellIndexType i) {
            assert( (i >= 0 ) and (static_cast<unsigned>(i) < this->movement_costs_.size()) );
            return this->movement_costs_[i];
        }
        
        void set_selection_strengths(const std::vector<float>& strengths) {
            this->selection_strengths_ = strengths;
        }
        void set_default_genotype(const FitnessFactors& genotype) {
            memcpy(this->default_genotype_, genotype, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));
        }         
        
        // --- fitness ---
        float calc_fitness(const Organism& organism, const FitnessFactors environment) const {
            const FitnessFactors& genotype = organism.genotype();        
            const FitnessFactorType * g = genotype;
            const FitnessFactorType * e = environment;
            std::vector<float>::const_iterator s = this->selection_strengths_.begin();
            float weighted_distance = 0.0;
            for (unsigned i = 0; i < this->num_fitness_factors_; ++i, ++g, ++e, ++s) {
                weighted_distance += pow((*e - *g), 2) * *s; // each distance weighted by selection strength
            }
            return exp(-weighted_distance);
        }                    
        
        // --- organism generation and reproduction ---
        Organism::Sex get_random_sex(float female_threshold=0.5) const {
            if (this->rng_.uniform_real() < female_threshold) {
                return Organism::Male;
            } else {
                return Organism::Female;
            }
        }  
        
        Organism new_organism() const {
            return Organism(this->index_, this->num_fitness_factors_, this->default_genotype_, this->get_random_sex());
        }
        
        Organism new_organism(const Organism& female, const Organism& male) const {
        	FitnessFactors offspring_genotype;
			this->compose_offspring_genotype(female.genotype(), male.genotype(), offspring_genotype);
            Organism organism(this->index_, this->num_fitness_factors_, offspring_genotype, this->get_random_sex());
            organism.inherit_genealogy(female);
            return organism;
        }        
        
        void compose_offspring_genotype(const FitnessFactors& female_genotype, 
                const FitnessFactors& male_genotype, 
                FitnessFactors & offspring_genotype) const {
				for (unsigned i = 0; i < this->num_fitness_factors_; ++i) {
					FitnessFactorType genotype_value = this->rng_.select(male_genotype[i], female_genotype[i]);  
					if (this->rng_.uniform_real() < this->mutation_rate_) {
						genotype_value += this->rng_.uniform_int(-this->max_mutation_size_,
															 this->max_mutation_size_);
					}
					offspring_genotype[i] = genotype_value;
				}
        }

    private:
        // declared as private (and undefined) to prevent copying/assignment
        Species(const Species& species);    
        const Species& operator=(const Species& species);
        
    private:        
        unsigned                    index_;                     // "slot" in cell's pop vector    
        std::string                 label_;                     // arbitrary identifier
        unsigned                    num_fitness_factors_;       // so genotypes of appropriate length can be composed                
        std::vector<float>          selection_strengths_;       // weighted_distance = distance / (sel. strength)
        float                       mutation_rate_;             // rate of mutations
        FitnessFactorType           max_mutation_size_;         // window "size" of mutations
        unsigned                    mean_reproductive_rate_;    // "base" reproductive rate
        unsigned                    reproductive_rate_mutation_size_;  // if reprod. rate evolves, size of step
        std::vector<int>            movement_costs_;            // the movement surface: the "cost" to enter into every cell on the landscape
        int                         movement_capacity_;         // number of cells per round an individual can move
        FitnessFactors              default_genotype_;          // genotype of individuals generated de novo        
        RandomNumberGenerator&      rng_;                       // rng to use

};
// Species
///////////////////////////////////////////////////////////////////////////////

} // gingko namespace

#endif
