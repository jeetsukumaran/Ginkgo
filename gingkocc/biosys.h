///////////////////////////////////////////////////////////////////////////////
//
// GINGKO Biogeographical and Evolution Simulator.
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

#if !defined(GINGKO_RANDOM_H)
#define GINGKO_RANDOM_H

#include <cassert>
#include <vector>

#include "gingko_defs.h"

namespace gingko {

class Organism;
class Species;

typedef std::vector<Species *>  SpeciesPointerVector;
typedef std::vector<Organism>   OrganismVector;
typedef int                     FitnessFactorType;
typedef FitnessFactorType       FitnessFactors[MAX_FITNESS_FACTORS];

///////////////////////////////////////////////////////////////////////////////
// Tracks an organisms pedigree.
//
class Pedigree {

	public:
	    
		Pedigree()
		: maternal_pedigree_(0L),
		  paternal_pedigree_(0L),
		  reference_count_(1)
		{}
		
		void set_parents(Pedigree * maternal, Pedigree * paternal) {
			this->maternal_pedigree_ = maternal;
			this->paternal_pedigree_ = paternal;
			if (maternal)
				maternal->increment_count();
			if (paternal)
				paternal->increment_count();
		}
		
		~Pedigree() {
			if (this->maternal_pedigree_)
				this->maternal_pedigree_->decrement_count();
			if (this->paternal_pedigree_)
				this->paternal_pedigree_->decrement_count();
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
		Pedigree * maternal_pedigree_;
		Pedigree * paternal_pedigree_;
		unsigned reference_count_;
		
}; 
// Pedigree
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
        
        Organism(unsigned species_index, const FitnessFactors& new_genotype, Organism::Sex new_sex) 
            : pedigree_node_(0L),
              species_index_(species_index),
              sex_(new_sex),
              fitness_(-1),
              expired_(false) {
            memcpy(this->genotype_, new_genotype, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));
		}
        
        //! Copy constructor.
        Organism(const Organism& ind)
        	: pedigree_node_(0L) {
            *this = ind;
        }

        ~Organism() {
            if (this->pedigree_node_)
            	this->pedigree_node_->decrement_count();
        }
        
        //! Assignment.
        const Organism& operator=(const Organism& ind) {
            if (this == &ind) {
                return *this;
            }
            this->species_index_ = ind.species_index_;
            memcpy(this->genotype_, ind.genotype_, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));
            this->sex_ = ind.sex_;
            this->fitness_ = ind.fitness_;
            this->expired_ = ind.expired_;
            if (this->pedigree_node_)
            	this->pedigree_node_->decrement_count();
            this->pedigree_node_ = ind.pedigree_node_;
            if (this->pedigree_node_)
            	this->pedigree_node_->increment_count();
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
				
		void set_parents(Pedigree* maternal, Pedigree* paternal) {
			assert(this->pedigree_node_ == 0L);
			this->pedigree_node_ = new Pedigree();
			this->pedigree_node_->set_parents(maternal, paternal);
		}
		
    private:
        Pedigree*       pedigree_node_;     // track the pedigree of this organism
        unsigned        species_index_;     // species
        FitnessFactors  genotype_;          // non-neutral genotype: maps to fitness phenotype
        Organism::Sex   sex_;               // male or female
        float           fitness_;           // cache this organism's fitness
        bool            expired_;           // flag an organism to be removed allowing for use of std::remove_if() and std::resize() or v.erase()
        
}; // Organism

} // gingko namespace

#endif
