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

#if !defined(GINGKO_CELL_H)
#define GINGKO_CELL_H

#include "biosys.h"
#include <iostream>

namespace gingko {

class Landscape;

/**
 * The fundamental and atomic spatial unit of the world.
 */
class Cell {

    public:
    
        // --- lifecycle and assignment ---
        
        /**
         * Constructs a cell in a particular Landscape of a particular World, 
         * at a specific position and with the specified number of active
         * environmental fitness factors.
         *
         * @param index     the index of this cell on the Landscape
         * @param num_environmental_factors 
         *                  the number of factors to be considered for fitness
         * @species         reference to the World species pool
         * @rng             reference to the World random number generator
         */
        Cell(CellIndexType index,
             unsigned num_environmental_factors,
             Landscape& landscape, 
             const SpeciesPointerVector& species, 
             RandomNumberGenerator& rng);
             
        /**
         * No-op destructor.
         */
        ~Cell() {};
        
        // --- geospatial ---
        
        /**
         * Returns index of this cell in the landscape.
         * @return index of this cell in the landscape
         */
        CellIndexType get_index() const {
            return this->index_;
        }

        // --- abiotic ---
        
        /** 
         * Returns the maximum number of organisms (of all species) that can
         * occupy this cell at the end of every generation.
         *
         * @return  maximum occupancy of the cell
         */
        unsigned long get_carrying_capacity() const {
            return this->carrying_capacity_;
        }  
        
        /** 
         * Sets the maximum number of organisms (of all species) that can
         * occupy this cell at the end of every generation.
         *
         * @param cc  maximum occupancy of the cell
         */        
        void set_carrying_capacity(unsigned long cc) {
            this->carrying_capacity_ = cc;
        }        
        
        /** 
         * Returns the specified environmental fitness factor for this cell. 
         *
         * @return  environmental fitness factor value
         */        
        FitnessFactorType get_environment_factor(unsigned idx) const {
            assert(idx < this->num_fitness_factors_);
            return this->environment_[idx];
        }           
                 
        /** 
         * Sets the specified environmental fitness factor for this cell. 
         *
         * @param  idx      the the environmental fitness factor to set
         * @param  value    the value to set it to
         */                       
        void set_environment_factor(unsigned idx, FitnessFactorType e) {
            assert(idx < this->num_fitness_factors_);
            this->environment_[idx] = e;
        }
        
        /** 
         * Returns the number of active environmental fitness factor. 
         *
         * @return  number of active environmental fitness factors
         */              
        unsigned get_num_environmental_factors() const {
            return this->num_fitness_factors_;
        } 
        
        // --- basic biotics ---
        
        /** 
         * Returns direct handle to the organisms of this cell. 
         *
         * @return  vector of organisms that occupy this cell.
         */          
        OrganismVector& organisms() {
            return this->organisms_;
        }
        
        /** 
         * Returns the number of organisms that occupy this cell. 
         *
         * @return  number of organisms tha occupy this cell.
         */          
        unsigned long num_organisms() const {
            return this->organisms_.size();
        }
        
        /**
         * Creates the specified number of organisms of the specified species 
         * and adds them to the population of this cell.
         *
         * @param species_index index of the pointer to the Species object in 
         *                      the World species pool
         * @param num           number of organisms to create
         */
        void generate_new_organisms(unsigned species_index, CellIndexType num);
        
        /**
         * Adds a single organism to the population of this cell.
         *
         * @param organism  the organism to be copied into this cell's
         *                  vector of organisms
         */
        void insert_organism(const Organism& organism) {
            this->organisms_.push_back(organism);
        }        
        
        //** Removes organisms flagged for removal from this cell's population */
        void purge_expired_organisms() {
            OrganismVector::iterator end_unexpired = std::remove_if(this->organisms_.begin(), 
                this->organisms_.end(), 
                std::mem_fun_ref(&Organism::is_expired));
            this->organisms_.erase(end_unexpired, this->organisms_.end());
        }        
    
        // --- primary biogeographical and evolutionary processes ---
        
        /** 
         * Implements reproduction, where organisms are paired, offspring are 
         * generated, parents removed from cell's population, and offspring
         * are inserted.
         */
        void reproduction();
        
        /**
         * Shuffles organisms around using brownian motion; organisms that
         * end up in another cell are flagged for removal from this cell and
         * copied to landscape-wide container for insertion into destination
         * cells.
         */
        void migration();
        
        /** 
         * Each organism in this cell is tested for survival given the
         * environment of the current cell, with the probability of survival
         * proportional to the fitness score of the organism.
         */
        void survival();
        
        /** 
         * If the carrying capacity of the cell is exceeded, the organisms of
         * this cell are sorted according to their fitness factor, and only
         * the first \f$K\f$ organisms are retained to reproduce in the
         * following generation, where \f$K\f$ is the carrying capacity of this
         * cell.
         */        
        void competition();        
        
        // --- supporting methods ---
        
        /** 
         * Given a species index, extracts pointers to male and female 
         * organisms of the specified species within this cell.
         *
         * @param species_index     index of species
         * @param organisms         vector of organisms (source)
         * @param female_ptrs       pointers to organisms in source vector
         *                          that are females of the specified species
         * @param male_ptrs         pointers to organisms in source vector
         *                          that are males of the specified species
         */         
        void extract_breeding_groups(unsigned species_index, 
            const OrganismVector& organisms,
            std::vector<const Organism *>& female_ptrs,
            std::vector<const Organism *>& male_ptrs) const;
        
    private:
        /** Copy constructor, disabled by private scoping and no definition */
        Cell(const Cell&);
        /** Assignment, disabled by private scoping and no definition */
        const Cell& operator=(const Cell&);
        
    private:
        /** index of this cell in the landscape */
        CellIndexType                           index_;
        /** maximum number of individual organisms that can occupy this cell */
        unsigned long                           carrying_capacity_;
        /** number of active fitness factors */
        unsigned                                num_fitness_factors_;
        /** vector of environmental fitness factor components in this cell */
        FitnessFactors                          environment_;
        /** collection of organisms occupying this cell */
        OrganismVector                          organisms_;
        /** reference to Landscape in which this cell is located */
        Landscape&                              landscape_;
        /** reference to pool of species in the World of this cell */
        const SpeciesPointerVector&             species_;
        /** reference to random number generator of the World of this cell */
        RandomNumberGenerator&                  rng_;

    private:        
        /** scratch space to sort females of a species for reproduction */
        static std::vector<const Organism *>    breeding_female_ptrs;
        /** scratch space to sort males of a species for reproduction */
        static std::vector<const Organism *>    breeding_male_ptrs;
        /** scratch space to parents while compose next generation */
        static OrganismVector                   previous_gen;

}; 
// Cell
///////////////////////////////////////////////////////////////////////////////

} // gingko namespace

#endif
