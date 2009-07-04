///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Biogeographical Evolution Simulator.
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

#if !defined(GINKGO_CELL_H)
#define GINKGO_CELL_H

#include "biosys.hpp"
#include <iostream>

namespace ginkgo {

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
         * @param index     the x geographic coordinate of this cell
         * @param index     the y geographic coordinate of this cell
         * @param num_fitness_factors 
         *                  the number of factors to be considered for fitness
         * @species         reference to the World species pool
         * @rng             reference to the World random number generator
         */
        Cell(CellIndexType index,
             CellIndexType x,
             CellIndexType y,
             unsigned num_fitness_factors,
             Landscape& landscape, 
             const SpeciesByLabel& species, 
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
        
        /**
         * Returns x-coordinate of this cell in the landscape.
         * @return x-coordinate of this cell in the landscape
         */
        CellIndexType get_x() const {
            return this->x_;
        }        
        
        /**
         * Returns y-coordinate of this cell in the landscape.
         * @return y-coordinate of this cell in the landscape
         */
        CellIndexType get_y() const {
            return this->y_;
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
        unsigned get_num_fitness_factors() const {
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
         * @return  number of organisms that occupy this cell.
         */          
        unsigned long num_organisms() const {
            return this->organisms_.size();
        }
        
        /** 
         * Gets the number of male and female organisms of a particular lineage 
         * that occupy this cell. 
         *
         * @param   species_ptr pointer to species for which to filter; if NULL
         *                      then not filtered by species
         * @param   num_females to be loaded with the number of females         
         * @param   num_males   to be loaded with the number of males
         *
         * @return  number of male organisms that occupy this cell.
         */          
        void num_organisms(Species * species_ptr, unsigned long& num_females, unsigned long& num_males) const;       
        
        /** 
         * Returns the number of organisms of a particular species that occupy 
         * this cell. 
         *
         * @param   species_ptr pointer to species
         * @return              number of invidiuals of a particular species 
         *                      that occupy this cell
         */          
        unsigned long num_organisms(Species * sp_ptr) const;       
        
        /**
         * Creates the specified number of organisms of the specified species 
         * and adds them to the population of this cell.
         *
         * @param species       pointer to Species object
         * @param num           number of organisms to create
         */
        void generate_new_organisms(Species * sp, unsigned long num);
        
        
        /**
         * Generates a new population of organisms, and runs reproductive cycles
         * under random/panmictic mating and fixed population size for specified 
         * generations. Useful for "bootstrapping" a seed or colonizing 
         * population of organsms. Typically, you would specify a fairly small 
         * popultion size (e.g. N = 100 - 1000), and run this process for 
         * 10N-20N generations, to ensure coalescence of all the final 
         * generation to a single ancestor.
         *
         * @param sp                       pointer to species
         * @param final_pop_size           number of individuals to sample from 
         *                                 the ancestral population (n)
         * @param ancestral_pop_size       size of ancestral population (N; 0 => N = n )
         * @param ancestral_generations    number of generations to run (0 => 10N)
         */
        void generate_new_population(Species * sp, 
            unsigned long final_pop_size,
            unsigned long ancestral_pop_size = 0,      
            unsigned long ancestral_generations = 0);
        
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
                        
        // --- sampling for tree-building ---

        /** 
         * Adds pointers to organisms of a particular species to the given 
         * vector. If num_organisms is 0 all organisms of that species are 
         * added, otherwise limited to num_organisms, sampled at random. If 
         * num_organisms exceeds the number of organisms of 
         * given species in the cell, then all the organisms are returned.
         * @param sp            pointer to Species object
         * @param samples       vector of pointers to organisms to add to         
         * @param num_organisms number of organisms (0=all)
         */
        void sample_organisms(Species * sp_ptr,
            std::vector<const Organism *>& samples, unsigned long num_organisms);            
        
    private:        
        
        // --- supporting methods ---    
        
        /** 
         * Given a species index, extracts pointers to male and female 
         * organisms of the specified species (or any species, if none specified) 
         * within this cell.
         *
         * @param sp_ptr            pointer to Species object; if not NULL, 
         *                          only organisms of this species will be sampled
         * @param organisms         vector of organisms (source)
         * @param female_ptrs       pointers to organisms in source vector
         *                          that are females of the specified species
         * @param male_ptrs         pointers to organisms in source vector
         *                          that are males of the specified species
         */         
        void extract_breeding_groups(Species * sp_ptr, 
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
        /** x-coordinate of this cell in the landscape */
        CellIndexType                           x_;             
        /** y-coordinate of this cell in the landscape */ 
        CellIndexType                           y_;   
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
        const SpeciesByLabel&                   species_;
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

} // ginkgo namespace

#endif