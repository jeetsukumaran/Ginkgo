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

#if !defined(WORLD_H)
#define WORLD_H

#include "gingko_defs.h"
#include "randgen.h"
#include "cell.h"
#include "landscape.h"

namespace gingko {

/**
 * Meta-framework that binds everything together.
 */
class World {

    public:
    
        // --- lifecycle --

        /** 
         * Constructs a World with a given RNG seed.
         * @param seed  seed for the random number generator
         */
        World(unsigned long seed);
        
        /** 
         * Destructor, destroys Species and frees memory allocated to Species 
         * objects. 
         */
        ~World();
        
        // --- access and mutation ---

        /**
         * Returns reference to this World's RandomNumberGenerator object.
         * @return reference to this World's RandomNumberGenerator object
         */
        RandomNumberGenerator& rng() {
            return this->rng_;
        }  
        
        /**
         * Returns reference to this World's Landscape object.
         * @return reference to this World's Landscape object
         */        
        Landscape& landscape() {
            return this->landscape_;
        }
        
        /**
         * Returns number of active fitness factors.
         * @return number of active fitness factors
         */        
        unsigned get_num_fitness_factors() const {
            return this->num_fitness_factors_;
        }
        
        /**
         * Sets number of active fitness factors.
         * @param num_fitness_factors number of active fitness factors
         */      
        void set_num_fitness_factors(unsigned num_fitness_factors) {
            this->num_fitness_factors_ = num_fitness_factors;
        }        
        
        /**
         * Build a landscape of the specified spatial and environmental 
         * dimensions.
         *
         * @param size_x                the size of the Landscape from a
         *                              geospatial perspective in the 
         *                              x-dimension
         * @param size_y                the size of the Landscape from a
         *                              geospatial perspective in the 
         *                              y-dimension
         * @param num_fitness_factors   the number of fitness factors
         */
        void generate_landscape(CellIndexType size_x, CellIndexType size_y, unsigned num_fitness_factors);
        
        /**
         * Globally set individual cell carrying capacity.
         *
         * @param carrying_capacity the maximum number of organisms that can
         *                          occupy each cell at the end of every 
         *                          generation
         */
        void set_cell_carrying_capacity(unsigned long carrying_capacity) {
            for (CellIndexType i = 0; i < this->landscape_.size(); ++i) {
                this->landscape_[i].set_carrying_capacity(carrying_capacity);
            }        
        }
        
        /**
         * Set the costs for entering particular cells on the landscape for 
         * a particular species.
         *
         * @param species_index index of pointer to the Species object in the
         *                      Species pool of the Landscape/World
         * @param costs         vector of costs for entering a cell, with 
         *                      cost for cell \f$i\f$ in the landscape given
         *                      by element \f$i\f$ in the costs vector
         */
        void set_species_movement_costs(unsigned species_index, const std::vector<int>& costs) {
            assert(species_index < this->species_.size());
            assert(costs.size() == static_cast<unsigned long>(this->landscape_.size()));
            this->species_[species_index]->set_movement_costs(costs);
        }
        
        /**
         * Sets the weight for each fitness factor for a particular species.
         *
         * @param species_index index of pointer to the Species object in the
         *                      Species pool of the Landscape/World
         * @param strengths     vector coeffecients to the Gaussian distance
         *                      equation used to evaluate fitness
         */
        void set_species_selection_strengths(unsigned species_index, const std::vector<float>& strengths) {
            assert(species_index < this->species_.size());        
            assert(strengths.size() == this->num_fitness_factors_);
            this->species_[species_index]->set_selection_strengths(strengths);
        }
        
        /**
         * Sets the default genotypic (inheritable) component of fitness for a 
         * organisms of the given species when generated de novo.
         *
         * @param species_index index of pointer to the Species object in the
         *                      Species pool of the Landscape/World
         * @param genotype      vector genotypic fitness values for a new 
         *                      organism of the given species
         */        
        void set_species_default_genotype(unsigned species_index, const FitnessFactors& genotype) {
            assert(species_index < this->species_.size());        
            this->species_[species_index]->set_default_genotype(genotype);
        }        
                                
        // --- setup, initialization and seeding ---
        
        /**
         * Instantiates a new species on the heap and inserts into species 
         * pool.
         *
         * @param label     unique label identifying species
         */
        Species& new_species(const char *label);        
        
        /**
         * Generates specified number of new Organism objects of the specified 
         * Species, and inserts them into a Cell specified by its geospatial
         * coordinates.
         *
         * @param x              the geospatial x-coordinate of the Cell into 
         *                       which the new Organism objects will added
         * @param y              the geospatial y-coordinate of the Cell into 
         *                       which the new Organism objects will added
         * @param species_index  index of pointer to the Species object in the
         *                       Species pool of the Landscape/World
         * @param size           number of new Organism objects to generate
         */
        void seed_population(CellIndexType x, CellIndexType y, unsigned species_index, unsigned long size);
        
        /**
         * Generates specified number of new Organism objects of the specified 
         * Species, and inserts them into a Cell specified by its vector index.
         *
         * @param cell_index     the vector index of the Cell into which the 
         *                       new Organism objects will be added
         * @param species_index  index of pointer to the Species object in the
         *                       Species pool of the Landscape/World
         * @param size           number of new Organism objects to generate
         */        
        void seed_population(CellIndexType cell_index, unsigned species_index, unsigned long size);
        
        // --- simulation cycles ---
        
        /**
         * A single cycle or generation of the simulation, including the 
         * reproduction, migration, survival and competition phases.
         */
        void cycle();
        
        /**
         * Run multiple cycles or generations of the simulation.
         */
        void run(unsigned long num_generations);
        
    private:
    
        /** Collection of pointers to the Species objects of this World. */
        SpeciesPointerVector                species_;
        /** The RandomNumberGenerator that is used by all objects of this World. */
        RandomNumberGenerator               rng_;
        /** The geospatial framework of this World. */
        Landscape                           landscape_;
        /** The number of dimensions to the fitness function. */
        unsigned                            num_fitness_factors_;
        /** Tracks the number of generations that have been run. */
        unsigned long                       current_generation_;                

    private:
        /** Disabled copy constructor. */
		World(const World &);
		/** Disabled assignment operator. */
		World & operator=(const World &);
    
}; 

} // gingko namespace

#endif
