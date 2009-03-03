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

#if !defined(WORLD_H)
#define WORLD_H

#include "gingko_defs.h"

namespace gingko {

///////////////////////////////////////////////////////////////////////////////	
//! The world.
class World {

    public:
    
        // --- lifecycle --

        //World();
        World(unsigned long seed);
        ~World();
        
        // --- access and mutation ---

        RandomNumberGenerator& rng() {
            return this->rng_;
        }  
        Landscape& landscape() {
            return this->landscape_;
        }
        unsigned get_num_fitness_factors() const {
            return this->num_fitness_factors_;
        }
        void set_num_fitness_factors(unsigned num_fitness_factors) {
            this->num_fitness_factors_ = num_fitness_factors;
        }        
        
        // --- initialization and set up ---
        void generate_landscape(CellIndexType size_x, CellIndexType size_y, unsigned num_environmental_factors);
        
        // actual implementation will load this from a file, as will
        // cell environment and species travel costs
        void set_cell_carrying_capacity(unsigned long carrying_capacity) {
            for (CellIndexType i = 0; i < this->landscape_.size(); ++i) {
                this->landscape_[i].set_carrying_capacity(carrying_capacity);
            }        
        }
        
        // --- species configuration ---
        void set_species_movement_costs(unsigned species_index, const std::vector<int>& costs) {
            assert(species_index < this->species_pool_.size());
            assert(costs.size() == static_cast<unsigned long>(this->landscape_.size()));
            this->species_pool_[species_index]->set_movement_costs(costs);
        }
        void set_species_selection_strengths(unsigned species_index, const std::vector<float>& strengths) {
            assert(species_index < this->species_pool_.size());        
            assert(strengths.size() == this->num_fitness_factors_);
            this->species_pool_[species_index]->set_selection_strengths(strengths);
        }
        void set_species_default_genotype(unsigned species_index, const FitnessFactors& genotype) {
            assert(species_index < this->species_pool_.size());        
            assert(genotype.size() == MAX_FITNESS_FACTORS);
            this->species_pool_[species_index]->set_default_genotype(genotype);
        }        
                                
        // to kick start
        Species& new_species(const char *label);        
        void seed_population(CellIndexType x, CellIndexType y, unsigned species_index, CellIndexType size);
        void seed_population(CellIndexType cell_index, unsigned species_index, CellIndexType size);
        
        // --- simulation cycles ---
        void cycle();
        void run(int num_generations);
        
    private:
    

        SpeciesPointerVector                         species_pool_;
        RandomNumberGenerator               rng_;
        Landscape                           landscape_;        
        unsigned                            num_fitness_factors_;
        CellIndexType                       current_generation_;                

		World(); // declare but do not define -- require seed
		World(const World &); // declare but do not define -- not copyable
		World & operator=(const World &); // declare but do not define -- not copyable
    
}; 
// World
///////////////////////////////////////////////////////////////////////////////

} // gingko namespace

#endif