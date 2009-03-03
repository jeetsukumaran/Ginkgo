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

///////////////////////////////////////////////////////////////////////////////
//! The fundamental atomic spatial unit of the world.
class Cell {
    public:
    
        // --- lifecycle and assignment ---
        Cell(CellIndexType index,
             unsigned num_environmental_factors,
             Landscape& landscape, 
             const SpeciesPointerVector& species, 
             RandomNumberGenerator& rng);
        ~Cell() {};
        
        // --- geospatial ---
        CellIndexType get_index() const {
            return this->index_;
        }

        // --- abiotic ---
        unsigned long get_carrying_capacity() const {
            return this->carrying_capacity_;
        }        
        void set_carrying_capacity(unsigned long cc) {
            this->carrying_capacity_ = cc;
        }        
        void set_environment_factor(unsigned idx, FitnessFactorType e) {
            assert(idx < this->num_fitness_factors_);
            this->environment_[idx] = e;
        }
        FitnessFactorType get_environment_factor(unsigned idx) const {
            assert(idx < this->num_fitness_factors_);
            return this->environment_[idx];
        }
        unsigned get_num_environmental_factors() const {
            return this->num_fitness_factors_;
        }            
        
        // --- basic biotics ---
        
        CellIndexType num_organisms() const {
            return this->organisms_.size();
        }
        
        void generate_new_organisms(unsigned species_index, CellIndexType num);
        
        void insert_organism(const Organism& organism) {
            this->organisms_.push_back(organism);
        }        
        
        void purge_expired_organisms() {
            OrganismVector::iterator end_unexpired = std::remove_if(this->organisms_.begin(), 
                this->organisms_.end(), 
                std::mem_fun_ref(&Organism::is_expired));
            this->organisms_.erase(end_unexpired, this->organisms_.end());
        }        
    
        // --- primary biogeographical and evolutionary processes ---
        
        void reproduction();
        void migration();
        void survival();
        void competition();        
        
        // --- supporting methods ---
        
        void extract_breeding_groups(unsigned species_index, 
            const OrganismVector& organisms,
            std::vector<const Organism *>& female_ptrs,
            std::vector<const Organism *>& male_ptrs) const;
        
    private:
        // disable copying/assignment
        const Cell& operator=(const Cell& cell);
        Cell(const Cell& cell);
        
    private:        
        CellIndexType                           index_;                     // cell index
        unsigned long                           carrying_capacity_;         // max # ind     
        unsigned                                num_fitness_factors_;       // number of environmental factors
        FitnessFactors                          environment_;               // environmental factors
        OrganismVector                          organisms_;                 // the individual organisms of this biota
        
        Landscape&                              landscape_;                 // host landscape
        const SpeciesPointerVector&             species_;                   // species pool
        RandomNumberGenerator&                  rng_;                       // random number generator

    private:        
        static std::vector<const Organism *>    breeding_female_ptrs;       // scratch space for breeding
        static std::vector<const Organism *>    breeding_male_ptrs;         // scratch space for breeding
        static OrganismVector                   previous_gen;               // scratch space for next gen

}; 
// Cell
///////////////////////////////////////////////////////////////////////////////

} // gingko namespace

#endif
