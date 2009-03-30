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

#include "biosys.hpp"

using namespace gingko;

///////////////////////////////////////////////////////////////////////////////	
// Organism

Species& Organism::species() const {
    return *this->species_;
} 

///////////////////////////////////////////////////////////////////////////////	
// Species

// --- lifecycle and assignment ---

Species::Species(const std::string& label, 
                 unsigned num_fitness_factors,
                 RandomNumberGenerator& rng) 
    : label_(label),
      num_fitness_factors_(num_fitness_factors),
      rng_(rng),
      organism_label_index_(0) {
    this->mutation_rate_ = 0.1;
    this->max_mutation_size_ = 1;
    this->mean_reproductive_rate_ = 6;
    this->reproductive_rate_mutation_size_ = 1;
    this->selection_strengths_.assign(this->num_fitness_factors_, 1);
    memset(this->default_genotypic_fitness_factors_, 0, 
        MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));    
    this->movement_capacity_ = 1;
}
