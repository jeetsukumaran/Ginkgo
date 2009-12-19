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

#include "biosys.hpp"
#include "cell.hpp"

using namespace ginkgo;

///////////////////////////////////////////////////////////////////////////////
// Organism

Species& Organism::species() const {
    return *this->species_;
}

///////////////////////////////////////////////////////////////////////////////
// Species

// --- lifecycle and assignment ---

Species::Species(const std::string& label,
                 unsigned num_fitness_traits,
                 float global_selection_strength,
                 RandomNumberGenerator& rng)
    : label_(label),
      num_fitness_traits_(num_fitness_traits),
      global_selection_strength_(global_selection_strength),
      rng_(rng),
      organism_label_index_(0) {
    this->mean_reproductive_rate_ = 6;
    this->selection_weights_.assign(this->num_fitness_traits_, 1);
    this->fitness_trait_inheritance_sd_.assign(this->num_fitness_traits_, 0.7071068);
    memset(this->default_heritable_fitness_traits_, 0,
        MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));
    this->movement_capacity_ = 1;
}

///////////////////////////////////////////////////////////////////////////////
// Specialization of std::swap when dealing with Organisms (for efficiency)

namespace std
{
    /**
     * Template specialization of std::swap() for Organism objects.
     */
    template<>
    void swap(ginkgo::Organism& o1, ginkgo::Organism& o2) {
        o1.swap(o2);
    }
}
