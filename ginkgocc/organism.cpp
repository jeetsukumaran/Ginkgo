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

#include "organism.hpp"

namespace ginkgo
{

// key for sets/maps: sorts by fitness, or random if equal fitness
bool compare_organism_fitness(const Organism * o1, const Organism * o2) {
    float f1 = o1->get_fitness();
    float f2 = o2->get_fitness();
    assert(f1 >= 0);
    assert(f2 >= 0);
    if (f1 == f2) {
        RandomNumberGenerator& rng = RandomNumberGenerator::get_instance();
        return rng.uniform_01() < 0.5;
    } else {
        return f1 < f2;
    }
}

} // namespace ginkgo

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
