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

#include <sstream>
#include "organism.hpp"
#include "world.hpp"

namespace ginkgo
{

// singleton instance
OrganismMemoryManager OrganismMemoryManager::instance_;

// allocate a new Organism slot
Organism* OrganismMemoryManager::allocate(){
 	if (this->free_organism_ptrs_.empty()) {

 	    World& world = World::get_instance();
 	    Logger& logger = world.logger();
 	    std::ostringstream o;
 	    o << "Organism memory manager: allocating new block of " << this->block_size_ << " organism slots.";
        logger.debug(o.str());

 		this->organism_pool_.push_back(std::vector<Organism>());
 		std::vector<Organism> & last_pool_element = *(this->organism_pool_.rbegin());
 		last_pool_element.resize(this->block_size_);
 		for (std::vector<Organism>::iterator i = last_pool_element.begin();
                i != last_pool_element.end();
 		        ++i) {
 			this->free_organism_ptrs_.push(&(*i));
        }

 	    std::ostringstream o2;
 	    o2 << "Organism memory manager: pool size is now " << this->organism_pool_.size() << " blocks, for a total of ";
 	    o2 << (this->block_size_ * this->organism_pool_.size()) << " organism slots.";
        logger.debug(o2.str());

 		return this->allocate();
 	}
 	else {
 		Organism* n = this->free_organism_ptrs_.top();
 		this->free_organism_ptrs_.pop();
 		return n;
 	}
}

// return Organism slot to pool
void OrganismMemoryManager::deallocate(Organism *org) {
	org->clear();
	this->free_organism_ptrs_.push(org);
}

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
