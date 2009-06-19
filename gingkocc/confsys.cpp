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

#include <set>
#include <string>
#include <vector>
#include <ctime>

#include "gingko_defs.hpp"
#include "confsys.hpp"
#include "textutil.hpp"
#include "convert.hpp"
#include "world.hpp"
#include "biosys.hpp"
#include "asciigrid.hpp"
#include "filesys.hpp"
#include "xmlParser.h"

namespace gingko {
namespace confsys {


///////////////////////////////////////////////////////////////////////////////
// Client code should call one of the following to configure World objects.

World& configure_world(World& world, const char * conf_fpath) {
    std::string s(conf_fpath);
    return configure_world(world, s);
}

World& configure_world(World& world, const std::string& conf_fpath) {
    std::ifstream f(conf_fpath.c_str());
    confsys_detail::ConfigurationFile cf(conf_fpath);   
    cf.configure(world); 
    return world;
}

///////////////////////////////////////////////////////////////////////////////
// Supporting Classes and Constructs

namespace confsys_detail {


///////////////////////////////////////////////////////////////////////////////
// ConfigurationFile

ConfigurationFile::ConfigurationFile(const char * fpath) {
    this->open(fpath);
//     bool success = this->xml_.Load(fpath);
//     if (!success) {
//         std::ostringstream msg;
//         msg << "invalid source \"" << fpath << "\"";
//         throw ConfigurationIOError(msg.str());    
//     }
}

ConfigurationFile::ConfigurationFile(const std::string& fpath) {
    this->open(fpath.c_str());
//     bool success = this->xml_.Load(fpath);
//     if (!success) {
//         std::ostringstream msg;
//         msg << "invalid source \"" << fpath << "\"";
//         throw ConfigurationIOError(msg.str());    
//     }
}

ConfigurationFile::~ConfigurationFile() { }

void ConfigurationFile::open(const char * fpath) {
    this->xml_ = XMLNode::openFileHelper(fpath,"gingko");
}

void ConfigurationFile::process_world(World& world) {

    XmlElementType world_node = this->xml_.getChildNode("world");
    if (world_node.isEmpty()) {
        throw ConfigurationSyntaxError("world element is missing from configuration file");
    }
    
    world.set_label( this->get_attribute<std::string>(world_node, "label", "GingkoWorld") );
    world.set_random_seed( this->get_attribute<unsigned long>(world_node, "random_seed", time(0)) );
    world.set_generations_to_run( this->get_attribute<unsigned long>(world_node, "num_gens") );

    unsigned fitness_dim = this->get_attribute<unsigned>(world_node, "fitness_dimensions");    
    if (fitness_dim > MAX_FITNESS_FACTORS) {
        std::ostringstream s;
        s << "maximum number of fitness factors allowed is " << MAX_FITNESS_FACTORS;        
        s << ", but requested " << fitness_dim;
        throw ConfigurationError(s.str());
    }
    world.set_num_fitness_factors(fitness_dim);
    
    world.set_fitness_factor_grain(this->get_attribute<unsigned>(world_node, "fitness_grain", 1));     
    
    
    std::string produce_final = this->get_attribute<std::string>(world_node, "suppress_final_output", "False");
    if (produce_final == "True") {
        world.set_produce_final_output(true);
    } else {
        world.set_produce_final_output(false);
    }
    
    world.generate_landscape( this->get_attribute<CellIndexType>(world_node, "x_range"), 
                              this->get_attribute<CellIndexType>(world_node, "y_range") );
}

void ConfigurationFile::process_biota(World& world) {

    XmlElementType bio_node = this->xml_.getChildNode("biota");
    if (bio_node.isEmpty()) {
        throw ConfigurationSyntaxError("biota element is missing from configuration file");
    }
    
    for (unsigned i = 0; i < bio_node.nChildNode("lineage"); ++i) {
        XmlElementType lnode = bio_node.getChildNode("lineage", i);
        std::string lid = this->get_attribute<std::string>(lnode, "id");
        if (world.has_species(lid)) {
            throw ConfigurationError("lineage \"" + lid + "\" defined multiple times");
        }
        Species& lineage = world.new_species(lid);
        
        
        
//     sp.set_default_genotypic_fitness_factors(this->default_genotypic_fitness_factors_);
//     sp.set_mutation_rate(this->mutation_rate_);
//     sp.set_max_mutation_size(this->max_mutation_size_);
//     sp.set_mean_reproductive_rate(this->mean_reproductive_rate_);
//     sp.set_reproductive_rate_mutation_size(this->reproductive_rate_mutation_size_);
//     sp.set_movement_capacity(this->movement_capacity_);
//     sp.set_movement_probability(this->movement_probability_);
    }
    
}

void ConfigurationFile::configure(World& world) {
    this->process_world(world);
    this->process_biota(world);
}
    
} // confsys_detail

} // namespace confsys

} // namespace gingko
