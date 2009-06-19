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
//         msg << "invalineage_id source \"" << fpath << "\"";
//         throw ConfigurationIOError(msg.str());    
//     }
}

ConfigurationFile::ConfigurationFile(const std::string& fpath) {
    this->open(fpath.c_str());
//     bool success = this->xml_.Load(fpath);
//     if (!success) {
//         std::ostringstream msg;
//         msg << "invalineage_id source \"" << fpath << "\"";
//         throw ConfigurationIOError(msg.str());    
//     }
}

ConfigurationFile::~ConfigurationFile() { }

void ConfigurationFile::open(const char * fpath) {
    this->xml_ = XMLNode::openFileHelper(fpath,"gingko");
}

XmlElementType ConfigurationFile::get_child_node(XmlElementType& current_node, const char * node_name, bool required) {
    XmlElementType cnode = this->xml_.getChildNode(node_name);
    if (cnode.isEmpty() && required) {
        std::ostringstream msg;
        msg << "mandatory element \"" << node_name << "\" is missing from configuration file";
        throw ConfigurationSyntaxError(msg.str());
    }
    return cnode;
}

void ConfigurationFile::process_world(World& world) {

    XmlElementType world_node = this->get_child_node(this->xml_, "world");
    world.set_label( this->get_attribute<std::string>(world_node, "label", "GingkoWorld") );
    world.set_random_seed( this->get_attribute<unsigned long>(world_node, "random_seed", time(0)) );
    world.set_generations_to_run( this->get_attribute<unsigned long>(world_node, "num_gens") );
    try {
        world.set_global_cell_carrying_capacity(this->get_attribute<unsigned long>(world_node, "default_cell_carrying_capacity"));
    } catch (ConfigurationIncompleteError& c) {
        // do nothing
    }

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
    XmlElementType bio_node = this->xml_.getChildNode("world").getChildNode("biota");
    if (bio_node.isEmpty()) {
        throw ConfigurationSyntaxError("biota element is missing from configuration file");
    }
    
    for (unsigned i = 0; i < bio_node.nChildNode("lineage"); ++i) {
        XmlElementType lineage_node = bio_node.getChildNode("lineage", i);
        this->process_lineage(lineage_node, world);
    }    
}

void ConfigurationFile::process_lineage(XmlElementType& lineage_node, World& world) {

    // lineage name / id
    std::string lineage_id = this->get_attribute<std::string>(lineage_node, "id");
    if (world.has_species(lineage_id)) {
        throw ConfigurationError("lineage \"" + lineage_id + "\" defined multiple times");
    }
    Species& lineage = world.new_species(lineage_id);
    
    // genotypic fitness factor
    XmlElementType gtf_node = this->get_child_node(lineage_node, "genotypicFitness", false);
    if (!gtf_node.isEmpty()) {
        std::vector<int> gff = this->get_element_vector<int>(gtf_node);
        if (gff.size() != lineage.get_num_fitness_factors()) {
            std::ostringstream msg;
            msg << "expecting " << lineage.get_num_fitness_factors();
            msg << " default genotypic fitness factors, but found ";
            msg << gff.size() << " instead";
            throw ConfigurationError(msg.str());            
        }
        lineage.set_default_genotypic_fitness_factors(gff);
    }            
   
    // selection weights
    XmlElementType sw_node = this->get_child_node(lineage_node, "selectionWeights", false);
    if (!sw_node.isEmpty()) {
        std::vector<float> sw = this->get_element_vector<float>(sw_node);
        if (sw.size() != lineage.get_num_fitness_factors()) {
            std::ostringstream msg;
            msg << "expecting " << lineage.get_num_fitness_factors();
            msg << " default selection weights, but found ";
            msg << sw.size() << " instead";
            throw ConfigurationError(msg.str());            
        }
        lineage.set_selection_weights(sw);
    }            
          
    lineage.set_mutation_rate(this->get_child_node_scalar<float>(lineage_node, "genotypicFitnessMutationRate", 0.0));
    lineage.set_max_mutation_size(this->get_child_node_scalar<float>(lineage_node, "genotypicFitnessMutationSize", 0.0));
    lineage.set_mean_reproductive_rate(this->get_child_node_scalar<float>(lineage_node, "fecundity", 16));
//         lineage.set_reproductive_rate_mutation_size(this->get_child_node_scalar<float>(lineage_node, "fecundityMutationRate", 0.0));
    lineage.set_movement_probability(this->get_child_node_scalar<float>(lineage_node, "movementProbability", 1.0));
    lineage.set_movement_capacity(this->get_child_node_scalar<unsigned>(lineage_node, "movementCapacity", 1));
    
    // seed populations
    XmlElementType seed_pops = lineage_node.getChildNode("seedPopulations");
    if (seed_pops.isEmpty()) {
        throw ConfigurationError("no seed populations defined for lineage \"" + lineage_id + "\"");   
    }
    for (unsigned i = 0; i < seed_pops.nChildNode("seedPopulation"); ++i) {
        XmlElementType pop_node = seed_pops.getChildNode("seedPopulation", i);
        CellIndexType x = this->get_attribute<CellIndexType>(pop_node, "x");
        if (x > world.landscape().size_x()-1) {
            std::ostringstream msg;
            msg << "maximum x-coordinate on landscape is " << world.landscape().size_x();
            msg << " but seed population position for species " << lineage_id << " specifies x-coordinate of " << x;
            throw ConfigurationError(msg.str());
        }        
        CellIndexType y = this->get_attribute<CellIndexType>(pop_node, "y");
        if (y > world.landscape().size_y()-1) {
            std::ostringstream msg;
            msg << "maximum y-coordinate on landscape is " << world.landscape().size_y();
            msg << " but seed population position for species " << lineage_id << " specifies y-coordinate of " << y;
            throw ConfigurationError(msg.str());
        }
        CellIndexType cell_index = world.landscape().xy_to_index(x, y);
        unsigned long size = this->get_attribute<CellIndexType>(pop_node, "size");
        unsigned long ancestral_pop_size = this->get_child_node_scalar<unsigned long>(pop_node, "ancestralPopulationSize");
        unsigned long ancestral_generations = this->get_child_node_scalar<unsigned long>(pop_node, "ancestralGenerations");
        world.add_seed_population(cell_index, lineage_id, size, ancestral_pop_size, ancestral_generations);
    }    
}

void ConfigurationFile::process_environments(World& world) {
    XmlElementType environs = this->xml_.getChildNode("world").getChildNode("environments");
    if (!environs.isEmpty()) {
        for (unsigned i = 0; i < environs.nChildNode("enviroment"); ++i) {
            XmlElementType env_node = environs.getChildNode("enviroment", i);
        }            
    }
}

void ConfigurationFile::configure(World& world) {
    this->process_world(world);
    this->process_biota(world);
}
    
} // confsys_detail

} // namespace confsys

} // namespace gingko
