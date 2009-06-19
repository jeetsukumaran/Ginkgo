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
#include "Markup.h"

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
    bool success = this->xml_.Load(fpath);
    if (!success) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw ConfigurationIOError(msg.str());    
    }
}

ConfigurationFile::ConfigurationFile(const std::string& fpath) {
    bool success = this->xml_.Load(fpath);
    if (!success) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw ConfigurationIOError(msg.str());    
    }
}

ConfigurationFile::~ConfigurationFile() { }

void ConfigurationFile::process_world(World& world) {
    if (!this->to_world_element()) {
        throw ConfigurationSyntaxError("world element is missing from configuration file");
    }
    
    world.set_label( this->get_attribute<std::string>("label", "GingkoWorld") );
    world.set_random_seed( this->get_attribute<unsigned long>("random_seed", time(0)) );
    world.set_generations_to_run( this->get_attribute<unsigned long>("num_gens") );

    unsigned fitness_dim = this->get_attribute<unsigned>("fitness_dimensions");    
    if (fitness_dim > MAX_FITNESS_FACTORS) {
        std::ostringstream s;
        s << "maximum number of fitness factors allowed is " << MAX_FITNESS_FACTORS;        
        s << ", but requested " << fitness_dim;
        throw ConfigurationError(s.str());
    }
    world.set_num_fitness_factors(fitness_dim);
    
    world.set_fitness_factor_grain(this->get_attribute<unsigned>("fitness_grain", 1));     
    
    
    std::string produce_final = this->get_attribute<std::string>("suppress_final_output", "False");
    if (produce_final == "True") {
        world.set_produce_final_output(true);
    } else {
        world.set_produce_final_output(false);
    }
    
    world.generate_landscape( this->get_attribute<CellIndexType>("x_range"), 
                              this->get_attribute<CellIndexType>("y_range") );
}

void ConfigurationFile::configure(World& world) {
    this->process_world(world);
}
    
} // confsys_detail

} // namespace confsys

} // namespace gingko
