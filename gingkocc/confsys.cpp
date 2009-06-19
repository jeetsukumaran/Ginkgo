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

void ConfigurationFile::configure(World& world) {
    
    if (!this->to_world_element()) {
    
    }
    
//         <world label="gingko1" 
//            x_range = "40" 
//            y_range = "50" 
//            num_gens = "20001" 
//            fitness_dimensions = "5" 
//            fitness_grain ="1"
//            suppress_final_output = "True"
//            random_seed ="2718281828">
}
    
} // confsys_detail

} // namespace confsys

} // namespace gingko
