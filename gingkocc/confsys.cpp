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

/* PARSE LOGIC:
 *
 * text => structured text [ConfigurationBlock] => structured values [XXXConf]
 *
 * (1) Data file parsed into [ConfigurationBlock] objects.
 * (2) Each [ConfigurationBlock] object = structured string representation of 
 *     data (e.g., name, type, and dictionary mapping string to strings).
 * (3) Each [ConfigurationBlock] object then is mapped into a corresponding
 *     WorldConfigurator, SpeciesConf, etc. objects, which take the strings and convert
 *     them into values of the appropriate type.
 * (4) The configure_world() functions then take the structured values and 
 *     populate/configure the World object correspondingly.
 */

///////////////////////////////////////////////////////////////////////////////
// Client code should call one of the following to configure World objects.

World& configure_world(World& world, std::istream& conf_src) {
    confsys_detail::ConfigurationFile cf(conf_src);
    cf.configure(world);
    return world;
}

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

ConfigurationFile::ConfigurationFile(std::istream& src)
        : src_(src) {
    if (not this->src_) {
        throw ConfigurationIOError("invalid source stream");
    }
}

ConfigurationFile::ConfigurationFile(const char * fpath)
        : config_filepath_(fpath),
          fsrc_(fpath),
          src_(fsrc_) {
    if (not this->src_) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw ConfigurationIOError(msg.str());
    }          
}

ConfigurationFile::ConfigurationFile(const std::string& fpath) 
        : config_filepath_(fpath),
          fsrc_(fpath.c_str()),
          src_(fsrc_) {
    if (not this->src_) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw ConfigurationIOError(msg.str());
    }
}

ConfigurationFile::~ConfigurationFile() { }

void ConfigurationFile::configure(World& world) {
    assert(this->src_);
    unsigned num_worlds = 0;
    std::set<std::string> species_labels;
    std::set<unsigned long> generations;
}
    
} // confsys_detail

} // namespace confsys

} // namespace gingko
