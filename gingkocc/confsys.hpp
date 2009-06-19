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


#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <map>
#include <utility>
#include <stdexcept>
#include <vector>

#include "world.hpp"
#include "convert.hpp"
#include "Markup.h"

#if !defined(GINGKO_CONFSYS_H)
#define GINGKO_CONFSYS_H

namespace gingko {
namespace confsys {

///////////////////////////////////////////////////////////////////////////////
// Client code should call one of the following to configure World objects.

/**
 * Build/populate a World object according to a given configuration source file.
 * @param   world       World object to configure
 * @param   conf_fpath  path of configuration file
 * @return              World object
 */
World& configure_world(World& world, const char * conf_fpath);

/**
 * Build/populate a World object according to a given configuration source file.
 * @param   world       World object to configure 
 * @param   conf_fpath  path of configuration file
 * @return              World object
 */
World& configure_world(World& world, const std::string& conf_fpath);

///////////////////////////////////////////////////////////////////////////////
// ERRORS AND EXCEPTIONS

/**
 * General configuration error.
 */
class ConfigurationError : public std::runtime_error {
    public:
        ConfigurationError(const char * msg) : std::runtime_error(msg) {}
        ConfigurationError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

/**
 * General i/o error.
 */
class ConfigurationIOError : public ConfigurationError {
    public:
        ConfigurationIOError(const char * msg) : ConfigurationError(msg) {}
        ConfigurationIOError(const std::string& msg) : ConfigurationError(msg) {}
};

/**
 * General configuration format error.
 */
class ConfigurationSyntaxError : public ConfigurationError {
    public:
        ConfigurationSyntaxError(const char * msg) : ConfigurationError(msg) {}
        ConfigurationSyntaxError(const std::string& msg) : ConfigurationError(msg) {}
};

/**
 * General configuration content error.
 */
class ConfigurationIncompleteError : public ConfigurationError {
    public:
        ConfigurationIncompleteError(const char * msg) : ConfigurationError(msg) {}
        ConfigurationIncompleteError(const std::string& msg) : ConfigurationError(msg) {}
};

///////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION DETAILS

namespace confsys_detail {

typedef CMarkup XmlElementType;

/**
 * Track distribution of organisms, either for seed populations or sampling 
 * regimes.
 */
struct OrganismDistribution {
    
    public:
        OrganismDistribution()
            : num_organisms_per_cell(0) { }

    public:
        std::string     species_label;
        unsigned long   num_organisms_per_cell;
        unsigned long   ancestral_population_size;
        unsigned long   ancestral_generations;
        std::vector<CellIndexType>  x;
        std::vector<CellIndexType>  y;
        

}; // OrganismDistribution

/**
 * Encapsulates parsing of a configuration file, and populating of WorldConfigurator,
 * SpeciesConf, GenerationConf, etc. objects.
 */
class ConfigurationFile {

    public:
        
        /**
         * Initializes metadata and binds to source file.
         *
         * @param fpath filepath of data source
         */
        ConfigurationFile(const char * fpath);
        
        /**
         * Initializes metadata and binds to source file.
         *
         * @param fpath filepath of data source
         */
        ConfigurationFile(const std::string& fpath);        
        
        /**
         * Default no-op destructor.
         */
        ~ConfigurationFile();    
        
        /**
         * Parses the configuration file, loading data into blocks.
         */
        void configure(World& world);
        
    private:
    
        bool to_world_element() {
            this->xml_.ResetPos();
            this->xml_.IntoElem(); // at GINGKO
            this->xml_.FindChildElem("gingko");
            if (this->xml_.FindChildElem("world")) {
                this->xml_.IntoElem();
                return true;
            } else {
                return false;
            }
        }
        
        bool to_biota_element() {
            to_world_element();
            if (this->xml_.FindChildElem("biota")) {
                this->xml_.IntoElem();
                return true;
            } else {
                return false;
            }
        }
        
    private: 
        
        /** Path to configuration file. */
        std::string         config_filepath_;
        
        /** XML parser. */
        XmlElementType      xml_;      
        
};  

    
} // confsys_detail
} // namespace confsys
} // namespace gingko

#endif
