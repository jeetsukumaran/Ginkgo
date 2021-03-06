///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Phylogeographical Evolution Simulator.
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
#include "xmlParser.h"

#if !defined(GINKGO_CONFSYS_H)
#define GINKGO_CONFSYS_H

namespace ginkgo {
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

typedef XMLNode XmlElementType;

/**
 * Track distribution of organisms, either for seed populations or sampling
 * regimes.
 */
struct OrganismDistribution {

    public:
        OrganismDistribution()
            : num_organisms_per_cell(0) { }

    public:
        std::string                 species_label;
        PopulationCountType         num_organisms_per_cell;
        PopulationCountType         ancestral_population_size;
        GenerationCountType         ancestral_generations;
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

        void open(const char *fpath);

        bool has_attribute(const XmlElementType& xml, const char * attr_name) const {
            const char * attr_value = xml.getAttribute(attr_name);
            return attr_value != NULL;
        }

        template <typename T>
        T get_attribute(const XmlElementType& xml, const char * attr_name) const {
            const char * attr_value = xml.getAttribute(attr_name);
            if (attr_value == NULL) {
                std::ostringstream msg;
                msg << "mandatory attribute \"" << attr_name << "\" missing for element \"" << xml.getName() << "\"";
                throw ConfigurationIncompleteError(msg.str());
            }
            return convert::to_scalar<T>(attr_value);
        }

        template <typename T>
        T get_attribute(const XmlElementType& xml, const char * attr_name, T default_value) const {
            const char * attr_value = xml.getAttribute(attr_name);
            if (attr_value == NULL) {
                return default_value;
            }
            return convert::to_scalar<T>(attr_value);
        }

        bool get_attribute_bool(const XmlElementType& xml, const char * attr_name, bool default_value) const {
            const char * attr_value = xml.getAttribute(attr_name);
            if (attr_value == NULL) {
                return default_value;
            }
            return convert::to_bool(attr_value);
        }

        template <typename T>
        T get_child_node_scalar(XmlElementType& xml, const char * node_name) const {
            XmlElementType cnode = xml.getChildNode(node_name);
            if (cnode.isEmpty()) {
                std::ostringstream msg;
                msg << "mandatory sub-element \"" << node_name << "\" missing for element \"" << xml.getName() << "\"";
                throw ConfigurationIncompleteError(msg.str());
            }
            std::ostringstream raw;
            for (int i = 0; i < cnode.nText(); ++i) {
                raw << cnode.getText(i);
            }
            return convert::to_scalar<T>(raw.str());
        }

        template <typename T>
        T get_child_node_scalar(XmlElementType& xml, const char * node_name, T default_value) const {
            XmlElementType cnode = xml.getChildNode(node_name);
            if (cnode.isEmpty()) {
                return default_value;
            }
            std::ostringstream raw;
            for (int i = 0; i < cnode.nText(); ++i) {
                raw << cnode.getText(i);
            }
            return convert::to_scalar<T>(raw.str());
        }

        template <typename T>
        T get_element_scalar(XmlElementType& xml) {
            std::ostringstream raw;
            for (int i = 0; i < xml.nText(); ++i) {
                raw << xml.getText(i);
            }
            return convert::to_scalar<T>(raw.str());
        }

        template <typename T>
        std::vector<T> get_element_vector(XmlElementType& xml) {
            std::ostringstream raw;
            for (int i = 0; i < xml.nText(); ++i) {
                raw << xml.getText(i);
            }
            return convert::to_vector_on_any<T>(raw.str(), " \t\r\n", true);
        }

        CellIndexType get_validated_cell_index(CellIndexType x, CellIndexType y, World& world, const char * item_desc);

        template <class T>
        std::string get_validated_grid_path(const std::string& grid_path, const World& world);

        XmlElementType get_child_node(XmlElementType& current_node, const char * node_name, bool required=true);

        void parse_meta(World& world);
        void parse_system(World& world);
        void parse_landscape(World& world);
        void parse_lineages(World& world);
        void parse_lineage(XmlElementType& lnode, World& world);
        void parse_initialization(World& world);
        void parse_environments(World& world);
        void parse_jump_dispersals(World& world);
        void parse_migration_trackings(World& world);
        MigrationTrackingRegime build_migration_tracking_regime(World& world,
                const XmlElementType& mig_track_element);
        void parse_samplings(World& world);

        EnvironmentSettings parse_environment_settings(World& world, const XmlElementType& env_node);
        CellIndexType parse_cell_index_from_node(World& world, XmlElementType& node);

    private:

        /** Path to configuration file. */
        std::string         config_filepath_;

        /** XML parser. */
        XmlElementType      ginkgo_root_;

};


} // confsys_detail
} // namespace confsys
} // namespace ginkgo

#endif
