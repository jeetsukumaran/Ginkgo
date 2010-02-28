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

#include <set>
#include <string>
#include <vector>
#include <ctime>

#include "ginkgo_defs.hpp"
#include "confsys.hpp"
#include "textutil.hpp"
#include "convert.hpp"
#include "world.hpp"
#include "organism.hpp"
#include "population.hpp"
#include "species.hpp"
#include "asciigrid.hpp"
#include "filesys.hpp"
#include "xmlParser.h"

namespace ginkgo {
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

double log_poisson_pdf(unsigned long k, long double mean) {
    double lp = 0.0;
    for (unsigned long i = 1 ; i <= k ; ++i)
        lp -= log((double)i);
    lp += ((double)k)*log(mean) - mean;
    return lp;
}

double poisson_pdf(unsigned long k, long double mean) {
    return exp(log_poisson_pdf(k, mean));
}

///////////////////////////////////////////////////////////////////////////////
// ConfigurationFile

ConfigurationFile::ConfigurationFile(const char * fpath) {
    this->open(fpath);
//     bool success = this->ginkgo_root_.Load(fpath);
//     if (!success) {
//         std::ostringstream msg;
//         msg << "invalineage_id source \"" << fpath << "\"";
//         throw ConfigurationIOError(msg.str());
//     }
}

ConfigurationFile::ConfigurationFile(const std::string& fpath) {
    this->open(fpath.c_str());
//     bool success = this->ginkgo_root_.Load(fpath);
//     if (!success) {
//         std::ostringstream msg;
//         msg << "invalineage_id source \"" << fpath << "\"";
//         throw ConfigurationIOError(msg.str());
//     }
}

ConfigurationFile::~ConfigurationFile() { }

void ConfigurationFile::open(const char * fpath) {
    this->config_filepath_ = fpath;
    this->ginkgo_root_ = XMLNode::openFileHelper(fpath,"ginkgo");
}

void ConfigurationFile::configure(World& world) {
    this->parse_meta(world);
    this->parse_system(world);
    this->parse_landscape(world);
    this->parse_lineages(world);
    this->parse_initialization(world);
    this->parse_environments(world);
//    this->parse_dispersals(world);
    this->parse_samplings(world);
}

void ConfigurationFile::parse_meta(World& world) {
    world.set_title( this->get_attribute<std::string>(this->ginkgo_root_, "title", "GinkgoWorld") );
    world.set_log_frequency(this->get_attribute<unsigned>(this->ginkgo_root_, "log_frequency", 10));
    world.set_allow_multifurcations(this->get_attribute_bool(this->ginkgo_root_, "multifurcating_trees", true));
    world.set_produce_final_output(this->get_attribute_bool(this->ginkgo_root_, "final_output", false));
    world.set_produce_full_complement_diploid_trees(this->get_attribute_bool(this->ginkgo_root_, "full_complement_diploid_trees", false));
}

void ConfigurationFile::parse_system(World& world) {
    XmlElementType sys_node = this->get_child_node(this->ginkgo_root_, "system");
    world.set_random_seed( this->get_child_node_scalar<unsigned long>(sys_node, "random_seed", time(0)) );
    world.set_global_selection_strength(this->get_child_node_scalar<float>(sys_node, "global_selection_strength", DEFAULT_GLOBAL_SELECTION_STRENGTH));
    unsigned fitness_dim = this->get_child_node_scalar<unsigned>(sys_node, "fitness_dimensions");
    if (fitness_dim > MAX_FITNESS_TRAITS) {
        std::ostringstream s;
        s << "maximum number of fitness dimensions allowed is " << MAX_FITNESS_TRAITS;
        s << ", but requested " << fitness_dim;
        throw ConfigurationError(s.str());
    }
    world.set_num_fitness_traits(fitness_dim);
    world.set_generations_to_run( this->get_child_node_scalar<GenerationCountType>(sys_node, "ngens") );
}

void ConfigurationFile::parse_landscape(World& world) {
    XmlElementType landscape_node = this->get_child_node(this->ginkgo_root_, "landscape");
    world.generate_landscape( this->get_attribute<CellIndexType>(landscape_node, "ncols"),
                              this->get_attribute<CellIndexType>(landscape_node, "nrows") );
    world.set_global_cell_carrying_capacity(this->get_child_node_scalar<PopulationCountType>(landscape_node, "default_cell_carrying_capacity", 0));
}

void ConfigurationFile::parse_lineages(World& world) {
    XmlElementType bio_node = this->ginkgo_root_.getChildNode("lineages");
    if (bio_node.isEmpty()) {
        throw ConfigurationSyntaxError("lineages element is missing from configuration file");
    }

    for (int i = 0; i < bio_node.nChildNode("lineage"); ++i) {
        XmlElementType lineage_node = bio_node.getChildNode("lineage", i);
        this->parse_lineage(lineage_node, world);
    }
}

void ConfigurationFile::parse_lineage(XmlElementType& lineage_node, World& world) {

    // lineage name / id
    std::string lineage_id = this->get_attribute<std::string>(lineage_node, "id");
    if (world.has_species(lineage_id)) {
        throw ConfigurationError("lineage \"" + lineage_id + "\" defined multiple times");
    }
    Species& lineage = world.new_species(lineage_id);

    // genotypic fitness factor
    XmlElementType gtf_node = this->get_child_node(lineage_node, "fitness_trait_default_genotypes", false);
    if (!gtf_node.isEmpty()) {
        std::vector<FitnessTraitType> gff = this->get_element_vector<FitnessTraitType>(gtf_node);
        if (gff.size() != lineage.get_num_fitness_traits()) {
            std::ostringstream msg;
            msg << "expecting " << lineage.get_num_fitness_traits();
            msg << " default genotypic fitness trait values, but found ";
            msg << gff.size() << " instead";
            throw ConfigurationError(msg.str());
        }
        lineage.set_default_fitness_trait_genotypes(gff);
    }

    // selection weights
    XmlElementType sw_node = this->get_child_node(lineage_node, "fitness_trait_relative_selection_weights", false);
    if (!sw_node.isEmpty()) {
        std::vector<float> sw = this->get_element_vector<float>(sw_node);
        if (sw.size() != lineage.get_num_fitness_traits()) {
            std::ostringstream msg;
            msg << "expecting " << lineage.get_num_fitness_traits();
            msg << " selection weights, but found ";
            msg << sw.size() << " instead";
            throw ConfigurationError(msg.str());
        }
        lineage.set_selection_weights(sw);
    }

    // fitness sd
//    XmlElementType fsd_node = this->get_child_node(lineage_node, "fitnessTraitInheritanceStdDev", false);
//    if (!fsd_node.isEmpty()) {
//        std::vector<float> sd = this->get_element_vector<float>(fsd_node);
//        if (sd.size() != lineage.get_num_fitness_traits()) {
//            std::ostringstream msg;
//            msg << "expecting " << lineage.get_num_fitness_traits();
//            msg << " fitness inheritance standard deviations, but found ";
//            msg << sd.size() << " instead";
//            throw ConfigurationError(msg.str());
//        }
//        lineage.set_fitness_trait_inheritance_sd(sd);
//    }

    lineage.set_mean_reproductive_rate(this->get_child_node_scalar<unsigned>(lineage_node, "fecundity", 16));

    XmlElementType mc_node = this->get_child_node(lineage_node, "movement_capacity", false);
    if (!mc_node.isEmpty()) {
        std::string dist = this->get_attribute<std::string>(mc_node, "distribution");
        if (dist == "constant") {
            lineage.set_movement_capacity_fixed(this->get_element_scalar<MovementCountType>(mc_node));
        } else if (dist == "user") {
            std::vector<float> mov_probs = this->get_element_vector<float>(mc_node);
            lineage.set_movement_capacity_probabilities(mov_probs);
        } else if (dist == "poisson") {
            int poisson_mean = this->get_element_scalar<MovementCountType>(mc_node);
            std::vector<float> mov_probs;
            mov_probs.reserve(MAX_MOVEMENT_COUNT+1);
            for (MovementCountType i = 0; i <= MAX_MOVEMENT_COUNT; ++i) {
                mov_probs.push_back(poisson_pdf(i, poisson_mean));
            }
            lineage.set_movement_capacity_probabilities(mov_probs);
        } else {
            throw ConfigurationError("movement_capacity for '" + lineage_id + "': 'distribution' attribute must be one of: 'constant', 'poisson', or 'user'");
        }
    } else {
        lineage.set_movement_capacity_fixed(1);
    }

}

void ConfigurationFile::parse_initialization(World& world) {
    InitializationRegime  initialization_regime;
    XmlElementType initialization = this->get_child_node(this->ginkgo_root_, "initialization", true);
    initialization_regime.max_cycles = this->get_attribute<GenerationCountType>(initialization, "max_cycles", 0);
    XmlElementType env_node = this->get_child_node(initialization, "environment", false);

    if (!env_node.isEmpty()) {
        initialization_regime.environment = this->parse_environment_settings(world, env_node);
    }
    XmlElementType cell_pops = this->get_child_node(initialization, "populations", true);
    for (int i = 0; i < cell_pops.nChildNode("cell"); ++i) {
        XmlElementType cell_node = cell_pops.getChildNode("cell", i);
        CellIndexType cell_index = this->parse_cell_index_from_node(world, cell_node);
        for (int i = 0; i < cell_node.nChildNode("population"); ++i) {
            XmlElementType cell_pop_node = cell_node.getChildNode("population", i);
            std::string lineage_id = this->get_attribute<std::string>(cell_pop_node, "lineage");
            PopulationCountType pop_size = this->get_attribute<PopulationCountType>(cell_pop_node, "size");
            if (not world.has_species(lineage_id)) {
                throw ConfigurationError("initialization: lineage \"" + lineage_id + "\" not defined");
            }
            Species * species_ptr = world.species_registry()[lineage_id];
            initialization_regime.cell_populations[cell_index][species_ptr] = pop_size;
        }
    }
    world.set_initialization_regime(initialization_regime);
}

void ConfigurationFile::parse_environments(World& world) {
    XmlElementType environs = this->ginkgo_root_.getChildNode("environments");
    if (!environs.isEmpty()) {
        for (int i = 0; i < environs.nChildNode("environment"); ++i) {
            XmlElementType env_node = environs.getChildNode("environment", i);
            GenerationCountType gen = this->get_attribute<GenerationCountType>(env_node, "gen");
            EnvironmentSettings environment_settings = this->parse_environment_settings(world, env_node);
            world.add_environment_settings(gen, environment_settings);
        }
    }
}

//void ConfigurationFile::parse_dispersals(World& world) {
//    XmlElementType dispersals = this->ginkgo_root_.getChildNode("world").getChildNode("dispersals");
//    if (!dispersals.isEmpty()) {
//        for (int i = 0; i < dispersals.nChildNode("dispersal"); ++i) {
//            XmlElementType disp_node = dispersals.getChildNode("dispersal", i);
//            GenerationCountType gen = this->get_attribute<GenerationCountType>(disp_node, "gen");
//            DispersalEvent disp_event;
//            std::ostringstream item_desc;
//            item_desc << "dispersal event " << i+1 << " in generation " << gen;
//            disp_event.source = this->get_validated_cell_index(this->get_attribute<CellIndexType>(disp_node, "from_x"),
//                    this->get_attribute<CellIndexType>(disp_node, "from_y"),
//                    world,
//                    item_desc.str().c_str());
//            disp_event.destination = this->get_validated_cell_index(this->get_attribute<CellIndexType>(disp_node, "to_x"),
//                    this->get_attribute<CellIndexType>(disp_node, "to_y"),
//                    world,
//                    item_desc.str().c_str());
//            disp_event.probability = this->get_child_node_scalar<float>(disp_node, "probability", 1.0);
//            std::string lineage_id = this->get_child_node_scalar<std::string>(disp_node, "lineage", "");
//            if (lineage_id.size() > 0) {
//                if (not world.has_species(lineage_id)) {
//                    throw ConfigurationError("dispersal: lineage \"" + lineage_id + "\" not defined");
//                }
//                disp_event.species_ptr = world.species_registry()[lineage_id];
//            } else {
//                disp_event.species_ptr = NULL;
//            }
//            world.add_dispersal_event(gen, disp_event);
//        }
//    }
//}

void ConfigurationFile::parse_samplings(World& world) {
    XmlElementType samplings = this->ginkgo_root_.getChildNode("samples");
    if (samplings.isEmpty()) {
        return;
    }
    for (int i = 0; i < samplings.nChildNode("sample"); ++i) {

        // common
        XmlElementType sample_node = samplings.getChildNode("sample", i);
        GenerationCountType gen = this->get_attribute<GenerationCountType>(sample_node, "gen");
        std::string lineage_id = this->get_attribute<std::string>(sample_node, "lineage");
        if (not world.has_species(lineage_id)) {
            throw ConfigurationError("sample: lineage \"" + lineage_id + "\" not defined");
        }

        // occurrence
        world.add_occurrence_sampling(gen, world.species_registry()[lineage_id]);

        // skip to next if trees not wanted
        if (!this->get_attribute_bool(sample_node, "trees", true) ) {
            continue;
        }

        // process tree sampling directive
        SamplingRegime world_sampling_regime;
        world_sampling_regime.species_ptr = world.species_registry()[lineage_id];
        std::string label = this->get_attribute<std::string>(sample_node, "label", "");
        if (label.size() > 0) {
            world_sampling_regime.label = label;
        }
        world_sampling_regime.num_organisms_per_cell = this->get_child_node_scalar<PopulationCountType>(sample_node, "individuals_per_cell", 0);
        XmlElementType cell_nodes = sample_node.getChildNode("cells");
        if (cell_nodes.isEmpty()) {
            continue; // all cells assumed by default
        }
        for (int j = 0; j < cell_nodes.nChildNode("cell"); ++j) {
            XmlElementType sample_cell_node = samplings.getChildNode("cell", i);
            CellIndexType cell_index = this->parse_cell_index_from_node(world, sample_cell_node);
            world_sampling_regime.cell_indexes.insert(cell_index);
        }
        world.add_tree_sampling(gen, world_sampling_regime);
    }
}

EnvironmentSettings ConfigurationFile::parse_environment_settings(World& world, const XmlElementType& env_node) {
    EnvironmentSettings environment_settings;
    for (int j = 0; j < env_node.nChildNode(); ++j) {
        XmlElementType sub_node = env_node.getChildNode(j);
        std::string node_name = sub_node.getName();
        if (node_name == "carrying_capacity") {
            environment_settings.carrying_capacity = this->get_validated_grid_path<PopulationCountType>(this->get_element_scalar<std::string>(sub_node), world);
        } else if (node_name == "fitness_trait_optima") {
            unsigned eidx = this->get_attribute<unsigned>(sub_node, "trait");
            if (eidx > world.get_num_fitness_traits()) {
                std::ostringstream msg;
                msg << "invalid fitness trait index: " << eidx;
                msg << " (maximum valid index is " << world.get_num_fitness_traits();
                msg << ", given " << world.get_num_fitness_traits() << " defined fitness traits)";
                throw ConfigurationError(msg.str());
            }
            std::string gridfile = this->get_validated_grid_path<FitnessTraitType>(this->get_element_scalar<std::string>(sub_node), world);
            environment_settings.fitness_trait_optima.insert(std::make_pair(eidx, gridfile));
        } else if (node_name == "movement_costs") {
            std::string lineage_id = this->get_attribute<std::string>(sub_node, "lineage");
            if (not world.has_species(lineage_id)) {
                throw ConfigurationError("movement costs: lineage \"" + lineage_id + "\" not defined");
            }
            Species * lineage = world.species_registry()[lineage_id];
            std::string gridfile = this->get_validated_grid_path<MovementCountType>(this->get_element_scalar<std::string>(sub_node), world);
            environment_settings.movement_costs.insert(std::make_pair(lineage, gridfile));
        }
//        else if (node_name == "movement_probabilities") {
//            std::string lineage_id = this->get_attribute<std::string>(sub_node, "lineage");
//            if (not world.has_species(lineage_id)) {
//                throw ConfigurationError("movement probabilities: lineage \"" + lineage_id + "\" not defined");
//            }
//            Species * lineage = world.species_registry()[lineage_id];
//            std::string gridfile = this->get_validated_grid_path<float>(this->get_element_scalar<std::string>(sub_node), world);
//            environment_settings.movement_probabilities.insert(std::make_pair(lineage, gridfile));
//        }
    }
    return environment_settings;
}

CellIndexType ConfigurationFile::parse_cell_index_from_node(World& world, XmlElementType& cell_node) {
        bool has_x = this->has_attribute(cell_node, "x");
        bool has_y = this->has_attribute(cell_node, "y");
        bool has_index = this->has_attribute(cell_node, "index");
        if ( (has_x or has_y) and has_index ) {
            throw ConfigurationError("cannot specify cell by both coordinates ('x'/'y') and indexes ('index')");
        } else if (has_x and not has_y) {
            throw ConfigurationError("'y' coordinate not specified for cell");
        } else if (has_y and not has_x)  {
            throw ConfigurationError("'x' coordinate not specified for cell");
        } else if (has_x and has_y) {
            CellIndexType x = this->get_attribute<CellIndexType>(cell_node, "x");
            CellIndexType y = this->get_attribute<CellIndexType>(cell_node, "y");
            return this->get_validated_cell_index(x, y, world, "cell coordinate");
        } else if (has_index) {
            CellIndexType cell_index = this->get_attribute<CellIndexType>(cell_node, "index");
            if (cell_index >= world.landscape().size()) {
                std::ostringstream msg;
                msg << "seed population: maximum cell index is " << world.landscape().size() - 1;
                msg << " (0-based indexing), but cell index of " << cell_index << " specified";
                throw ConfigurationError(msg.str());
            }
            return cell_index;
        } else {
            throw ConfigurationError("Logical switch failure processing cell");
        }
}

template <class T>
std::string ConfigurationFile::get_validated_grid_path(const std::string& grid_path, const World& world) {
    std::string top_dir = filesys::get_path_parent(this->config_filepath_);
    std::string full_grid_path;
    if (filesys::is_abs_path(grid_path) or top_dir.size() == 0) {
        full_grid_path = grid_path;
    } else {
        full_grid_path = filesys::compose_path(top_dir, grid_path);
    }
    try {
        asciigrid::AsciiGrid<T> grid(full_grid_path);
        std::vector<T> values = grid.get_cell_values();
        if (values.size() != world.size()) {
            std::ostringstream msg;
            msg << "landscape has " << world.size() << " cells, ";
            msg << "but grid \"" << full_grid_path << "\" describes " << values.size() << " cells";
            throw ConfigurationError(msg.str());
        }
        return full_grid_path;
    } catch (asciigrid::AsciiGridIOError e) {
        throw ConfigurationIOError("I/O error reading grid \"" + full_grid_path + "\": " + e.what());
    } catch (asciigrid::AsciiGridFormatError e) {
        throw ConfigurationIOError("format error reading grid \"" + full_grid_path + "\": " + e.what());
    }
}

CellIndexType ConfigurationFile::get_validated_cell_index(CellIndexType x,
        CellIndexType y,
        World& world,
        const char * item_desc) {
    if (x > world.landscape().size_x() - 1) {
        std::ostringstream msg;
        msg << "maximum x-coordinate on landscape is " << world.landscape().size_x() - 1;
        msg << " (0-based indexing), but x-coordinate of " << x << " specified";
        if (item_desc != NULL) {
            msg << " for " << item_desc;
        }
        throw ConfigurationError(msg.str());
    }
    if (y > world.landscape().size_y() - 1) {
        std::ostringstream msg;
        msg << "maximum y-coordinate on landscape is " << world.landscape().size_y() - 1;
        msg << " (0-based indexing), but y-coordinate of " << y << " specified";
        if (item_desc != NULL) {
            msg << " for " << item_desc;
        }
        throw ConfigurationError(msg.str());
    }
    return world.landscape().xy_to_index(x, y);
}

XmlElementType ConfigurationFile::get_child_node(XmlElementType& current_node, const char * node_name, bool required) {
    XmlElementType cnode = current_node.getChildNode(node_name);
    if (cnode.isEmpty() && required) {
        std::ostringstream msg;
        msg << "mandatory element \"" << node_name << "\" is missing from configuration file";
        throw ConfigurationSyntaxError(msg.str());
    }
    return cnode;
}

} // confsys_detail

} // namespace confsys

} // namespace ginkgo
