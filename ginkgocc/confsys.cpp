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
#include "biosys.hpp"
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
    this->config_filepath_ = fpath;
    this->xml_ = XMLNode::openFileHelper(fpath,"ginkgo");
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

void ConfigurationFile::process_world(World& world) {

    XmlElementType world_node = this->get_child_node(this->xml_, "world");
    world.set_label( this->get_attribute<std::string>(world_node, "label", "GinkgoWorld") );
    world.set_random_seed( this->get_attribute<unsigned long>(world_node, "random_seed", time(0)) );
    world.set_generations_to_run( this->get_attribute<GenerationCountType>(world_node, "num_gens") );
    unsigned fitness_dim = this->get_attribute<unsigned>(world_node, "num_fitness_traits");
    if (fitness_dim > MAX_FITNESS_TRAITS) {
        std::ostringstream s;
        s << "maximum number of fitness traits allowed is " << MAX_FITNESS_TRAITS;
        s << ", but requested " << fitness_dim;
        throw ConfigurationError(s.str());
    }
    world.set_num_fitness_traits(fitness_dim);
    world.set_global_selection_strength(this->get_attribute<float>(world_node, "global_selection_strength", DEFAULT_GLOBAL_SELECTION_STRENGTH));
    world.set_allow_multifurcations(this->get_attribute_bool(world_node, "multifurcating_trees", true));
    world.set_produce_final_output(this->get_attribute_bool(world_node, "final_output", false));
    world.set_produce_full_complement_diploid_trees(this->get_attribute_bool(world_node, "full_complement_diploid_trees", false));
    world.generate_landscape( this->get_attribute<CellIndexType>(world_node, "x_range"),
                              this->get_attribute<CellIndexType>(world_node, "y_range") );
    world.set_global_cell_carrying_capacity(this->get_attribute<PopulationCountType>(world_node, "default_cell_carrying_capacity", 0));
    world.set_log_frequency(this->get_attribute<unsigned>(world_node, "log_frequency", 10));
}

void ConfigurationFile::process_biota(World& world) {
    XmlElementType bio_node = this->xml_.getChildNode("world").getChildNode("biota");
    if (bio_node.isEmpty()) {
        throw ConfigurationSyntaxError("biota element is missing from configuration file");
    }

    for (int i = 0; i < bio_node.nChildNode("lineage"); ++i) {
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
    XmlElementType gtf_node = this->get_child_node(lineage_node, "fitnessTraitDefaultGenotype", false);
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
    XmlElementType sw_node = this->get_child_node(lineage_node, "fitnessTraitRelativeSelectionWeights", false);
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
    XmlElementType fsd_node = this->get_child_node(lineage_node, "fitnessTraitInheritanceStdDev", false);
    if (!fsd_node.isEmpty()) {
        std::vector<float> sd = this->get_element_vector<float>(fsd_node);
        if (sd.size() != lineage.get_num_fitness_traits()) {
            std::ostringstream msg;
            msg << "expecting " << lineage.get_num_fitness_traits();
            msg << " fitness inheritance standard deviations, but found ";
            msg << sd.size() << " instead";
            throw ConfigurationError(msg.str());
        }
        lineage.set_fitness_trait_inheritance_sd(sd);
    }

    lineage.set_mean_reproductive_rate(this->get_child_node_scalar<unsigned>(lineage_node, "fecundity", 16));
    lineage.set_movement_capacity(this->get_child_node_scalar<MovementCountType>(lineage_node, "movementCapacity", 1));

    // seed populations
    XmlElementType seed_pops = lineage_node.getChildNode("seedPopulations");
    if (seed_pops.isEmpty()) {
        throw ConfigurationError("no seed populations defined for lineage \"" + lineage_id + "\"");
    }
    for (int i = 0; i < seed_pops.nChildNode("seedPopulation"); ++i) {
        XmlElementType pop_node = seed_pops.getChildNode("seedPopulation", i);
        std::ostringstream item_desc;
        item_desc << "seed population " << i+1 << " for lineage \"" << lineage_id << "\"";
        CellIndexType cell_index = this->get_validated_cell_index(this->get_attribute<CellIndexType>(pop_node, "x"),
                                                                  this->get_attribute<CellIndexType>(pop_node, "y"),
                                                                  world,
                                                                  item_desc.str().c_str());
        PopulationCountType size = this->get_attribute<PopulationCountType>(pop_node, "size");
        PopulationCountType ancestral_pop_size = this->get_child_node_scalar<PopulationCountType>(pop_node, "ancestralPopulationSize");
        GenerationCountType ancestral_generations = this->get_child_node_scalar<GenerationCountType>(pop_node, "ancestralGenerations");
        world.add_seed_population(cell_index, &lineage, size, ancestral_pop_size, ancestral_generations);
    }
}

void ConfigurationFile::process_environments(World& world) {
    XmlElementType environs = this->xml_.getChildNode("world").getChildNode("environments");
    if (!environs.isEmpty()) {
        for (int i = 0; i < environs.nChildNode("environment"); ++i) {
            XmlElementType env_node = environs.getChildNode("environment", i);
            GenerationCountType gen = this->get_attribute<GenerationCountType>(env_node, "gen");
            WorldSettings world_settings;
            for (int j = 0; j < env_node.nChildNode(); ++j) {
                XmlElementType sub_node = env_node.getChildNode(j);
                std::string node_name = sub_node.getName();
                if ( node_name == "carryingCapacity") {
                    world_settings.carrying_capacity = this->get_validated_grid_path<PopulationCountType>(this->get_element_scalar<std::string>(sub_node), world);
                } else if (node_name == "fitnessTraitOptima") {
                    unsigned eidx = this->get_attribute<unsigned>(sub_node, "trait");
                    if (eidx > world.get_num_fitness_traits()) {
                        std::ostringstream msg;
                        msg << "invalid fitness trait index: " << eidx;
                        msg << " (maximum valid index is " << world.get_num_fitness_traits();
                        msg << ", given " << world.get_num_fitness_traits() << " defined fitness traits)";
                        throw ConfigurationError(msg.str());
                    }
                    std::string gridfile = this->get_validated_grid_path<FitnessTraitType>(this->get_element_scalar<std::string>(sub_node), world);
                    world_settings.fitness_trait_optima.insert(std::make_pair(eidx, gridfile));
                } else if (node_name == "movementCosts") {
                    std::string lineage_id = this->get_attribute<std::string>(sub_node, "lineage");
                    if (not world.has_species(lineage_id)) {
                        throw ConfigurationError("movement costs: lineage \"" + lineage_id + "\" not defined");
                    }
                    Species * lineage = world.get_species_ptr(lineage_id);
                    std::string gridfile = this->get_validated_grid_path<MovementCountType>(this->get_element_scalar<std::string>(sub_node), world);
                    world_settings.movement_costs.insert(std::make_pair(lineage, gridfile));
                } else if (node_name == "movementProbablities") {
                    std::string lineage_id = this->get_attribute<std::string>(sub_node, "lineage");
                    if (not world.has_species(lineage_id)) {
                        throw ConfigurationError("movement probabilities: lineage \"" + lineage_id + "\" not defined");
                    }
                    Species * lineage = world.get_species_ptr(lineage_id);
                    std::string gridfile = this->get_validated_grid_path<float>(this->get_element_scalar<std::string>(sub_node), world);
                    world_settings.movement_probabilities.insert(std::make_pair(lineage, gridfile));
                }
            }
            world.add_world_settings(gen, world_settings);
        }
    }
}

void ConfigurationFile::process_dispersals(World& world) {
    XmlElementType dispersals = this->xml_.getChildNode("world").getChildNode("dispersals");
    if (!dispersals.isEmpty()) {
        for (int i = 0; i < dispersals.nChildNode("dispersal"); ++i) {
            XmlElementType disp_node = dispersals.getChildNode("dispersal", i);
            GenerationCountType gen = this->get_attribute<GenerationCountType>(disp_node, "gen");
            DispersalEvent disp_event;
            std::ostringstream item_desc;
            item_desc << "dispersal event " << i+1 << " in generation " << gen;
            disp_event.source = this->get_validated_cell_index(this->get_attribute<CellIndexType>(disp_node, "from_x"),
                    this->get_attribute<CellIndexType>(disp_node, "from_y"),
                    world,
                    item_desc.str().c_str());
            disp_event.destination = this->get_validated_cell_index(this->get_attribute<CellIndexType>(disp_node, "to_x"),
                    this->get_attribute<CellIndexType>(disp_node, "to_y"),
                    world,
                    item_desc.str().c_str());
            disp_event.probability = this->get_child_node_scalar<float>(disp_node, "probability", 1.0);
            std::string lineage_id = this->get_child_node_scalar<std::string>(disp_node, "lineage", "");
            if (lineage_id.size() > 0) {
                if (not world.has_species(lineage_id)) {
                    throw ConfigurationError("dispersal: lineage \"" + lineage_id + "\" not defined");
                }
                disp_event.species_ptr = world.get_species_ptr(lineage_id);
            } else {
                disp_event.species_ptr = NULL;
            }
            world.add_dispersal_event(gen, disp_event);
        }
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

void ConfigurationFile::process_samplings(World& world) {
    XmlElementType samplings = this->xml_.getChildNode("world").getChildNode("samples");
    if (!samplings.isEmpty()) {
        for (int i = 0; i < samplings.nChildNode(); ++i) {
            XmlElementType snode = samplings.getChildNode(i);
            std::string node_name = snode.getName();
            if (node_name == "sample") {

                GenerationCountType gen = this->get_attribute<GenerationCountType>(snode, "gen");
                std::string lineage_id = this->get_attribute<std::string>(snode, "lineage");
                if (not world.has_species(lineage_id)) {
                    throw ConfigurationError("occurrence sample: lineage \"" + lineage_id + "\" not defined");
                }
                world.add_occurrence_sampling(gen, world.get_species_ptr(lineage_id));

                SamplingRegime world_sampling_regime;
                world_sampling_regime.species_ptr = world.get_species_ptr(lineage_id);
                std::string label = this->get_attribute<std::string>(snode, "label", "");
                if (label.size() > 0) {
                    world_sampling_regime.label = label;
                }
                world_sampling_regime.num_organisms_per_cell = this->get_child_node_scalar<PopulationCountType>(snode, "individualsPerCell", 0);
                XmlElementType cells_node = snode.getChildNode("cells");
                if (!cells_node.isEmpty()) {
                    std::ostringstream raw;
                    for (int i = 0; i < cells_node.nText(); ++i) {
                        raw << cells_node.getText(i);
                    }
                    std::string cells_desc = raw.str();
                    if (cells_desc.size() > 0) {
                        std::vector<std::string> cells_vec = textutil::split_on_any(cells_desc, " \r\n\t", 0, false);
                        for (std::vector<std::string>::iterator ci = cells_vec.begin(); ci != cells_vec.end(); ++ci) {
                            std::vector<std::string> xy = textutil::split(*ci, ",", 0, false);
                            if (xy.size() < 2) {
                                throw ConfigurationError("tree sample cell position: missing coordinate");
                            } else if (xy.size() > 2) {
                                throw ConfigurationError("tree sample cell position: too many coordinates");
                            }
                            CellIndexType x = convert::to_scalar<CellIndexType>(xy[0]);
                            CellIndexType y = convert::to_scalar<CellIndexType>(xy[1]);
                            CellIndexType cell_index = this->get_validated_cell_index(x, y, world, "sampling coordinate");
                            world_sampling_regime.cell_indexes.insert(cell_index);
                        }
                    }
                }
                world.add_tree_sampling(gen, world_sampling_regime);
            }
        }
    }
}

void ConfigurationFile::configure(World& world) {
    this->process_world(world);
    this->process_biota(world);
    this->process_environments(world);
    this->process_dispersals(world);
    this->process_samplings(world);
}

} // confsys_detail

} // namespace confsys

} // namespace ginkgo
