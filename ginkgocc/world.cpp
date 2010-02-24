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

#include <iostream>
#include <string>
#include <map>
#include <utility>
#include <ctime>
#include <fstream>
#include <iomanip>

#include "asciigrid.hpp"
#include "world.hpp"
#include "filesys.hpp"
#include "tree.hpp"
#include "convert.hpp"

#if defined(MEMCHECK)
#include "memcheck.hpp"
#endif

using namespace ginkgo;

// singleton instance
World World::instance_;

// constructor
World::World()
    : rng_(RandomNumberGenerator::get_instance()),
      landscape_(*this),
      num_fitness_traits_(1),
      global_selection_strength_(DEFAULT_GLOBAL_SELECTION_STRENGTH),
      generations_to_run_(0),
      current_generation_(0),
      log_frequency_(10),
      is_log_to_screen_(true),
      allow_multifurcations_(true),
      is_produce_final_output_(true),
      is_produce_full_complement_diploid_trees_(false) {
    this->current_generation_ = 0;
}

World::~World() {}

// --- initialization and set up ---

// Creates a new landscape.
void World::generate_landscape(CellIndexType size_x, CellIndexType size_y) {
    if ((size_x * size_y) > MAX_LANDSCAPE_SIZE) {
//        throw LandscapeSizeError("Maximum landscape size is " +  UCHAR_MAX  + "cells, but requested " + size_x + " x " + size_y + " = "  + (size_x*size_y) + " cells")
        std::ostringstream msg;
        msg << "Maximum landscape size allowed is " <<  MAX_LANDSCAPE_SIZE << ", but requested landscape has " << size_x << " x " << size_y << " = "  << (size_x*size_y) << " cells";
        throw LandscapeSizeError(msg.str());
    }
    this->landscape_.generate(size_x, size_y, this->num_fitness_traits_);
}

// Adds a new species definition to this world.
Species& World::new_species(const std::string& label) {
    Species* sp = this->species_registry_.new_species(label);
    sp->set_num_fitness_traits(this->num_fitness_traits_);
    sp->set_global_selection_strength(this->global_selection_strength_);
    std::vector<MovementCountType> default_movement_costs(this->landscape_.size(), 1);
    sp->set_movement_costs(default_movement_costs);
    std::vector<float> default_movement_probabilities(this->landscape_.size(), 1.0);
    sp->set_movement_probabilities(default_movement_probabilities);
    return *sp;
}

// Populates the cell cell_index with organisms of the given species.
void World::generate_seed_population(CellIndexType cell_index,
        Species * species_ptr,
        PopulationCountType pop_size,
        PopulationCountType ancestral_pop_size,
        GenerationCountType ancestral_generations) {
    std::ostringstream pre_msg;
    pre_msg << "[Generation " << this->current_generation_ << "] ";
    pre_msg << "Bootstrapping seed population of species " << species_ptr->get_label() << ": ";
    pre_msg << ancestral_pop_size << " individuals for " << ancestral_generations << " generations.";
    this->log_info(pre_msg.str());

    this->landscape_.at(cell_index).generate_new_population(species_ptr,
            pop_size,
            ancestral_pop_size,
            ancestral_generations);

    std::ostringstream post_msg;
    post_msg << "[Generation " << this->current_generation_ << "] ";
    post_msg << "Seeding population ";
    post_msg << "of species " << species_ptr->get_label()  << " ";
    post_msg << "in (" << this->landscape_.index_to_x(cell_index) <<  "," << this->landscape_.index_to_y(cell_index) << "): ";

    PopulationCountType num_females = 0;
    PopulationCountType num_males = 0;
    this->landscape_.at(cell_index).num_organisms(species_ptr, num_females, num_males);
    post_msg << num_females << " females and " << num_males << " males ";
    post_msg << "drawn from an ancestral population of " << ancestral_pop_size << " ";
    post_msg << "after " << ancestral_generations << " generations.";
    this->log_info(post_msg.str());
}

// --- event handlers ---

void World::add_world_settings(GenerationCountType generation, const WorldSettings& world_settings) {
    this->world_settings_[generation] = world_settings;
}

void World::add_dispersal_event(GenerationCountType generation, const DispersalEvent& dispersal_event) {
    this->log_error("Stochastic dispersal currently disabled");
    exit(1);
    this->dispersal_events_.insert(std::make_pair(generation, dispersal_event));
}


void World::add_tree_sampling(GenerationCountType generation, const SamplingRegime& sampling_regime) {
    this->tree_samples_.insert(std::make_pair(generation, sampling_regime));
}

void World::add_occurrence_sampling(GenerationCountType generation, Species * species_ptr) {
    this->occurrence_samples_.insert(std::make_pair(generation, species_ptr));
}

void World::add_seed_population(CellIndexType cell_index,
        Species * species_ptr,
        PopulationCountType pop_size,
        PopulationCountType ancestral_pop_size,
        GenerationCountType ancestral_generations) {
    this->seed_populations_.push_back( SeedPopulation(cell_index, species_ptr, pop_size, ancestral_pop_size, ancestral_generations) );
}


// --- simulation cycles ---

void World::cycle() {

// Results in inflated population: the migrants get distributed after the
// competition phase, resulting in a artificially (>> carrying capacity)
// boosted population when entering the next generation's reproduction phase.
// This leads to a standing population at the end of each generation sometimes
// an order or more of magnitude above the carrying capacity of the cell.
//     for (CellIndexType i = this->landscape_.size()-1; i >= 0; --i) {
//         this->landscape_[i].reproduction();
//         this->landscape_[i].migration();
//         this->landscape_[i].survival();
//         this->landscape_[i].competition();
//     }
//     this->landscape_.process_migrants();

    if ( this->current_generation_ % this->log_frequency_ == 0) {
        std::ostringstream gen;
        gen << "Generation " << this->current_generation_ << " life-cycle running.";
        this->log_info(gen.str());
    }
//     this->log_detail("Reproduction/migration phase.");
    for (CellIndexType i = 0; i < this->landscape_.size(); ++i) {
        this->landscape_[i].reproduction();
        this->landscape_[i].migration();
    }
//     this->log_detail("Processing migrants.");
    this->landscape_.process_migrants();
//     this->log_detail("Survival/competition phase.");
    for (CellIndexType i = 0; i < this->landscape_.size(); ++i) {
        this->landscape_[i].survival();
        this->landscape_[i].competition();
    }
//     this->log_detail("Generation life-cycle completed.");
    ++this->current_generation_;
}

void World::run() {
    this->open_logs();
    this->log_configuration();

#if defined(DEBUG)
    this->log_info("RUNNING DEBUG-MODE BUILD.");
#endif

#if defined(MEMCHECK)
    this->log_info("UNRELEASED NODE MEMORY WILL BE LOGGED.");
#endif

    this->log_info("Starting simulation.");

    // startup
    for (std::vector<SeedPopulation>::iterator spi = this->seed_populations_.begin();
            spi != this->seed_populations_.end();
            ++spi) {
        SeedPopulation& sp = *spi;
        this->generate_seed_population(sp.cell_index,
            sp.species_ptr,
            sp.pop_size,
            sp.ancestral_pop_size,
            sp.ancestral_generations);
    }

    while (this->current_generation_ <= this->generations_to_run_) {

        // clear organism labels
        for (SpeciesRegistry::iterator spi = this->species_registry_.begin();
                spi != this->species_registry_.end();
                ++spi) {
            (*spi)->clear_organism_labels();
        }

        // clear output filename stems
        this->output_filenames_.clear();

        // save occurrence data requested in this generation
        this->process_occurrence_samplings();

        // build trees requested in this generation
        this->process_tree_samplings();

        // process world changes
        this->process_world_settings();

        // process dispersal events
        // this->process_dispersal_events();

        // run the life cycle
        this->cycle();
    }

    if (this->is_produce_final_output_) {
        this->log_info("Saving final set of occurrences and trees for all species.");
        std::set<CellIndexType> cell_indexes;
        for (CellIndexType i = 0; i < this->landscape_.size(); ++i) {
            cell_indexes.insert(i);
        }
        for (SpeciesRegistry::iterator spi = this->species_registry_.begin();
                spi != this->species_registry_.end();
                ++spi) {
            this->save_occurrences(*spi);
            this->save_trees(*spi, 0, cell_indexes, "COMPLETE");
        }
    }

#if defined(MEMCHECK)
    run_final_cleanup_and_memory_check();
#endif

    this->log_info("Ending simulation.");
}

void World::process_world_settings() {
    std::map<GenerationCountType, WorldSettings>::iterator wi = this->world_settings_.find(this->current_generation_);
    if (wi == this->world_settings_.end()) {
        return;
    }
    if (wi->second.carrying_capacity.size() != 0) {
        this->log_info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Setting carrying capacity: \"" + wi->second.carrying_capacity + "\".");
        asciigrid::AsciiGrid<PopulationCountType> grid(wi->second.carrying_capacity);
        this->landscape_.set_carrying_capacities(grid.get_cell_values());
    }
    if (wi->second.fitness_trait_optima.size() != 0) {
        for (std::map<unsigned, std::string>::iterator ei = wi->second.fitness_trait_optima.begin();
                 ei != wi->second.fitness_trait_optima.end();
                 ++ei) {
            std::ostringstream msg;
            msg << "[Generation " << this->current_generation_ << "] Setting optima for fitness trait " <<  ei->first <<  ": \"" <<  ei->second <<  "\"";
            this->log_info(msg.str());
            asciigrid::AsciiGrid<FitnessTraitType> grid(ei->second);
            this->landscape_.set_environment(ei->first, grid.get_cell_values());
        }
    }
    if (wi->second.movement_costs.size() != 0) {
        for (std::map<Species *, std::string>::iterator mi = wi->second.movement_costs.begin();
                 mi != wi->second.movement_costs.end();
                 ++mi) {
            std::ostringstream msg;
            msg << "[Generation " << this->current_generation_ << "] Setting movement costs for species " <<  mi->first->get_label() <<  ": \"" <<  mi->second <<  "\"";
            this->log_info(msg.str());
            asciigrid::AsciiGrid<MovementCountType> grid(mi->second);
            mi->first->set_movement_costs(grid.get_cell_values());
        }
    }
    if (wi->second.movement_probabilities.size() != 0) {
        for (std::map<Species *, std::string>::iterator mi = wi->second.movement_probabilities.begin();
                 mi != wi->second.movement_probabilities.end();
                 ++mi) {
            std::ostringstream msg;
            msg << "[Generation " << this->current_generation_ << "] Setting movement probabilities for species " <<  mi->first->get_label() <<  ": \"" <<  mi->second <<  "\"";
            this->log_info(msg.str());
            asciigrid::AsciiGrid<float> grid(mi->second);
            mi->first->set_movement_probabilities(grid.get_cell_values());
        }
    }
}

void World::process_dispersal_events() {
    this->log_error("Stochastic dispersal currently disabled");
    exit(1);
//    typedef std::multimap<GenerationCountType, DispersalEvent> gen_disp_t;
//    typedef std::pair<gen_disp_t::iterator, gen_disp_t::iterator> gen_disp_iter_pair_t;
//    gen_disp_iter_pair_t this_gen_dispersals = this->dispersal_events_.equal_range(this->current_generation_);
//    for (gen_disp_t::iterator di = this_gen_dispersals.first; di != this_gen_dispersals.second; ++di) {
//        DispersalEvent& de = di->second;
//        std::ostringstream msg;
//        msg << "[Generation " << this->current_generation_ << "] Dispersal";
//        if (de.species_ptr != NULL) {
//            msg << " of " << de.species_ptr->get_label();
//        }
//        msg << " from (" << this->landscape_.index_to_x(de.source) << "," << this->landscape_.index_to_y(de.source) << ")";
//        msg << " to (" << this->landscape_.index_to_x(de.destination) << "," << this->landscape_.index_to_y(de.destination) << ")";
//        msg << " with probability " << de.probability << ": ";
//
//        std::vector<const Organism *> organisms;
//        this->landscape_[de.source].sample_organisms(de.species_ptr, organisms, de.num_organisms);
//        this->landscape_.clear_migrants();
//        unsigned long num_males = 0;
//        unsigned long num_females = 0;
//        for (std::vector<const Organism *>::iterator oi = organisms.begin();  oi != organisms.end(); ++oi) {
//            if (this->rng().uniform_01() > de.probability) {
//                Organism* og = const_cast<Organism*>(*oi);
//                if (og->is_male()) {
//                    ++num_males;
//                } else {
//                    ++num_females;
//                }
//                this->landscape_.add_migrant(*og, de.destination);
//                og->set_expired();
//            }
//        }
//        msg << " " << num_females << " females and " << num_males << " males.";
//        this->landscape_[de.source].purge_expired_organisms();
//        this->landscape_.process_migrants();
//        this->log_info(msg.str());
//    }
}

void World::process_tree_samplings() {
    typedef std::multimap<GenerationCountType, SamplingRegime> gen_sample_t;
    typedef std::pair<gen_sample_t::iterator, gen_sample_t::iterator> gen_sample_iter_pair_t;
    gen_sample_iter_pair_t this_gen_samples = this->tree_samples_.equal_range(this->current_generation_);
    for (gen_sample_t::iterator i = this_gen_samples.first; i != this_gen_samples.second; ++i) {
        this->save_trees(i->second.species_ptr, i->second.num_organisms_per_cell, i->second.cell_indexes, i->second.label);
    }
}

void World::process_occurrence_samplings() {
    typedef std::multimap<GenerationCountType, Species *> gen_sample_t;
    typedef std::pair<gen_sample_t::iterator, gen_sample_t::iterator> gen_sample_iter_pair_t;
    gen_sample_iter_pair_t this_gen_samples = this->occurrence_samples_.equal_range(this->current_generation_);
    for (gen_sample_t::iterator i = this_gen_samples.first; i != this_gen_samples.second; ++i) {
        this->save_occurrences(i->second);
    }
}

#if defined(MEMCHECK)
    // --- memory check ---
    void World::run_final_cleanup_and_memory_check() {
        for (CellIndexType i = 0; i < this->landscape_.size(); ++i) {
            this->landscape_.at(i).organisms().clear();
        }
        UNRELEASED_NODES_LOG.open((this->get_output_filename_stem() + ".unreleased.log").c_str());
        UNRELEASED_NODES_LOG << "### This file logs all zombie GenealogyNode objects as they get destroyed ###" << std::endl;
    }
#endif

// --- logging and output ---

void World::save_occurrences(Species * species_ptr ) {
    this->log_info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Saving occurrence data for species " + species_ptr->get_label() + ".");
    std::vector<PopulationCountType> counts;
    this->landscape_.count_organisms(species_ptr, counts);
    std::ofstream occs;
    this->open_ofstream(occs,
        this->compose_output_filename(species_ptr->get_label(), "occurrences", "asc"));
    asciigrid::write_grid(counts, this->landscape_.size_x(), this->landscape_.size_y(), occs);
}

void World::write_nexus_header(Species * sp_ptr,
        const std::vector<const Organism *>& organisms,
        std::ostream& out,
        bool add_diploid_allele_extensions) {
    out << "#NEXUS\n\n";
    out << "BEGIN TAXA;\n";
    if (add_diploid_allele_extensions) {
        out << "    DIMENSIONS NTAX=" << organisms.size()*2 << ";\n";
    } else {
        out << "    DIMENSIONS NTAX=" << organisms.size() << ";\n";
    }
    out << "    TAXLABELS\n";
    for (std::vector<const Organism *>::const_iterator oi = organisms.begin();
            oi != organisms.end();
            ++oi) {
        if (add_diploid_allele_extensions) {
            out << "        " << sp_ptr->get_organism_label(*oi) << "_a1" << "\n";
            out << "        " << sp_ptr->get_organism_label(*oi) << "_a2" << "\n";
        } else {
            out << "        " << sp_ptr->get_organism_label(*oi) << "\n";
        }
    }
    out << "    ;\n";
    out << "END;\n\n";
}

void World::write_traits(Species * sp_ptr,
                const std::vector<const Organism *>& organisms,
                std::ostream& out) {
    // write data
    this->write_nexus_header(sp_ptr, organisms, out);
    out << "BEGIN CHARACTERS;\n";
    unsigned int nchar = this->num_fitness_traits_ + 1;
    out << "    DIMENSIONS NCHAR=" << nchar << ";\n";
    out << "    FORMAT DATATYPE=CONTINUOUS ITEMS=(STATES);\n";
    out << "    CHARLABELS";
    for (unsigned int i=0; i < this->num_fitness_traits_; ++i) {
        out << " Trait_" << std::setw(2) << std::setfill('0') << i;
    }
    out << " Fitness;\n";

    // the things I do for OCD-driven tight formatting ...
    long max_label_len = 0;
    for (std::vector<const Organism *>::const_iterator oi = organisms.begin();
            oi != organisms.end();
            ++oi) {
        long sz = sp_ptr->get_organism_label(*oi).size();
        if (sz > max_label_len) {
            max_label_len = sz;
        }
    }
    out << "        " << std::setw(max_label_len) << std::setfill(' ') << "" << "    ";
    for (unsigned int i=0; i < this->num_fitness_traits_; ++i) {
        out << " " << std::setw(9) << std::setfill(' ') << "[Trait_" << std::setw(2) << std::setfill('0') << i << "]";
    }
    out << " " << std::setw(12) << std::setfill(' ') << "[Fitness]";
    out << std::endl;

    out << "    MATRIX\n";

    for (std::vector<const Organism *>::const_iterator oi = organisms.begin();
            oi != organisms.end();
            ++oi) {
        out << "        " << std::left << std::setw(max_label_len) << sp_ptr->get_organism_label(*oi) << "    ";
        for (unsigned int i=0; i < this->num_fitness_traits_; ++i) {
            out << " " << std::right << std::setw(12) << std::setfill(' ') << (**oi).get_fitness_trait_genotype(i);
        }
        out << " " << std::setw(12) << std::setfill(' ') << (**oi).get_fitness();
        out << "\n";
    }
    out << ";\n";
    out << "END;\n\n";
}

void World::write_haploid_tree(Species * sp_ptr,
                const std::vector<const Organism *>& organisms,
                std::ostream& out) {
    Tree tree(&this->landscape_, this->allow_multifurcations_);
    try {
        // build tree
        for (std::vector<const Organism *>::const_iterator oi = organisms.begin();
                oi != organisms.end();
                ++oi) {
            tree.add_leaf((*oi)->get_haploid_node(), &sp_ptr->get_organism_label(*oi));
        }

        // write tree
        this->write_nexus_header(sp_ptr, organisms, out);
        out << "BEGIN TREES;\n";
        out << "    TREE HaploidLocus = [&R] ";
        tree.write_newick_tree(out);
        out << ";\n";
        out << "END;\n\n";
    } catch (const TreeStructureMultipleRootError& e) {
        this->log_tree_multiple_root_error(sp_ptr->get_label(), organisms.size());
    }
}

void World::write_diploid2_trees(Species * sp_ptr,
                const std::vector<const Organism *>& organisms,
                std::ostream& out) {
    try {
        this->write_nexus_header(sp_ptr, organisms, out, true);
        out << std::setfill(' '); // reset
        out << "BEGIN TREES;\n";
        for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_loci; ++i) {
            Tree tree(&this->landscape_, this->allow_multifurcations_);
            std::string allele1;
            std::string allele2;
            for (std::vector<const Organism *>::const_iterator oi = organisms.begin();
                    oi != organisms.end();
                    ++oi) {
                allele1 = sp_ptr->get_organism_label(*oi) + "_a1";
                allele2 = sp_ptr->get_organism_label(*oi) + "_a2";
                tree.add_leaf((*oi)->get_diploid_node1(i), &allele1);
                tree.add_leaf((*oi)->get_diploid_node2(i), &allele2);
            }
            out << "    TREE DiploidLocus" << std::setw(2) << std::setfill('0') << i+1 << " = [&R] ";
            tree.write_newick_tree(out);
            out << ";\n";
        }
        out << "END;\n\n";
        out << std::setfill(' '); // reset
    } catch (const TreeStructureMultipleRootError& e) {
        this->log_tree_multiple_root_error(sp_ptr->get_label(), organisms.size());
    }
}

void World::write_diploid1_trees(Species * sp_ptr,
        const std::vector<const Organism *>& organisms,
        std::ostream& out) {
    try {
        this->write_nexus_header(sp_ptr, organisms, out);
        out << std::setfill(' '); // reset
        out << "BEGIN TREES;\n";
        for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_loci; ++i) {
            Tree tree(&this->landscape_, this->allow_multifurcations_);
            for (std::vector<const Organism *>::const_iterator oi = organisms.begin();
                    oi != organisms.end();
                    ++oi) {
                tree.add_leaf((*oi)->get_diploid_random_node(i, this->rng_), &sp_ptr->get_organism_label(*oi));
            }
            out << "    TREE DiploidLocus" << std::setw(2) << std::setfill('0') << i+1 << " = [&R] ";
            tree.write_newick_tree(out);
            out << ";\n";
        }
        out << "END;\n\n";
        out << std::setfill(' '); // reset
    } catch (const TreeStructureMultipleRootError& e) {
        this->log_tree_multiple_root_error(sp_ptr->get_label(), organisms.size());
    }
}

void World::log_tree_multiple_root_error(const std::string& species_label, unsigned long num_taxa) {
    std::ostringstream msg;
    msg << "tree for sample of organisms of species " << species_label;
    msg << " in generation " << this->current_generation_;
    msg << " could not be built due to multiple root nodes";
    msg << " (organism sample size = " << num_taxa << ").";
    this->log_error(msg.str());
}

void World::save_trees(Species * sp_ptr,
                PopulationCountType num_organisms_per_cell,
                const std::set<CellIndexType>& cell_indexes,
                const std::string& label) {

    std::ostringstream msg;
    msg << "[Generation " << this->current_generation_ << "] ";
    msg << "Sampling organisms of species " << sp_ptr->get_label() << " (";
    if (num_organisms_per_cell == 0) {
        msg << "all individuals per cell,";
    } else {
        msg << num_organisms_per_cell << " individuals per cell,";
    }
    if (cell_indexes.size() == 0) {
        msg << " all cells).";
    } else {
        msg << " from " << cell_indexes.size() << " cells).";
    }
    this->log_info(msg.str());

    std::vector<const Organism *> organisms;
    this->landscape_.sample_organisms(sp_ptr, num_organisms_per_cell, cell_indexes, organisms);
    std::string num_samples = convert::to_scalar<std::string>(organisms.size());

    if (organisms.size() == 0) {
        this->log_error("no organisms found in sample: aborting tree building");
        return;
    }

    this->log_info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Building single allele diploid loci trees for " + num_samples + " organisms (" + num_samples + " leaves per tree).");
    std::ofstream combined_trees;
    this->open_ofstream(combined_trees,
        this->compose_output_filename(sp_ptr->get_label(), label, "diploid1.tre"));
    this->write_diploid1_trees(sp_ptr, organisms, combined_trees);

    this->log_info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Building haploid locus tree for " + num_samples + " organisms (" + num_samples + " leaves per tree).");
    std::ofstream haploid_trees;
    this->open_ofstream(haploid_trees,
        this->compose_output_filename(sp_ptr->get_label(), label, "haploid.tre"));
    this->write_haploid_tree(sp_ptr, organisms, haploid_trees);

    if (this->is_produce_full_complement_diploid_trees_) {
        this->log_info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Building full complement diploid loci trees for " + num_samples + " organisms (" + convert::to_scalar<std::string>(organisms.size()*2) + " leaves per tree).");
        std::ofstream diploid_trees;
        this->open_ofstream(diploid_trees,
            this->compose_output_filename(sp_ptr->get_label(), label, "diploid2.tre"));
        this->write_diploid2_trees(sp_ptr, organisms, diploid_trees);
    }

    this->log_info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Building traits matrix for " + num_samples + " organisms (" + num_samples + " leaves per tree).");
    std::ofstream traits;
    this->open_ofstream(traits,
        this->compose_output_filename(sp_ptr->get_label(), label, "traits.nex"));
    this->write_traits(sp_ptr, organisms, traits);

}

void World::open_ofstream(std::ofstream& out, const std::string& fpath) {
    std::string full_fpath = filesys::compose_path(this->output_dir_, fpath);
    out.open(full_fpath.c_str());
    if (not out) {
        throw WorldIOError("cannot open file \"" + full_fpath + "\" for output");
    }
}

std::string World::get_output_filename_stem() {
    if (this->output_filename_stem_.size() == 0) {
        if (this->label_.size() == 0) {
            this->output_filename_stem_ += "ginkgorun";
        } else {
            this->output_filename_stem_ += this->label_;
        }
        if (this->replicate_id_.size() > 0) {
            this->output_filename_stem_ += this->replicate_id_;
        }
    }
    return this->output_filename_stem_;
}

std::string World::compose_output_filename(const std::string& species_label,
        const std::string& additional,
        const std::string& extension) {
    std::ostringstream f;
    f << this->get_output_filename_stem();
    f << "_G" << std::setw(8) << std::setfill('0') <<  this->current_generation_;
    f << "_" << species_label;
    if (additional.size() > 0) {
        f << "_" << additional;
    }
    unsigned index = 0;
    std::string candidate_name = f.str() + "." + extension;
    std::set<std::string>::iterator fname_found = this->output_filenames_.find(candidate_name);
    while (fname_found != this->output_filenames_.end()) {
        ++index;
        candidate_name = f.str() + "-" + convert::to_scalar<std::string>(index) + "." + extension;
        fname_found = this->output_filenames_.find(candidate_name);
    }
    this->output_filenames_.insert(candidate_name);
    return candidate_name;
}


void World::log_configuration() {
    std::ofstream out;
    this->open_ofstream(out, this->get_output_filename_stem() + ".conf.log");
    out <<  "GINKGO CONFIGURATION LOG " << this->get_timestamp() << std::endl;

    out << std::endl;
    out << "*** LOGGING ***" << std::endl;
    out << "Configuration log: " << this->get_output_filename_stem() + ".conf.log" << std::endl;
    out << "Output log: " << this->get_output_filename_stem() + ".out.log" << std::endl;
    out << "Error log: " << this->get_output_filename_stem() + ".err.log" << std::endl;

#if defined(MEMCHECK)
    out << "!!! RUNNING MEMORY CHECKS !!!" << std::endl;
    out << "Memory log: " << this->get_output_filename_stem() + ".mem.log" << std::endl;
#endif

#if defined(MINI)
    out << "!!! RUNNING MINI MOD !!!" << std::endl;
    out << "Landscape size limited to " << MAX_LANDSCAPE_SIZE << " cells" << std::endl;
#endif

    out << std::endl;
    out << "*** WORLD ***" << std::endl;
    out << "Label: " << this->label_ << std::endl;
    out << "Random seed: " << this->rng_.get_seed() << std::endl;
    out << "Maximum landscape size: " << MAX_LANDSCAPE_SIZE << " cells" << std::endl;
    out << "Number of fitness traits: " << this->num_fitness_traits_ << std::endl;
    out << "Generations to run: " << this->generations_to_run_ << std::endl;
    out << "Output directory: " << this->output_dir_ << std::endl;
    out << "Replicate ID: " << this->replicate_id_ << std::endl;
    out << "Log frequency: " << this->log_frequency_ << std::endl;
    out << "Allow multifurcating trees: ";
    if (this->allow_multifurcations_) {
        out << "yes" << std::endl;
    } else {
        out << "no" << std::endl;
    }
    out << "Full complement diploid trees: ";
    if (this->is_produce_full_complement_diploid_trees_) {
        out << "yes" << std::endl;
    } else {
        out << "no" << std::endl;
    }
    out << "Final output: ";
    if (this->is_produce_final_output_) {
        out << "yes" << std::endl;
    } else {
        out << "no" << std::endl;
    }

    out << std::endl;
    out << "*** LANDSCAPE ***" << std::endl;
    out << "Rows (X-dimension): " << this->landscape_.size_x() << std::endl;
    out << "Columns (Y-dimension): " << this->landscape_.size_y() << std::endl;
    out << "Minumum cell index: 0" << std::endl;
    out << "Maximum cell index: " << this->landscape_.size()-1 << std::endl;
    out << "Landscape cell index limit: " << MAX_LANDSCAPE_SIZE << std::endl;
    if (this->landscape_.get_origin_upper_left()) {
        out << "Landscape grid origin: upper-left" << std::endl;
    } else {
        out << "Landscape grid origin: lower-left" << std::endl;
    }

    out << std::endl;

    out << "Cell indexes: " << std::endl;
    this->landscape_.debug_dump_cell_indexes(out);

    out << "Cell coordinates: " << std::endl;
    this->landscape_.debug_dump_cell_xy(out);

    out << "Default carrying capacity: " << std::endl;
    this->landscape_.debug_dump_carrying_capacity(out);

    out << std::endl;
    out << "Global selection strength: " << this->global_selection_strength_ << std::endl;
    if (this->global_selection_strength_ <= 0) {
        out << "WARNING: GLOBAL SELECTION STRENGTH IS 0: ENVIRONMENT HAS NO EFFECT!" << std::endl;
    }

    out << std::endl;
    out << "*** LINEAGES ***" << std::endl;
    out << "(" << this->species_registry_.size() << " lineages specified)" << std::endl;
    unsigned i = 0;
    for (SpeciesRegistry::iterator spi = this->species_registry_.begin(); spi != this->species_registry_.end(); ++spi) {
        i += 1;
        Species& lineage = **spi;
        out << std::endl;
        out << "\"" << lineage.get_label() << "\"" << std::endl;
        out << "     Number of fitness traits: " << lineage.get_num_fitness_traits() << std::endl;
        out << "    Global selection strength: " << this->global_selection_strength_ << std::endl;
        out << " Normalized selection weights: ";
        std::vector<float> sw = lineage.get_selection_weights();
        for (std::vector<float>::iterator swi = sw.begin(); swi != sw.end(); ++swi) {
            if ((swi - sw.begin()) > 0) {
                out << " ";
            }
            out << *swi;
        }
        out << std::endl;
        out << "     Default heritable fitness trait values: ";
        std::vector<FitnessTraitType> g = lineage.get_default_fitness_trait_genotypes();
        for (std::vector<FitnessTraitType>::iterator gi = g.begin(); gi != g.end(); ++gi) {
            if ((gi - g.begin()) > 0) {
                out << " ";
            }
            out << *gi;
        }
        out << std::endl;
        out << "     Fecundity: " <<  lineage.get_mean_reproductive_rate() << std::endl;
        if (lineage.is_fixed_movement_capacity()) {
            out << "     Movement capacity (fixed): " << lineage.get_movement_capacity() << std::endl;
        } else {
            out << "     Movement capacity (variable): ";
            const std::vector<float>& mc_probs = lineage.movement_capacity_probabilities();
            for (unsigned int i = 0; i < mc_probs.size(); ++i) {
                if (i > 0) {
                    out << ", ";
                }
                out << i << " (p=" << mc_probs[i] << ")";
            }
        }
    }

    out << std::endl;
    out << "*** ENVIRONMENTS ***" << std::endl;
    out << "(" << this->world_settings_.size() << " generations with environmental changes specified)" << std::endl;

    for (std::map<GenerationCountType, WorldSettings>::iterator wi = this->world_settings_.begin(); wi != this->world_settings_.end(); ++wi) {
        out << "\n[ENVIRONMENT: GENERATION " << wi->first << "]" << std::endl;
        if (wi->second.carrying_capacity.size() != 0) {
            out << "    Carrying-capacity: \"" << wi->second.carrying_capacity << "\"" << std::endl;
        }
        if (wi->second.fitness_trait_optima.size() != 0) {
            for (std::map<unsigned, std::string>::iterator ei = wi->second.fitness_trait_optima.begin();
                     ei != wi->second.fitness_trait_optima.end();
                     ++ei) {
                out << "    Fitness trait #" << (ei->first + 1) << " optima: \"" << ei->second << "\"" << std::endl;
            }
        }
        if (wi->second.movement_costs.size() != 0) {
            for (std::map<Species *, std::string>::iterator mi = wi->second.movement_costs.begin();
                     mi != wi->second.movement_costs.end();
                     ++mi) {
                out << "    Movement costs for lineage \"" << mi->first->get_label() << "\": \"" << mi->second << "\"" << std::endl;
            }
        }
        if (wi->second.movement_probabilities.size() != 0) {
            for (std::map<Species *, std::string>::iterator mi = wi->second.movement_probabilities.begin();
                     mi != wi->second.movement_probabilities.end();
                     ++mi) {
                out << "    Movement probablities for lineage \"" << mi->first->get_label() << "\": \"" << mi->second << "\"" << std::endl;
            }
        }
    }

    out << std::endl;
    out << "*** DISPERSALS ***" << std::endl;
    out << "(" << this->dispersal_events_.size() << " dispersal events specified)" << std::endl;
    out << std::endl;
    for (std::map<GenerationCountType, DispersalEvent>::iterator di = this->dispersal_events_.begin(); di != this->dispersal_events_.end(); ++di) {
        DispersalEvent& de = di->second;
        out << "    GENERATION " << di->first << ": " << "Dispersal";
        if (de.species_ptr != NULL) {
            out << " of " << de.species_ptr->get_label();
        }
        out << " from (" << this->landscape_.index_to_x(de.source) << "," << this->landscape_.index_to_y(de.source) << ")";
        out << " to (" << this->landscape_.index_to_x(de.destination) << "," << this->landscape_.index_to_y(de.destination) << ")";
        out << " with probability " << de.probability << "." << std::endl;
    }

    out << std::endl;
    out << "*** OCCURRENCE SAMPLES ***";
    out << std::endl;
    for (std::multimap<GenerationCountType, Species *>::iterator oci = this->occurrence_samples_.begin();
            oci != this->occurrence_samples_.end();
            ++oci) {
        out << "    GENERATION " << oci->first << ": " << oci->second->get_label() << std::endl;
    }

    out << std::endl;
    out << "*** TREE SAMPLES ***";
    out << std::endl;
    for (std::multimap<GenerationCountType, SamplingRegime>::iterator tci = this->tree_samples_.begin();
            tci != this->tree_samples_.end();
            ++tci) {
        SamplingRegime& sr = tci->second;
        out << "    GENERATION " << tci->first << ":";
        if (sr.label.size() > 0) {
            out << " [sample \"" << sr.label << "\"]";
        }
        out << " Lineage \"" << sr.species_ptr->get_label() << "\"";
        if (sr.num_organisms_per_cell == 0) {
            out << ", all individuals";
        } else {
            out << ", " << sr.num_organisms_per_cell << " individuals";
        }
        if (sr.cell_indexes.size() == 0) {
            out << " from each cell, all cells";
        } else {
            out << " from each cell, from following cells: ";
            for (std::set<CellIndexType>::iterator ci = sr.cell_indexes.begin(); ci != sr.cell_indexes.end(); ++ci) {
                out << *ci << " (" << this->landscape_.index_to_x(*ci) << "," << this->landscape_.index_to_y(*ci) << ") ";
            }
        }
        out << std::endl;
    }
    out.close();
}

void World::open_logs() {
    if (not this->infos_.is_open()) {
        this->open_ofstream(this->infos_, this->get_output_filename_stem() + ".out.log");
    }
    if (not this->errs_.is_open()) {
        this->open_ofstream(this->errs_, this->get_output_filename_stem() + ".err.log");
    }
}

void World::close_logs() {
    if (this->infos_.is_open()) {
        this->infos_.close();
    }
    if (this->errs_.is_open()) {
        this->errs_.close();
    }
}

std::string World::get_timestamp() {
    time_t rawtime;
    time ( &rawtime );
    struct tm * timeinfo = localtime ( &rawtime );
    char buffer[80];
    strftime (buffer,80,"%Y-%m-%d %H:%M:%S",timeinfo);
    return "[" + std::string(buffer) + "]";
}

void World::log_detail(const std::string& message) {
    assert(this->infos_);
    this->infos_ << this->get_timestamp();
    this->infos_ << " [Generation ";
    this->infos_ << this->current_generation_;
    this->infos_ << "] ";
    this->infos_ << message << std::endl;
}

void World::log_info(const std::string& message) {
    assert(this->infos_);
    std::string ts = this->get_timestamp();
    if (this->is_log_to_screen_) {
        std::cout << ts << " " << message << std::endl;
    }
    this->infos_ << ts << " " << message << std::endl;
}

void World::log_error(const std::string& message) {
    assert(this->errs_);
    assert(this->infos_);
    std::string ts = this->get_timestamp();
    if (this->is_log_to_screen_) {
        std::cerr << ts << " ERROR: " << message << std::endl;
    }
    this->errs_ << ts << " ERROR: " << message << std::endl;
    this->infos_ << ts << " ERROR: " << message << std::endl;
}


