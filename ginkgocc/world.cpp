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
#include <string>
#include <map>
#include <utility>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "asciigrid.hpp"
#include "world.hpp"
#include "filesys.hpp"
#include "tree.hpp"
#include "convert.hpp"
#include "randgen.hpp"

#if defined(MEMCHECK)
#include "memcheck.hpp"
#endif

using namespace ginkgo;

///////////////////////////////////////////////////////////////////////////////
// RecurringAction

RecurringAction::RecurringAction(GenerationCountType start_gen, GenerationCountType end_gen)
        : start_gen_(start_gen),
          end_gen_(end_gen) {
}

bool RecurringAction::is_active(GenerationCountType current_gen) {
    return (current_gen >= this->start_gen_) && (current_gen <= this->end_gen_);
}

bool RecurringAction::is_expired(GenerationCountType current_gen) {
    return current_gen > this->end_gen_;
}

// RecurringAction
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// JumpDispersalRegime

JumpDispersalRegime::JumpDispersalRegime(GenerationCountType start_gen,
        GenerationCountType end_gen,
        Species * sp_ptr,
        float probability,
        CellIndexType src_cell,
        CellIndexType dest_cell)
    : RecurringAction(start_gen, end_gen),
      species_ptr_(sp_ptr),
      probability_(probability),
      src_cell_(src_cell),
      dest_cell_(dest_cell) {
}

///////////////////////////////////////////////////////////////////////////////
// MigrationTrackingRegime

MigrationTrackingRegime::MigrationTrackingRegime(
        const Landscape& landscape,
        GenerationCountType start_gen,
        GenerationCountType end_gen,
        Species * sp_ptr)
    : RecurringAction(start_gen, end_gen),
      landscape_(landscape),
      species_ptr_(sp_ptr),
      num_gens_counted_(0),
      sum_of_proportions_(landscape.size(), CellOrganismProvenanceProportions(landscape.size(), 0)) {
}

void MigrationTrackingRegime::log_provenances(const LandscapeOrganismProvenanceProportions& landscape_organism_provenances) {
    assert(landscape_organism_provenances.size() == this->sum_of_proportions_.size());
    for (CellIndexType i = 0; i < this->sum_of_proportions_.size(); ++i) {
        CellOrganismProvenanceProportions copp = landscape_organism_provenances[i];
        assert(copp.size() == this->sum_of_proportions_[i].size());
        for (CellIndexType j = 0; j < this->sum_of_proportions_.size(); ++j) {
            this->sum_of_proportions_[i][j] += landscape_organism_provenances[i][j];
        }
    }
    ++this->num_gens_counted_;
}

LandscapeOrganismProvenanceProportions MigrationTrackingRegime::calc_mean_proportions() const {
    assert(this->num_gens_counted_ > 0);
    LandscapeOrganismProvenanceProportions mean_proportions(this->landscape_.size(), CellOrganismProvenanceProportions(this->landscape_.size(), 0));
    for (CellIndexType i = 0; i < this->sum_of_proportions_.size(); ++i) {
        for (CellIndexType j = 0; j < this->sum_of_proportions_.size(); ++j) {
            mean_proportions[i][j] = this->sum_of_proportions_[i][j] / this->num_gens_counted_;
        }
    }
    return mean_proportions;
}

float MigrationTrackingRegime::total_influx_rate(CellIndexType cur_cell_index,
        const LandscapeOrganismProvenanceProportions& mean_proportions) const {
    float total_influx = 0.0;
    for (CellIndexType j = 0; j < mean_proportions[cur_cell_index].size(); ++j) {
        if (j != cur_cell_index) {
            total_influx += mean_proportions[cur_cell_index][j];
        }
    }
    return total_influx;
}

void MigrationTrackingRegime::write(std::ofstream& out, const std::string& separator) const {
    LandscapeOrganismProvenanceProportions mean_proportions = this->calc_mean_proportions();

    // header
    out << "cell.index";
    out << separator << "influx.rate";
    out << separator << "(x,y)";
    for (CellIndexType i = 0; i < mean_proportions.size(); ++i) {
        out << separator << "(" << this->landscape_.index_to_x(i) << "," << this->landscape_.index_to_y(i) << ")";
    }
    out << std::endl;

    // main rows
    for (CellIndexType i = 0; i < mean_proportions.size(); ++i) {
        out << i;
        out << separator << this->total_influx_rate(i, mean_proportions);
        out << separator << "(" << this->landscape_.index_to_x(i) << "," << this->landscape_.index_to_y(i) << ")";
        for (CellIndexType j = 0; j < mean_proportions.size(); ++j) {
            out << separator << mean_proportions[i][j];
        }
        out << std::endl;
    }
}


// MigrationTrackingRegime
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// World

// lifecycle ##################################################################

// singleton instance
World World::instance_;

// constructor
World::World()
    : species_registry_(SpeciesRegistry::get_instance()),
      rng_(RandomNumberGenerator::get_instance()),
      num_fitness_traits_(1),
      global_selection_strength_(DEFAULT_GLOBAL_SELECTION_STRENGTH),
      generations_to_run_(0),
      current_generation_(0),
      logger_init_(false),
      log_frequency_(10),
      is_log_to_screen_(true),
      allow_multifurcations_(true),
      is_produce_final_output_(true),
      is_produce_full_complement_diploid_trees_(false) {
    this->current_generation_ = 0;
}

World::~World() {}

// set up #####################################################################

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
//    std::vector<float> default_movement_probabilities(this->landscape_.size(), 1.0);
//    sp->set_movement_probabilities(default_movement_probabilities);
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
    this->logger_.info(pre_msg.str());

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
    this->logger_.info(post_msg.str());
}

// event handlers #############################################################

void World::add_environment_settings(GenerationCountType generation, const EnvironmentSettings& environment_settings) {
    this->environment_settings_[generation] = environment_settings;
}

void World::add_pre_reproduction_migration_tracker(const MigrationTrackingRegime& migration_tracking_regime) {
    this->pre_reproduction_migration_tracking_regimes_.push_back(migration_tracking_regime);
}

void World::add_post_dispersal_migration_tracker(const MigrationTrackingRegime& migration_tracking_regime) {
    this->post_dispersal_migration_tracking_regimes_.push_back(migration_tracking_regime);
}

void World::add_jump_dispersal_regime(const JumpDispersalRegime& jump_dispersal) {
    this->jump_dispersal_regimes_.push_back(jump_dispersal);
}

void World::add_tree_sampling(GenerationCountType generation, const SamplingRegime& sampling_regime) {
    this->tree_samples_.insert(std::make_pair(generation, sampling_regime));
}

void World::add_occurrence_sampling(GenerationCountType generation, Species * species_ptr) {
    this->occurrence_samples_.insert(std::make_pair(generation, species_ptr));
}

void World::set_initialization_regime(const InitializationRegime& initialization_regime) {
    this->initialization_regime_ = initialization_regime;
}

// simulation cycles ##########################################################

void World::run() {

    // ---- SETUP ----

    if (!this->logger_init_) {
        this->init_logger();
        this->logger_.hide_elapsed_time();
    }
    this->log_configuration();

#if defined(DEBUG)
    this->logger_.info("RUNNING DEBUG-MODE BUILD.");
#endif

#if defined(MEMCHECK)
    this->logger_.info("UNRELEASED NODE MEMORY WILL BE LOGGED.");
#endif

    // ---- INITIALIZATION ----

    this->logger_.reset_elapsed_time("I+");
    this->logger_.show_elapsed_time();
    this->logger_.info("Starting initialization cycles.");
    this->run_initialization_cycles();
    this->logger_.info("Initialization cycles ended.");

    // ---- MAIN RUN ----

    this->logger_.reset_elapsed_time("M+");
    this->logger_.info("Starting main run cycles.");
    this->run_main_cycles();
    this->logger_.info("Ending main run cycles.");

    // ---- POST-RUN ----

    if (this->is_produce_final_output_) {
        this->logger_.info("Saving final set of occurrences and trees for all species.");
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

}

void World::run_initialization_cycles() {

    // set initialization environment
        this->set_world_environment(this->initialization_regime_.environment, "[Initialization]");

    // create ancestral alleles for each locus for each species, recording
    // where the first locus for each species is
    std::vector<GenealogyNode *>  ancestral_genes;
    std::map<Species *, unsigned> species_ancestral_gene_index_start;
    for (SpeciesRegistry::iterator spi = this->species_registry_.begin();
            spi != this->species_registry_.end();
            ++spi) {
        species_ancestral_gene_index_start.insert(std::make_pair(*spi, ancestral_genes.size()));
        // haploid
        GenealogyNode * haploid = new GenealogyNode();
        haploid->increment_count();
        ancestral_genes.push_back(haploid);
        // diploid 1 and 2
        for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_loci; ++i) {
            GenealogyNode * diploid = new GenealogyNode();
            diploid->increment_count();
            ancestral_genes.push_back(diploid);
        }
    }

    // generate organisms, linking to ancestral alleles
    for (InitializationRegime::CellSpeciesPopulationMapType::iterator cpm = this->initialization_regime_.cell_populations.begin();
            cpm != this->initialization_regime_.cell_populations.end();
            ++cpm) {
        CellIndexType cell_index = cpm->first;
        Cell& cell = this->landscape_[cell_index];
        InitializationRegime::SpeciesPopulationSizeMapType& sp_pop_sizes = cpm->second;
        for (SpeciesRegistry::iterator spi = this->species_registry_.begin();
                spi != this->species_registry_.end();
                ++spi) {
            Species * sp = *spi;
            InitializationRegime::SpeciesPopulationSizeMapType::iterator sppi = sp_pop_sizes.find(sp);
            if (sppi != sp_pop_sizes.end()) {
                PopulationCountType pop_size = sppi->second;
                cell.generate_new_organisms(sp, pop_size);
                BreedingPopulation& pop = cell.population(sp);
                for (BreedingPopulation::iterator bpi = pop.begin();
                        bpi != pop.end();
                        ++bpi) {
                    Organism * o = *bpi;
                    unsigned int ancestor_gene_idx = species_ancestral_gene_index_start[sp];
                    assert( (ancestor_gene_idx + NUM_NEUTRAL_DIPLOID_loci - 1) < ancestral_genes.size() );
                    o->haploid_marker().link(ancestral_genes[ancestor_gene_idx]);
                    for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_loci; ++i) {
                        ++ancestor_gene_idx;
                        o->diploid_marker(i).link(ancestral_genes[ancestor_gene_idx],
                                ancestral_genes[ancestor_gene_idx]);
                    }
                }
            }
        }
    }

    // loop until all ancestral alleles have reference count of 3
    bool all_coalesced = false;
    unsigned long cycle_count = 0;
    while (!all_coalesced
            && (this->initialization_regime_.max_cycles < 0
                    || cycle_count < static_cast<GenerationCountType>(this->initialization_regime_.max_cycles) )) {
        cycle_count += 1;
        this->run_life_cycle();
        unsigned num_uncoalesced = 0;
        for (std::vector<GenealogyNode *>::iterator gni = ancestral_genes.begin();
                gni != ancestral_genes.end();
                ++gni) {
            if ( (*gni)->reference_count() > 3 ) {
                ++num_uncoalesced;
            }
        }
        if (num_uncoalesced == 0) {
            all_coalesced = true;
            std::ostringstream log_msg;
            log_msg << "Initialization cycle " << cycle_count;
            log_msg << ": all " << ancestral_genes.size() << " loci haved coalesced -- terminating initialization cycles.";
            this->logger_.info(log_msg.str());
        } else if ( this->initialization_regime_.max_cycles > 0
                && cycle_count >= static_cast<GenerationCountType>(this->initialization_regime_.max_cycles)) {
            std::ostringstream msg;
            msg << "Reached initialization cycle limit of " << this->initialization_regime_.max_cycles;
            msg << ": terminating initialization cycles";
            msg << " (" << num_uncoalesced << " of " << ancestral_genes.size() << " loci remaining to coalesce). ";
            this->logger_.info(msg.str());
        } else if ( cycle_count % this->log_frequency_ == 0 ) {
            std::ostringstream log_msg;
            log_msg << "Initialization cycle " << cycle_count;
            log_msg << ": " << num_uncoalesced << " of " << ancestral_genes.size() << " loci remaining to coalesce.";
            this->logger_.info(log_msg.str());
        }
    }

    // clean-up: decrement reference count self
    for (std::vector<GenealogyNode *>::iterator gni = ancestral_genes.begin();
            gni != ancestral_genes.end();
            ++gni) {
        (*gni)->decrement_count();
//        delete *gni;
    }

}

void World::run_main_cycles() {

    while (this->current_generation_ <= this->generations_to_run_) {

        // clear organism labels
        for (SpeciesRegistry::iterator spi = this->species_registry_.begin();
                spi != this->species_registry_.end();
                ++spi) {
            (*spi)->clear_organism_labels();
        }

        // process world changes
        this->process_environment_settings();

        // process dispersal events
        // this->process_dispersal_events();

        // run the life cycle
        if ( this->current_generation_ % this->log_frequency_ == 0) {
            std::ostringstream gen;
            gen << "Main cycle generation " << this->current_generation_ << " running.";
            this->logger_.info(gen.str());
        }
        this->run_life_cycle();

        // clear output filename stems
        this->output_filenames_.clear();

        // save occurrence data requested in this generation
        this->process_occurrence_samplings();

        // build trees requested in this generation
        this->process_tree_samplings();

        // increment generation count
        ++this->current_generation_;

    }

}

void World::run_life_cycle() {

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

    // get cell indexes in random order
    std::vector<CellIndexType> cell_indexes;
    cell_indexes.reserve(this->landscape_.size());
    for (CellIndexType i = 0; i != this->landscape_.size(); ++i) {
        cell_indexes.push_back(i);
    }
    RandomPointer cell_index_rp(this->rng_);
    std::random_shuffle(cell_indexes.begin(), cell_indexes.end(), cell_index_rp);

    // get jump dispersals in random order
    std::vector<JumpDispersalRegime *>  jump_dispersals_to_process;
    for (std::list<JumpDispersalRegime>::iterator jdi = this->jump_dispersal_regimes_.begin(); \
            jdi != this->jump_dispersal_regimes_.end();) {
        if (jdi->is_expired(this->current_generation_)) {
            jdi = this->jump_dispersal_regimes_.erase(jdi);
        } else {
            if (jdi->is_active(this->current_generation_)) {
                jump_dispersals_to_process.push_back(&(*jdi));
            }
            ++jdi;
        }
    }
    RandomPointer jump_dispersals_rp(this->rng_);
    std::random_shuffle(jump_dispersals_to_process.begin(), jump_dispersals_to_process.end(), jump_dispersals_rp);

    // pre-reproduction migration tracking
    this->process_migration_tracking(this->pre_reproduction_migration_tracking_regimes_, "dem-mig-");

    // reproduction
    for (std::vector<CellIndexType>::iterator i = cell_indexes.begin();
            i != cell_indexes.end();
            ++i) {
        this->landscape_[*i].reproduction();
    }

    // diffusion dispersal
    for (std::vector<CellIndexType>::iterator i = cell_indexes.begin();
            i != cell_indexes.end();
            ++i) {
        this->landscape_[*i].diffusion_dispersal();
    }

    // jump dispersal
    for (std::vector<JumpDispersalRegime *>::iterator jdi = jump_dispersals_to_process.begin();
            jdi != jump_dispersals_to_process.end();
            ++jdi) {
        JumpDispersalRegime& jd = *(*jdi);
        this->landscape_[jd.get_src_cell()].jump_dispersal(jd.get_species_ptr(), \
                jd.get_probability(), \
                jd.get_dest_cell());
    }

    // process migrants
    this->landscape_.process_migrants();

    // post-dispersal migration tracking
    this->process_migration_tracking(this->post_dispersal_migration_tracking_regimes_, "geo-mig-");

    // survival and competition
    for (std::vector<CellIndexType>::iterator i = cell_indexes.begin();
            i != cell_indexes.end();
            ++i) {
        this->landscape_[*i].survival();
        this->landscape_[*i].competition();
    }

}

void World::process_environment_settings() {
    std::map<GenerationCountType, EnvironmentSettings>::iterator wi = this->environment_settings_.find(this->current_generation_);
    if (wi == this->environment_settings_.end()) {
        return;
    }
    std::ostringstream log_leader;
    log_leader << "[Generation " << this->current_generation_  << "]";
    this->set_world_environment(wi->second, log_leader.str().c_str());
}

void World::set_world_environment(EnvironmentSettings& env, const char * log_leader) {
    if (env.carrying_capacity.size() != 0) {
        std::ostringstream msg;
        msg << log_leader << " Setting carrying capacity: \"" << env.carrying_capacity << "\".";
        this->logger_.info(msg.str());
        asciigrid::AsciiGrid<PopulationCountType> grid(env.carrying_capacity);
        this->landscape_.set_carrying_capacities(grid.get_cell_values());
    }
    if (env.fitness_trait_optima.size() != 0) {
        for (std::map<unsigned, std::string>::iterator ei = env.fitness_trait_optima.begin();
                 ei != env.fitness_trait_optima.end();
                 ++ei) {
            std::ostringstream msg;
            msg << log_leader << " Setting optima for fitness trait " <<  ei->first <<  ": \"" <<  ei->second <<  "\"";
            this->logger_.info(msg.str());
            asciigrid::AsciiGrid<FitnessTraitType> grid(ei->second);
            this->landscape_.set_environment(ei->first, grid.get_cell_values());
        }
    }
    if (env.movement_costs.size() != 0) {
        for (std::map<Species *, std::string>::iterator mi = env.movement_costs.begin();
                 mi != env.movement_costs.end();
                 ++mi) {
            std::ostringstream msg;
            msg << log_leader << " Setting movement costs for species " <<  mi->first->get_label() <<  ": \"" <<  mi->second <<  "\"";
            this->logger_.info(msg.str());
            asciigrid::AsciiGrid<MovementCountType> grid(mi->second);
            mi->first->set_movement_costs(grid.get_cell_values());
        }
    }
//    if (env.movement_probabilities.size() != 0) {
//        for (std::map<Species *, std::string>::iterator mi = env.movement_probabilities.begin();
//                 mi != env.movement_probabilities.end();
//                 ++mi) {
//            std::ostringstream msg;
//            msg << log_leader << " Setting movement probabilities for species " <<  mi->first->get_label() <<  ": \"" <<  mi->second <<  "\"";
//            this->logger_.info(msg.str());
//            asciigrid::AsciiGrid<float> grid(mi->second);
//            mi->first->set_movement_probabilities(grid.get_cell_values());
//        }
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

void World::process_migration_tracking(std::list<MigrationTrackingRegime>& migration_tracking_regimes, const char * title_prefix) {
    for (std::list<MigrationTrackingRegime>::iterator mti = migration_tracking_regimes.begin(); \
            mti != migration_tracking_regimes.end();) {
        if (mti->is_expired(this->current_generation_)) {
            std::ostringstream gen_span;
            gen_span << title_prefix;
            gen_span << std::setw(8) << std::setfill('0') << mti->get_start_gen();
            gen_span << "-";
            gen_span << std::setw(8) << std::setfill('0') << mti->get_end_gen();
            std::ofstream migration_tracking_output;
            this->open_ofstream(migration_tracking_output,
                this->compose_output_filename(mti->get_species_ptr()->get_label(), "", gen_span.str() + ".tsv"));
            mti->write(migration_tracking_output, "\t");
            mti = migration_tracking_regimes.erase(mti);
        } else {
            if (mti->is_active(this->current_generation_)) {
                LandscapeOrganismProvenanceProportions landscape_organism_provenances = this->landscape_.get_organism_provenances(mti->get_species_ptr());
                mti->log_provenances(landscape_organism_provenances);
            }
            ++mti;
        }
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

// results output #############################################################

void World::save_occurrences(Species * species_ptr ) {
    this->logger_.info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Saving occurrence data for species " + species_ptr->get_label() + ".");
    std::vector<PopulationCountType> counts;
    this->landscape_.count_organisms(species_ptr, counts);
    std::ofstream occs;
    this->open_ofstream(occs,
        this->compose_output_filename(species_ptr->get_label(), "", "occurrences.asc"));
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
    this->logger_.info(msg.str());

    std::vector<const Organism *> organisms;
    this->landscape_.sample_organisms(sp_ptr, num_organisms_per_cell, cell_indexes, organisms);
    std::string num_samples = convert::to_scalar<std::string>(organisms.size());

    if (organisms.size() == 0) {
        this->logger_.error("no organisms found in sample: aborting tree building");
        return;
    }

    this->logger_.info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Building single allele diploid loci trees for " + num_samples + " organisms (" + num_samples + " leaves per tree).");
    std::ofstream combined_trees;
    this->open_ofstream(combined_trees,
        this->compose_output_filename(sp_ptr->get_label(), label, "diploid1.tre"));
    this->write_diploid1_trees(sp_ptr, organisms, combined_trees);

    this->logger_.info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Building haploid locus tree for " + num_samples + " organisms (" + num_samples + " leaves per tree).");
    std::ofstream haploid_trees;
    this->open_ofstream(haploid_trees,
        this->compose_output_filename(sp_ptr->get_label(), label, "haploid.tre"));
    this->write_haploid_tree(sp_ptr, organisms, haploid_trees);

    if (this->is_produce_full_complement_diploid_trees_) {
        this->logger_.info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Building full complement diploid loci trees for " + num_samples + " organisms (" + convert::to_scalar<std::string>(organisms.size()*2) + " leaves per tree).");
        std::ofstream diploid_trees;
        this->open_ofstream(diploid_trees,
            this->compose_output_filename(sp_ptr->get_label(), label, "diploid2.tre"));
        this->write_diploid2_trees(sp_ptr, organisms, diploid_trees);
    }

    this->logger_.info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Building traits matrix for " + num_samples + " organisms (" + num_samples + " leaves per tree).");
    std::ofstream traits;
    this->open_ofstream(traits,
        this->compose_output_filename(sp_ptr->get_label(), label, "traits.nex"));
    this->write_traits(sp_ptr, organisms, traits);

}

// logging ####################################################################


void World::init_logger() {
    if (!this->logger_init_) {
        this->logger_.reset_elapsed_time();
        if (this->is_log_to_screen_) {
            this->logger_.add_handler(std::cout, LogHandler::LOG_INFO);
        }
        this->logger_.add_handler(filesys::compose_path(this->output_dir_, this->get_output_filename_stem() + ".run.log"), LogHandler::LOG_DEBUG);
        this->logger_.add_handler(filesys::compose_path(this->output_dir_, this->get_output_filename_stem() + ".err.log"), LogHandler::LOG_ERROR);
        this->logger_init_ = true;
    }
}

void World::log_tree_multiple_root_error(const std::string& species_label, unsigned long num_taxa) {
    std::ostringstream msg;
    msg << "tree for sample of organisms of species " << species_label;
    msg << " in generation " << this->current_generation_;
    msg << " could not be built due to multiple root nodes";
    msg << " (organism sample size = " << num_taxa << ").";
    this->logger_.error(msg.str());
}

// configuration logging ######################################################

void World::log_configuration() {
    std::ofstream out;
    this->open_ofstream(out, this->get_output_filename_stem() + ".conf.log");
    out <<  "GINKGO CONFIGURATION LOG " << this->logger_.get_current_timestamp() << std::endl;

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
    out << "*** GINKGO ***" << std::endl;
    out << "Title: " << this->title_ << std::endl;
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

    out << "*** SYSTEM ***" << std::endl;
    out << "Random seed: " << this->rng_.get_seed() << std::endl;
    out << "Number of fitness traits: " << this->num_fitness_traits_ << std::endl;
    out << "Global selection strength: " << this->global_selection_strength_ << std::endl;
    if (this->global_selection_strength_ == 0) {
        out << "WARNING: GLOBAL SELECTION STRENGTH IS 0: ENVIRONMENT HAS NO EFFECT!" << std::endl;
    }
    out << "Generations to run: " << this->generations_to_run_ << std::endl;
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
    out << "(" << this->environment_settings_.size() << " generations with environmental changes specified)" << std::endl;

    for (std::map<GenerationCountType, EnvironmentSettings>::iterator wi = this->environment_settings_.begin(); wi != this->environment_settings_.end(); ++wi) {
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
//        if (wi->second.movement_probabilities.size() != 0) {
//            for (std::map<Species *, std::string>::iterator mi = wi->second.movement_probabilities.begin();
//                     mi != wi->second.movement_probabilities.end();
//                     ++mi) {
//                out << "    Movement probablities for lineage \"" << mi->first->get_label() << "\": \"" << mi->second << "\"" << std::endl;
//            }
//        }
    }

   out << std::endl;
   out << "*** JUMP DISPERSALS ***" << std::endl;

    for (std::list<JumpDispersalRegime>::iterator jdi = this->jump_dispersal_regimes_.begin();
            jdi != this->jump_dispersal_regimes_.end();
            ++jdi)  {
        out << "Generations " << jdi->get_start_gen() << " to " << jdi->get_end_gen() << ": ";
        out << "dispersal of " << jdi->get_species_ptr()->get_label() << " ";
        out << " from (" << this->landscape_.index_to_x(jdi->get_src_cell()) << "," << this->landscape_.index_to_y(jdi->get_src_cell()) << ")";
        out << " to (" << this->landscape_.index_to_x(jdi->get_dest_cell()) << "," << this->landscape_.index_to_y(jdi->get_dest_cell()) << ")";
        out << " with probability " << jdi->get_probability() << "." << std::endl;
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

// file-handling helpers ######################################################

void World::open_ofstream(std::ofstream& out, const std::string& fpath) {
    std::string full_fpath = filesys::compose_path(this->output_dir_, fpath);
    out.open(full_fpath.c_str());
    if (!out) {
        throw WorldIOError("cannot open file \"" + full_fpath + "\" for output");
    }
}

std::string World::get_output_filename_stem() {
    if (this->output_filename_stem_.size() == 0) {
        if (this->title_.size() == 0) {
            this->output_filename_stem_ += "ginkgorun";
        } else {
            this->output_filename_stem_ += this->title_;
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

// World
///////////////////////////////////////////////////////////////////////////////

