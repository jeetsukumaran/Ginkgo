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

#if !defined(GINKGO_WORLD_H)
#define GINKGO_WORLD_H

#include <map>
#include <utility>
#include <istream>
#include <fstream>
#include <stdexcept>
#include <set>
#include <string>

#include "ginkgo_defs.hpp"
#include "randgen.hpp"
#include "cell.hpp"
#include "landscape.hpp"
#include "tree.hpp"
#include "logging.hpp"

namespace ginkgo {

///////////////////////////////////////////////////////////////////////////////
// Exceptions

/**
 * General i/o error.
 */
class WorldIOError : public std::runtime_error {
    public:
        WorldIOError(const char * msg) : std::runtime_error(msg) {}
        WorldIOError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

// Exceptions
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// RecurringAction

class RecurringAction {

    public:
        RecurringAction(GenerationCountType start_gen, GenerationCountType end_gen);

        bool is_active(GenerationCountType current_gen);
        bool is_completed(GenerationCountType current_gen);

        GenerationCountType get_start_gen() {
            return this->start_gen_;
        }
        void set_start_gen(GenerationCountType gen) {
            this->start_gen_ = gen;
        }

        GenerationCountType get_end_gen() {
            return this->end_gen_;
        }
        void set_end_gen(GenerationCountType gen) {
            this->end_gen_ = gen;
        }

    private:
        GenerationCountType     start_gen_;
        GenerationCountType     end_gen_;
};

// RecurringAction
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// JumpDispersalRegime

class JumpDispersalRegime : public RecurringAction {

    public:
        JumpDispersalRegime(GenerationCountType start_gen,
                GenerationCountType end_gen,
                Species * sp_ptr,
                float probability,
                CellIndexType src_cell,
                CellIndexType dest_cell);

        Species *  get_species_ptr() {
            return this->species_ptr_;
        }
        void set_species_ptr(Species *  sp_ptr) {
            this->species_ptr_ = sp_ptr;
        }

        float get_probability() {
            return this->probability_;
        }
        void set_probability(float prob) {
            this->probability_ = prob;
        }

        CellIndexType get_src_cell() {
            return this->src_cell_;
        }
        void set_src_cell(CellIndexType cell_idx) {
            this->src_cell_ = cell_idx;
        }

        CellIndexType get_dest_cell() {
            return this->dest_cell_;
        }
        void set_dest_cell(CellIndexType cell_idx) {
            this->dest_cell_ = cell_idx;
        }

    private:
        /** Pointer to species. */
        Species *                   species_ptr_;
        /** Probability of dispersal. */
        float                       probability_;
        /** Origin cell index. */
        CellIndexType               src_cell_;
        /** Destination cell index. */
        CellIndexType               dest_cell_;
};

// JumpDispersalRegime
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// supporting data

/**
 * Sampling regime, tracking the number of organisms and list of cells to
 * be sampled.
 */
struct SamplingRegime {

    public:
        SamplingRegime()
            : num_organisms_per_cell(0) { }

    public:
        /** Pointer to species. */
        Species *                   species_ptr;
        /**  Number of organisms from each cell to be sampled (0 = all). */
        PopulationCountType         num_organisms_per_cell;
        /** List of cell indexes to be sampled. */
        std::set<CellIndexType>     cell_indexes;
        /** Label prefix for file(s). */
        std::string                 label;
};

/**
 * Stochastic dispersal.
 */
struct DispersalEvent {

    public:
        DispersalEvent()
            : species_ptr(0),
              num_organisms(0),
              source(0),
              destination(0),
              probability(0) { }

    public:
        /** Pointer to species. */
        Species *                   species_ptr;
        /**  Number of organisms from cell to be sampled (0 = all). */
        PopulationCountType         num_organisms;
        /** Origin cell index. */
        CellIndexType               source;
        /** Destination cell index. */
        CellIndexType               destination;
        /** Probability of dispersal. */
        float                       probability;

};




/**
 * Changes to world geographical template at the start of a generation,
 * modelling climate change, changes in landscape etc.
 */
struct EnvironmentSettings {

    /** Path to grid setting carrying capacity. */
    std::string                             carrying_capacity;

    /**
     * Environmental regimes that need to be changed/set (expressed as factor
     * indexes mapped to ESRI ASCII Grid file paths).
     */
    std::map<unsigned, std::string>         fitness_trait_optima;

    /**
     * Movement costs that need to be changed/set. (expressed as species labels
     * mapped to ESRI ASCII Grid file paths).
     */
    std::map<Species *, std::string>        movement_costs;
//
//    /**
//     * Movement probabiities that need to be changed/set. (expressed as species labels
//     * mapped to ESRI ASCII Grid file paths).
//     */
//    std::map<Species *, std::string>        movement_probabilities;

}; // EnvironmentSettings


/**
 * Initialization of system.
 */
struct InitializationRegime {
    public:
        InitializationRegime() :  max_cycles(-1) { }
        typedef std::map<Species *, PopulationCountType>  SpeciesPopulationSizeMapType;
        typedef std::map<CellIndexType, SpeciesPopulationSizeMapType > CellSpeciesPopulationMapType;
        CellSpeciesPopulationMapType  cell_populations;
        EnvironmentSettings           environment;
        long                          max_cycles;
}; // initialize


// Supporting Data
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// World

/**
 * Meta-framework that binds everything together.
 */
class World {

    public:

        // lifecycle ##########################################################

        /**
         * Destructor, destroys Species and frees memory allocated to Species
         * objects.
         */
        ~World();

        // set up #############################################################

        /**
         * Build a landscape of the specified spatial and environmental
         * dimensions.
         *
         * @param size_x                the size of the Landscape from a
         *                              geospatial perspective in the
         *                              x-dimension
         * @param size_y                the size of the Landscape from a
         *                              geospatial perspective in the
         *                              y-dimension
         */
        void generate_landscape(CellIndexType size_x, CellIndexType size_y);

        /**
         * Instantiates a new species on the heap and inserts into species
         * pool.
         *
         * @param label     unique label identifying species
         */
        Species& new_species(const std::string& label);

        /**
         * Generates specified number of new Organism objects of the specified
         * Species, and inserts them into a Cell specified by its vector index.
         *
         * @param cell_index     the vector index of the Cell into which the
         *                       new Organism objects will be added
         * @param species_index  index of pointer to the Species object in the
         *                       Species pool of the Landscape/World
         * @param size           number of new Organism objects to generate
         * @param ancestral_pop_size       size of ancestral population of seed population (N; 0 => N = n )
         * @param ancestral_generations    number of generations in ancestral population  (0 => 10N)
         */
        void generate_seed_population(CellIndexType cell_index,
                Species * species_ptr,
                PopulationCountType pop_size,
                PopulationCountType ancestral_pop_size,
                GenerationCountType ancestral_generations);

        // event handlers #####################################################

        /**
         * Add a set of environments and settings that reconfigure the
         * world.
         *
         * @param   generation      generation number for this set of events
         *                          to be activated
         * @param   environment_settings  EnvironmentSettings data
         */
        void add_environment_settings(GenerationCountType generation, const EnvironmentSettings& environment_settings);

        /**
         * Add a jump dispersal regime.
         *
         * @param   jump_dispersal  a JumpDispersalRegime object describing
         *                          the regime.
         */
        void add_jump_dispersal_regime(const JumpDispersalRegime& jump_dispersal);

        /**
         * Add a directive to sample organisms and save a tree.
         *
         * @param   generation       generation number for this tree to be built
         * @param   sampling_regime  the leaves to add to the tree
         */
        void add_tree_sampling(GenerationCountType generation, const SamplingRegime& sampling_regime);

        /**
         * Add a directive to sample organisms and save their occurrence distribution.
         *
         * @param   generation       generation number for this occurrence be built
         * @param   species_ptr      pointer to species
         */
        void add_occurrence_sampling(GenerationCountType generation, Species * species_ptr);

        /**
         * Add a directive to seed a population.
         * @param cell_index     the vector index of the Cell into which the
         *                       new Organism objects will be added
         * @param species_index  index of pointer to the Species object in the
         *                       Species pool of the Landscape/World
         * @param size           number of new Organism objects to generate
         * @param ancestral_pop_size       size of ancestral population of seed population (N; 0 => N = n )
         * @param ancestral_generations    number of generations in ancestral population  (0 => 10N)
         */
//        void add_seed_population(CellIndexType cell_index,
//                Species * species_ptr,
//                PopulationCountType pop_size,
//                PopulationCountType ancestral_pop_size,
//                GenerationCountType ancestral_generations);

        // initialization #####################################################

        /**
         * Set up initialization regime.
         */
        void set_initialization_regime(const InitializationRegime& initialization_regime);

        // simulation cycles ##################################################

        /**
         * Run this simulation: initialization followded by main.
         */
        void run();

        /**
         * Run initialization cycles.
         */
        void run_initialization_cycles();

        /**
         * Run main simulation: set environment(s), single life-cycle, followed
         * by samplings.
         */
        void run_main_cycles();

        /**
         * A single cycle or generation of the simulation, including the
         * reproduction, migration, survival and competition phases.
         */
        void run_life_cycle();

        /**
         * Process world settings for the current generation.
         */
        void process_environment_settings();

        /**
         * Configure world according to given environment.
         */
        void set_world_environment(EnvironmentSettings& env, const char * log_leader);

        /**
         * Process tree building directives for the current generation.
         */
        void process_tree_samplings();

        /**
         * Process occurrence sampling directives for the current generation.
         */
        void process_occurrence_samplings();

#if defined(MEMCHECK)
        // --- memory check ---
        void run_final_cleanup_and_memory_check();
#endif

        // results output #####################################################

        /**
         * Save occurrence info for the current generation.
         * @param species_ptr   Pointer to species for which to save info.
         */
        void save_occurrences(Species * species_ptr);

        /**
         * Composes and writes a NEXUS file header.
         */
        void write_nexus_header(Species * sp_ptr,
                const std::vector<const Organism *>& organisms,
                std::ostream& out,
                bool add_diploid_allele_extensions=false);

        /**
         * Given a list of pointers to organisms, writes trait data.
         * @param   sp_ptr      pointer to species
         * @param   organisms   sample of organisms
         * @param   out         output stream
         */
        void write_traits(Species * sp_ptr,
                const std::vector<const Organism *>& organisms,
                std::ostream& out);

        /**
         * Given a list of pointers to organisms, builds and saves a tree
         * of the haploid locus allele to the given outputstream.
         * @param sp_ptr                    pointer to species
         * @param   organisms   sample of organisms
         * @param   out         output stream
         */
        void write_haploid_tree(Species * sp_ptr,
                const std::vector<const Organism *>& organisms,
                std::ostream& out);

        /**
         * Given a list of pointers to organisms, builds and saves trees
         * of the diploid locus alleles to the given outputstream.
         * @param   sp_ptr      pointer to species
         * @param   organisms   sample of organisms
         * @param   out         output stream
         */
        void write_diploid2_trees(Species * sp_ptr,
                const std::vector<const Organism *>& organisms,
                std::ostream& out);

        /**
         * Given a list of pointers to organisms, builds and saves trees
         * of the haploid and diploid locus alleles to the given outputstream,
         * with a single allele of each diploid complement chosen at random for
         * the diploid trees.
         * @param   sp_ptr      pointer to species
         * @param   organisms   sample of organisms
         * @param   out         output stream
         */
        void write_diploid1_trees(Species * sp_ptr,
                const std::vector<const Organism *>& organisms,
                std::ostream& out);

        /**
         * Samples organisms of specified species and specified cells of the
         * the landscape, and write out the corresponding trees.
         * Three tree files will be produced:
         *   (1) one for the haploid locus
         *   (2) one for all the diploid loci, with both alleles from each
         *       individual included
         *   (3) one for the haploid locus and all the diploid loci, with one
         *       allele sampled at random from each individual's diploid locus
         * @param sp_ptr                    pointer to species
         * @param num_organisms_per_cell    number of organisms of the given
         *                                  species per cell to sample
         *                                  (0 = sample all)
         * @param cell_indexes              list of cell indexes to sample
         * @param tree_filename             filename for trees
         */
        void save_trees(Species * sp_ptr,
                        PopulationCountType num_organisms_per_cell,
                        const std::set<CellIndexType>& cell_indexes,
                        const std::string& tree_filename_stem);

        // logging ############################################################

        /**
         * Configure the logger.
         */
        void init_logger();

        /**
         * Report coalescence failure/
         *
         * @param species_label     identifier of lineage/species that failed to coalesce
         * @param num_taxa          number of taxa
         */
        void log_tree_multiple_root_error(const std::string& species_label, unsigned long num_taxa);

        // configuration logging ##############################################

        /**
         * Dump out configuration file.
         */
        void log_configuration();


        // file-handling helpers ##############################################

        /**
         * Tries to open file, throwing exception if failed.
         * @param fpath     file path to open
         * @param ofstream  output file stream to use
         */
        void open_ofstream(std::ofstream& out, const std::string& fpath);

        /**
         * Returns label for output and reporting.
         * @return label for output and reporting
         */
        std::string get_output_filename_stem();


        /**
         * Composes and returns a name for a sampling (occurrence/tree) output
         * file.
         * @param   species_label   label for the species
         * @param   additional      any additional label info
         * @param   extension       file extension
         */
        std::string compose_output_filename(const std::string& species_label,
                const std::string& additional,
                const std::string& extension);

        // references to system-wide utilities ################################

        /**
         * Returns reference to this World's RandomNumberGenerator object.
         * @return reference to this World's RandomNumberGenerator object
         */
        RandomNumberGenerator& rng() {
            return this->rng_;
        }

        /**
         * Returns reference to this World's SpeciesRegistry object.
         * @return reference to this World's SpeciesRegistry object
         */
        SpeciesRegistry& species_registry() {
            return this->species_registry_;
        }

        /**
         * Returns reference to this World's Logger object.
         * @return reference to this World's Logger object
         */
        Logger& logger() {
            return this->logger_;
        }

        // accessors and mutators #############################################

        /**
         * Returns the current random number seed
         * return random number seed
         */
        unsigned long get_random_seed() const {
            return this->rng_.get_seed();
        }


        /**
         * Sets random number seed
         * @param seed random number seed
         */
        void set_random_seed(unsigned long seed) {
            this->rng_.set_seed(seed);
        }

        /**
         * Returns reference to this World's Landscape object.
         * @return reference to this World's Landscape object
         */
        Landscape& landscape() {
            return this->landscape_;
        }

        /**
         * Sets label for output and reporting.
         * @param label label for output and reporting
         */
        void set_title(std::string title) {
            this->title_ = title;
        }

        /**
         * Sets frequency of routine log messages in terms of numbers of
         * generations.
         * @param log_freq log output interval
         */
        void set_log_frequency(unsigned log_freq) {
            this->log_frequency_ = log_freq;
        }

        /**
         * Switch on/off final output generation.
         * @param val   <code>true</code> to produce final output.
         */
        void set_produce_final_output(bool val) {
            this->is_produce_final_output_ = val;
        }

        /**
         * Allow multifurcations in tree (otherwise 0-length branches will
         * be used to represent multifurcations.
         * @param val   <code>true</code> to produce trees with multifurcations.
         */
        void set_allow_multifurcations(bool val) {
            this->allow_multifurcations_ = val;
        }

        /**
         * Switch on/off full diploid tree generation.
         * @param val   <code>true</code> to produce trees tracking both alleles
         *              at each locus (by default, only haploid and subsampled
         *              diploid, i.e. one allele sampled at random from each
         *              diploid locus, produced).
         */
        void set_produce_full_complement_diploid_trees(bool val) {
            this->is_produce_full_complement_diploid_trees_ = val;
        }

        /**
         * Returns number of active fitness factors.
         * @return number of active fitness factors
         */
        unsigned get_num_fitness_traits() const {
            return this->num_fitness_traits_;
        }

        /**
         * Sets number of active fitness factors.
         * @param num_fitness_traits number of active fitness factors
         */
        void set_num_fitness_traits(unsigned num_fitness_traits) {
            this->num_fitness_traits_ = num_fitness_traits;
        }

        /**
         * Returns number of active fitness factors.
         * @return number of active fitness factors
         */
        float get_global_selection_strength() const {
            return this->global_selection_strength_;
        }

        /**
         * Sets number of active fitness factors.
         * @param global_selection_strength number of active fitness factors
         */
        void set_global_selection_strength(float global_selection_strength) {
            this->global_selection_strength_ = global_selection_strength;
        }

        /**
         * Sets the total number of generations to run.
         * @param ngens     number of generations to run
         */
        void set_generations_to_run(GenerationCountType ngens) {
            this->generations_to_run_ = ngens;
        }

        /**
         * Sets the output directory.
         * @param dir_path     existing directory path to which output files will
         *                     be written
         */
        void set_output_dir(const std::string& dir_path) {
            this->output_dir_ = dir_path;
        }

        /**
         * Sets the replicate id for a run (used to decorate output filenames).
         * @param rep_id     replicate indentification
         */
        void set_replicate_id(const std::string& rep_id) {
            this->replicate_id_ = rep_id;
        }

        /**
         * Returns number of cells in the landscape.
         * @param number of cells in the landscape
         */
        CellIndexType size() const {
            return static_cast<CellIndexType>(this->landscape_.size());
        }

        /**
         * Globally set individual cell carrying capacity.
         *
         * @param carrying_capacity the maximum number of organisms that can
         *                          occupy each cell at the end of every
         *                          generation
         */
        void set_global_cell_carrying_capacity(PopulationCountType carrying_capacity) {
            for (CellIndexType i = 0; i < this->landscape_.size(); ++i) {
                this->landscape_[i].set_carrying_capacity(carrying_capacity);
            }
        }

        ///////////////////////////////////////////////////////////////////////
        //// TO BE REFACTORED OUT
        ///////////////////////////////////////////////////////////////////////

        /**
         * Returns <code>true</code> if a species of the specified name/label
         * has been defined.
         * @return <code>true</code> if a species of the specified name/label
         * has been defined
         */
        bool has_species(const std::string& label) {
            return this->species_registry_.has_species(label);
        }


    ///////////////////////////////////////////////////////////////////////////
    // Singleton infrastructure

    public:
        static World& get_instance() {
            return World::instance_;
        }

    private:
        static World instance_;
        World();
        World(const World &);
        World & operator=(const World &);

    private:
        /** Name of this World (used for output files/reports). */
        std::string                             title_;
        /** Collection of pointers to the Species objects of this World. */
        SpeciesRegistry&                        species_registry_;
        /** The RandomNumberGenerator that is used by all objects of this World. */
        RandomNumberGenerator&                  rng_;
        /** The geospatial framework of this World. */
        Landscape                               landscape_;
        /** The number of dimensions to the fitness function. */
        unsigned                                num_fitness_traits_;
        /** The global strength of selection. */
        float                                   global_selection_strength_;
        /** Tracks the total number of generations to run. */
        GenerationCountType                     generations_to_run_;
        /** Tracks the number of generations that have been run. */
        GenerationCountType                     current_generation_;
        /** Output directory. */
        std::string                             output_dir_;
        /** Replicate id. */
        std::string                             replicate_id_;
        /** Output filename prefix. */
        std::string                             output_filename_stem_;
        /** Log control. */
        Logger                                  logger_;
        /** Flag that logger has been created and configured. */
        bool                                    logger_init_;
        /** Frequency (in # of generations) to log status. */
        unsigned                                log_frequency_;
        /** Duplicate log output to stdout/stderr? */
        bool                                    is_log_to_screen_;
        /** Allow multifurcations in tree (otherwise 0-length branches will be used to represent multifurcations. */
        bool                                    allow_multifurcations_;
        /** Produce final set of trees/occurrences, even if not requested? */
        bool                                    is_produce_final_output_;
        /** Produce full diploid trees (i.e., both alleles at each locus)? */
        bool                                    is_produce_full_complement_diploid_trees_;
        /** Initialization regime. */
        InitializationRegime                    initialization_regime_;
        /** Collection of events (key = generation #). */
        std::map<GenerationCountType, EnvironmentSettings>  environment_settings_;
        /** Collection of jump dispersal regimes */
        std::list<JumpDispersalRegime>          jump_dispersal_regimes_;
        /** Collection of tree building directives (key = generation #). */
        std::multimap<GenerationCountType, SamplingRegime> tree_samples_;
        /** Collection of occurrence description directives (key = generation #). */
        std::multimap<GenerationCountType, Species *> occurrence_samples_;
        /** Collection of seed population directives. */
        // std::vector<SeedPopulation>             seed_populations_;
        /** Track output filenames, so as to prevent clashes. */
        std::set<std::string>                   output_filenames_;

}; // World

// World
///////////////////////////////////////////////////////////////////////////////

} // ginkgo namespace

#endif
