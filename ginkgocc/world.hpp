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

namespace ginkgo {

/**
 * General i/o error.
 */
class WorldIOError : public std::runtime_error {
    public:
        WorldIOError(const char * msg) : std::runtime_error(msg) {}
        WorldIOError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

/**
 * Seed colonies.
 */
struct SeedPopulation {
    public:
        CellIndexType       cell_index;
        Species *           species_ptr;
        PopulationCountType pop_size;
        PopulationCountType ancestral_pop_size;
        GenerationCountType ancestral_generations;
    public:
        SeedPopulation()
            : cell_index(0),
              species_ptr(NULL),
              pop_size(0),
              ancestral_pop_size(0),
              ancestral_generations(0) { }

        SeedPopulation(CellIndexType cell_index,
            Species * species_ptr,
            PopulationCountType pop_size,
            PopulationCountType ancestral_pop_size,
            GenerationCountType ancestral_generations)
            : cell_index(cell_index),
              species_ptr(species_ptr),
              pop_size(pop_size),
              ancestral_pop_size(ancestral_pop_size),
              ancestral_generations(ancestral_generations) { }
};

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
        /** Label prefix for tree(s) / tree file(s). */
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
struct WorldSettings {

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

}; // WorldSettings

/**
 * Meta-framework that binds everything together.
 */
class World {

    public:

        // --- lifecycle --

        /**
         * Default constructor.
         */
        World();

        /**
         * Constructs a World with a given RNG seed.
         * @param seed  seed for the random number generator
         */
        World(unsigned long seed);

        /**
         * Destructor, destroys Species and frees memory allocated to Species
         * objects.
         */
        ~World();

        // --- access and mutation ---

        /**
         * Returns reference to this World's RandomNumberGenerator object.
         * @return reference to this World's RandomNumberGenerator object
         */
        RandomNumberGenerator& rng() {
            return this->rng_;
        }

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
         * Returns label for output and reporting.
         * @return label for output and reporting
         */
        std::string get_output_filename_stem();

        /**
         * Sets label for output and reporting.
         * @param label label for output and reporting
         */
        void set_label(std::string label) {
            this->label_ = label;
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
         * Returns number of cells in the landscape.
         * @param number of cells in the landscape
         */
        CellIndexType size() const {
            return static_cast<CellIndexType>(this->landscape_.size());
        }

        /**
         * Returns pointer to species if it exists.
         */
        Species * get_species_ptr(const std::string& label) {
            SpeciesByLabel::iterator sp = this->species_.find(label);
            assert(sp != this->species_.end());
            return (*sp).second;
        }

        /**
         * Returns <code>true</code> if a species of the specified name/label
         * has been defined.
         * @return <code>true</code> if a species of the specified name/label
         * has been defined
         */
        bool has_species(const std::string& label) {
            return (this->species_.find(label) != this->species_.end());
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

        /**
         * Set the costs for entering particular cells on the landscape for
         * a particular species.
         *
         * @param species_label  label of the Species object
         * @param costs         vector of costs for entering a cell, with
         *                      cost for cell \f$i\f$ in the landscape given
         *                      by element \f$i\f$ in the costs vector
         */
        void set_species_movement_costs(const std::string& species_label, const std::vector<MovementCountType>& costs) {
            assert(this->species_.find(species_label) != this->species_.end());
            assert(costs.size() == static_cast<CellIndexType>(this->landscape_.size()));
            this->species_[species_label]->set_movement_costs(costs);
        }

        /**
         * Set the costs for entering particular cells on the landscape for
         * a particular species.
         *
         * @param species_ptr   pointer to species object
         * @param costs         vector of costs for entering a cell, with
         *                      cost for cell \f$i\f$ in the landscape given
         *                      by element \f$i\f$ in the costs vector
         */
        void set_species_movement_costs(Species * species_ptr, const std::vector<MovementCountType>& costs) {
            assert(costs.size() == static_cast<CellIndexType>(this->landscape_.size()));
            species_ptr->set_movement_costs(costs);
        }

        /**
         * Sets the weight for each fitness factor for a particular species.
         *
         * @param species_label  label of the Species object
         * @param strengths     vector coeffecients to the Gaussian distance
         *                      equation used to evaluate fitness
         */
        void set_species_selection_weights(const std::string& species_label, const std::vector<float>& strengths) {
            assert(this->species_.find(species_label) != this->species_.end());
            assert(strengths.size() == this->num_fitness_traits_);
            this->species_[species_label]->set_selection_weights(strengths);
        }

        /**
         * Sets the weight for each fitness factor for a particular species.
         *
         * @param species_ptr   pointer to species object
         * @param strengths     vector coeffecients to the Gaussian distance
         *                      equation used to evaluate fitness
         */
        void set_species_selection_weights(Species * species_ptr, const std::vector<float>& strengths) {
            assert(strengths.size() == this->num_fitness_traits_);
            species_ptr->set_selection_weights(strengths);
        }


        /**
         * Sets the default genotypic (inheritable) component of fitness for a
         * organisms of the given species when generated de novo.
         *
         * @param species_label  label of the Species object
         * @param genotype      vector genotypic fitness values for a new
         *                      organism of the given species
         */
        void set_species_default_genotype(const std::string& species_label, const FitnessTraits& genotype) {
            assert(this->species_.find(species_label) != this->species_.end());
            this->species_[species_label]->set_default_fitness_trait_genotypes(genotype);
        }

        /**
         * Sets the default genotypic (inheritable) component of fitness for a
         * organisms of the given species when generated de novo.
         *
         * @param species_ptr   pointer to species object
         * @param genotype      vector genotypic fitness values for a new
         *                      organism of the given species
         */
        void set_species_default_genotype(Species * species_ptr, const FitnessTraits& genotype) {
            species_ptr->set_default_fitness_trait_genotypes(genotype);
        }

        // --- setup, initialization and seeding ---

        /**
         * Instantiates a new species on the heap and inserts into species
         * pool.
         *
         * @param label     unique label identifying species
         */
        Species& new_species(const std::string& label);

        /**
         * Generates specified number of new Organism objects of the specified
         * Species, and inserts them into a Cell specified by its geospatial
         * coordinates.
         *
         * @param x              the geospatial x-coordinate of the Cell into
         *                       which the new Organism objects will added
         * @param y              the geospatial y-coordinate of the Cell into
         *                       which the new Organism objects will added
         * @param species_label  label of the Species object
         * @param size           number of new Organism objects to generate
         */
//         void seed_population(CellIndexType x, CellIndexType y, const std::string& species_label, unsigned long size);

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

        // --- event handlers ---

        /**
         * Add a set of environments and settings that reconfigure the
         * world.
         *
         * @param   generation      generation number for this set of events
         *                          to be activated
         * @param   world_settings  WorldSettings data
         */
        void add_world_settings(GenerationCountType generation, const WorldSettings& world_settings);

        /**
         * Add a dispersal event.
         *
         * @param   generation      generation number for this set of events
         *                          to be activated
         * @param   dispersal_event descripion of event
         */
        void add_dispersal_event(GenerationCountType generation, const DispersalEvent& dispersal_event);

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
        void add_seed_population(CellIndexType cell_index,
                Species * species_ptr,
                PopulationCountType pop_size,
                PopulationCountType ancestral_pop_size,
                GenerationCountType ancestral_generations);

        // --- simulation cycles ---

        /**
         * A single cycle or generation of the simulation, including the
         * reproduction, migration, survival and competition phases.
         */
        void cycle();

        /**
         * Run multiple cycles or generations of the simulation.
         */
        void run();

        /**
         * Process world settings for the current generation.
         */
        void process_world_settings();

        /**
         * Process dispersal events for the current generation.
         */
        void process_dispersal_events();

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

        // --- logging and output ---

        /**
         * Dump out configuration file.
         */
        void log_configuration();

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
         * Given a list of pointers to organisms, writes trait data.
         * @param   sp_ptr      pointer to species
         * @param   organisms   sample of organisms
         * @param   out         output stream
         */
        void write_traits(Species * sp_ptr,
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

        void log_tree_multiple_root_error(const std::string& species_label, unsigned long num_taxa);

        /**
         * Tries to open file, throwing exception if failed.
         * @param fpath     file path to open
         * @param ofstream  output file stream to use
         */
        void open_ofstream(std::ofstream& out, const std::string& fpath);

        /**
         * Opens standard and error log file streams.
         */
        void open_logs();

        /**
         * Closes standard and error log file streams.
         */
        void close_logs();

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

        /**
         * Time stamp for log.
         * @return  formatted time stamp
         */
        std::string get_timestamp();

        /**
         * Write message to the general log stream (file only).
         * @param   message message to write
         */
        void log_detail(const std::string& message);

        /**
         * Write message to the general log stream (file and stdout).
         * @param   message message to write
         */
        void log_info(const std::string& message);

        /**
         * Write message to the error log stream (with duplicate to general
         * log stream).
         * @param   message message to write
         */
        void log_error(const std::string& message);

    private:
        /** Name of this World (used for output files/reports). */
        std::string                             label_;
        /** Collection of pointers to the Species objects of this World. */
        SpeciesByLabel                          species_;
        /** The RandomNumberGenerator that is used by all objects of this World. */
        RandomNumberGenerator                   rng_;
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
        /** Info log file stream. */
        std::ofstream                           infos_;
        /** Error log stream. */
        std::ofstream                           errs_;
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
        /** Collection of events (key = generation #). */
        std::map<GenerationCountType, WorldSettings>  world_settings_;
        /** Collection of dispersal events (key = generation #). */
        std::multimap<GenerationCountType, DispersalEvent>  dispersal_events_;
        /** Collection of tree building directives (key = generation #). */
        std::multimap<GenerationCountType, SamplingRegime> tree_samples_;
        /** Collection of occurrence description directives (key = generation #). */
        std::multimap<GenerationCountType, Species *> occurrence_samples_;
        /** Collection of seed population directives. */
        std::vector<SeedPopulation>             seed_populations_;
        /** Track output filenames, so as to prevent clashes. */
        std::set<std::string>                   output_filenames_;

    private:
        /** Disabled copy constructor. */
		World(const World &);
		/** Disabled assignment operator. */
		World & operator=(const World &);

}; // World


/**
 * Parses a configuration file, creating and returning a correspondingly
 * configured and populated World object.
 */
class WorldFactory {

    public:

        /**
         * Parses the given input stream and returns a World object based on
         * this.
         *
         * @param   src     reference to input stream describing the World
         * @return          a new World object
         */
        World build(std::istream& src);

}; // ConfigurationFileTokenizer


} // ginkgo namespace

#endif
