

#if !defined(GINKGO_SPECIES_H)
#define GINKGO_SPECIES_H

#include <algorithm>
#include <cassert>
#include <vector>
#include <iterator>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <cstring>
#include "ginkgo_defs.hpp"
#include "randgen.hpp"
#include "organism.hpp"

#if defined(MEMCHECK)
#include "memcheck.hpp"
#endif

namespace ginkgo {

/** collection of Species Pointers */
typedef std::map<std::string, Species *> SpeciesByLabel;

///////////////////////////////////////////////////////////////////////////////
// Species
/**
 * Processes and properties that determine the ecologies, reproduction
 * characteristics, fitness, etc. of organisms.
 */
class Species {

    public:

        // --- lifecycle and assignment ---

        /**
         * Constructor.
         *
         * @param label                     identifier string for this species
         * @param num_fitness_traits        number of factors in fitness function
         * @param global_selection_strength strength of selection
         * @param rng                       random number generator
         */
        Species(const std::string& label,
                unsigned num_fitness_traits,
                float global_selection_strength,
                RandomNumberGenerator& rng);

        /**
         * Destructor.
         */
        ~Species() {}

        // --- access and mutation ---

        /**
         * Returns label of this species.
         *
         * @return label of this species.
         */
        std::string get_label() const {
            return this->label_;
        }

        /**
         * Returns number of fitness factors.
         *
         * Returns the number of elements of the statically allocated
         * genotypic fitness factor array of each organism that will actually
         * be used.
         *
         * @return number of "active" fitness factors
         */
        unsigned get_num_fitness_traits() const {
            return this->num_fitness_traits_;
        }

        /**
         * Sets the number of fitness factors.
         *
         * Sets the number of elements of the statically allocated
         * genotypic fitness factor array of each organism that will actually
         * be used.
         *
         * @param i number of "active" fitness factors
         */
        void set_num_fitness_traits(unsigned num_fitness_traits) {
            this->num_fitness_traits_ = num_fitness_traits;
        }

        /**
         * Returns the mean number of offspring per female.
         *
         * @return mean number of offspring per female
         */
        PopulationCountType get_mean_reproductive_rate() const {
            return this->mean_reproductive_rate_;
        }

        /**
         * Sets the mean number of offspring per female.
         *
         * @param mean number of offspring per female
         */
        void set_mean_reproductive_rate(PopulationCountType rate) {
            this->mean_reproductive_rate_ = rate;
        }

        /**
         * Returns the movement currency or potential for an organism.
         *
         * @return the movement potential of an organism
         */
        MovementCountType get_movement_capacity() const {
            if (this->is_fixed_movement_capacity_) {
                return this->movement_capacity_;
            } else {
                return static_cast<MovementCountType>(this->rng_.weighted_index_choice(this->movement_capacity_probabilities_));
            }
        }

        bool is_fixed_movement_capacity() {
            return this->is_fixed_movement_capacity_;
        }

        const std::vector<float>& movement_capacity_probabilities() {
            return this->movement_capacity_probabilities_;
        }

        /**
         * Sets the movement currency or potential for an organism as a fixed
         * value.
         *
         * @param moves the movement potential of an organism
         */
        void set_movement_capacity_fixed(MovementCountType moves) {
            this->is_fixed_movement_capacity_ = true;
            this->movement_capacity_ = moves;
        }

        /**
         * Sets the movement currency or potential for an organism by specifying
         * the probability of different capacities.
         * For example, given a vector of {0.1, 0.1, 0.2, 0.2, 0.4}, this means
         * that movement capacity of 0 with 0.1 probability, 1 with 0.2
         * probability, 2 with 0.2 probability, and 3 with 0.4 probability.
         *
         * @param probs vector of probablities
         */
        void set_movement_capacity_probabilities(const std::vector<float>& probs) {
            this->is_fixed_movement_capacity_ = false;
            this->movement_capacity_probabilities_ = probs;
        }

        /**
         * Sets the cell-specific probability of movement of an organism of
         * this lineage on the landscape.
         *
         * @param movement_probabilities probability of movement
         */
        void set_movement_probabilities(const std::vector<float>& movement_probabilities) {
            this->movement_probabilities_ = movement_probabilities;
        }

        /**
         * Gets the probability of movement of an organism of this species in
         * a particular cell on the landscape.
         *
         * @param movement_probabilities probability of movement
         */
        float movement_probability(CellIndexType i) {
            assert( i < this->movement_probabilities_.size() );
            return this->movement_probabilities_[i];
        }

        /**
         * Defines the costs for entering cells in the landscape.
         *
         * Each element in the vector <code>costs</code> corresponds to a cell
         * of the same index on the landscape. The value of each element is
         * the 'cost' that will be 'paid' out of an organism's movement
         * capacity to enter the cell.
         *
         * @param costs vector of movement costs
         */
        void set_movement_costs(const std::vector<MovementCountType>& costs) {
            this->movement_costs_ = costs;
        }

        /**
         * Returns the cost to be paid by an organism to enter a particular
         * cell.
         * @param  i    index of the cell to be entered
         * @return      value to be deducted from organism's movement capacity
         */
        MovementCountType movement_cost(CellIndexType i) {
//            std::cout << i << ", " << this->movement_costs_.size() << std::endl;
            assert( i < this->movement_costs_.size() );
            return this->movement_costs_[i];
        }

        /**
         * Sets the strengths of individual fitness factors on the fitness as
         * calculated for this species.
         *
         * Each element in the vector is a coefficient for the fitness function
         * calculated using the corresponding elements of the environmental and
         * genotypic fitness factor vectors.
         *
         * @param  strengths    vector of coefficients to the multivariate
         *                      fitness function
         */
        void set_selection_weights(const std::vector<float>& strengths) {
            float sum = 0;
            for (std::vector<float>::const_iterator i = strengths.begin();
                    i != strengths.end();
                    ++i) {
                sum += *i;
            }
            this->selection_weights_.clear();
            this->selection_weights_.reserve(strengths.size());
            for (std::vector<float>::const_iterator i = strengths.begin();
                    i != strengths.end();
                    ++i) {
                this->selection_weights_.push_back(*i / sum);
            }
        }

        /**
         * Gets the strengths of individual fitness factors on the fitness as
         * calculated for this species.
         *
         * Each element in the vector is a coefficient for the fitness function
         * calculated using the corresponding elements of the environmental and
         * genotypic fitness factor vectors.
         *
         * @param  strengths    vector of coefficients to the multivariate
         *                      fitness function
         */
        std::vector<float> get_selection_weights() {
            return this->selection_weights_;
        }

        /**
         * Sets the strengths of individual fitness factors on the fitness as
         * calculated for this species.
         *
         * Each element in the vector is a coefficient for the fitness function
         * calculated using the corresponding elements of the environmental and
         * genotypic fitness factor vectors.
         *
         * @param  strengths    vector of coefficients to the multivariate
         *                      fitness function
         */
        void set_default_fitness_trait_genotypes(const FitnessTraits& fitness_traits) {
            memcpy(this->default_fitness_trait_genotypes_, fitness_traits, MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));
        }

        /**
         * Sets the strengths of individual fitness factors on the fitness as
         * calculated for this species.
         *
         * Each element in the vector is a coefficient for the fitness function
         * calculated using the corresponding elements of the environmental and
         * genotypic fitness factor vectors.
         *
         * @param  strengths    vector of coefficients to the multivariate
         *                      fitness function
         */
        void set_default_fitness_trait_genotypes(std::vector<FitnessTraitType> fitness_traits) {
            assert(fitness_traits.size() <= MAX_FITNESS_TRAITS);
            for (unsigned i = 0; i != fitness_traits.size(); ++i) {
                this->default_fitness_trait_genotypes_[i] = fitness_traits[i];
            }
        }

        std::vector<FitnessTraitType> get_default_fitness_trait_genotypes() {
            std::vector<FitnessTraitType> v;
            v.reserve(this->get_num_fitness_traits());
            for (unsigned i = 0; i < this->get_num_fitness_traits(); ++i) {
                v.push_back(this->default_fitness_trait_genotypes_[i]);
            }
            return v;
        }


        /**
         * Sets the strengths of individual fitness factors on the fitness as
         * calculated for this species.
         *
         * Each element in the vector is a coefficient for the fitness function
         * calculated using the corresponding elements of the environmental and
         * genotypic fitness factor vectors.
         *
         * @param  sd    vector of standard deviations
         */
        void set_fitness_trait_inheritance_sd(const std::vector<float>& sd) {
            this->fitness_trait_inheritance_sd_ = sd;
        }

        /**
         * Gets the strengths of individual fitness factors on the fitness as
         * calculated for this species.
         *
         * Each element in the vector is a coefficient for the fitness function
         * calculated using the corresponding elements of the environmental and
         * genotypic fitness factor vectors.
         *
         * @param  variances    vector of coefficients to the multivariate
         *                      fitness function
         */
        std::vector<float> get_fitness_trait_inheritance_sd() {
            return this->fitness_trait_inheritance_sd_;
        }

        /**
         * Returns the fitness for a particular organism in a particular
         * environment.
         *
         * The fitness score is calculated based on Fisher's geometric fitness
         * function: multivariate gaussian distance between elements of the
         * organism's genotypic fitness factors and the environment's fitness
         * factors, with each distance weighted by the species-specific
         * selection strength for that element.
         *
         * \f[
         *      F_{j,k} = \sum_{i=1}^{i=Q}{S_i * (E_i - G_i)^2}
         * \f]
         *
         * Where \f$F_{j,k}\f$, is the fitness for organism \f$j\f$ in cell
         * \f$k\f$, \f$Q\f$ is the number of fitness factors, \f$S_i\f$ is the
         * selection strength for this species for factor \f$i\f$, \f$E_i\f$ is
         * the \f$i^{th}\f$ environmental factor in cell \f$k\f$ and \f$G_i\f$
         * is then \f$i^{th}\f$ genotypic factor for organism \f$j\f$.
         *
         * @param  organism     the organism
         * @param  environment  the vector of fitness factors of the
         *                      environment of the organism
         * @return              the fitness score for the organism in the given
         *                      environment
         */
        float calc_fitness(const Organism * organism_ptr, const FitnessTraits fitness_trait_optimum_optima) const {
            const FitnessTraits& fitness_traits = organism_ptr->get_fitness_trait_genotypes();
            const FitnessTraitType * g = fitness_traits;
            const FitnessTraitType * e = fitness_trait_optimum_optima;
            std::vector<float>::const_iterator s = this->selection_weights_.begin();
            float weighted_distance = 0.0;
            float diff = 0.0;
            for (unsigned i = 0; i < this->num_fitness_traits_; ++i, ++g, ++e, ++s) {
                diff = *e - *g;
                weighted_distance += (diff * diff) * *s; // each distance weighted by selection strength
            }
            weighted_distance = this->global_selection_strength_ * weighted_distance;
            float fitness = exp(-weighted_distance);
            return fitness;
        }

        // --- organism generation and reproduction ---

        /**
         * Resets the count of the number of labels produced.
         */
        void clear_organism_labels() {
            this->organism_labels_.clear();
            this->organism_label_index_ = 0;
        }

        /**
         * Returns a new unique (from the the last resetting of the label
         * index) organism label.
         *
         * @param   x       x-coordinate of organism's location
         * @param   y       y-coordiante of organism's location
         * @return          OTU label
         */
        std::string compose_organism_label(CellIndexType i, CellIndexType x, CellIndexType y) {
            std::ostringstream label_ostr;
            label_ostr << this->label_;
            label_ostr << "_i" << i;
            label_ostr << "_x" << x;
            label_ostr << "_y" << y;
            label_ostr << "_" << this->organism_label_index_++;
            return label_ostr.str();
        }

        /**
         * Returns a label for the given organism.
         *
         * @param   organism_ptr    pointer to organism to be labelled
         * @return                  unique label for organism
         */
         const std::string& get_organism_label(const Organism * organism_ptr) {
            std::map<const Organism *, std::string>::const_iterator ol = this->organism_labels_.find(organism_ptr);
            if (ol == this->organism_labels_.end()) {
                throw std::runtime_error("label requested, but not assigned to organism");
            } else {
                return ol->second;
            }
         }

        /**
         * Returns a label for the given organism, creating a new one if
         * it has not already been assigned.
         *
         * @param   organism_ptr    pointer to organism to be labelled
         * @param   x           x-coordinate of organism's location
         * @param   y           y-coordinate of organism's location
         * @return              unique label for organism
         */
         const std::string& set_organism_label(const Organism * organism_ptr, CellIndexType i, CellIndexType x, CellIndexType y) {
            this->organism_labels_[organism_ptr] = this->compose_organism_label(i, x, y);
            std::map<const Organism *, std::string>::const_iterator ol = this->organism_labels_.find(organism_ptr);
            return ol->second;
         }

        /**
         * Erases a label assigned to an organism.
         *
         * @param       organism_ptr label mapping to be erased
         */
         void erase_organism_label(const Organism * organism_ptr) {
            this->organism_labels_.erase(organism_ptr);
         }

        /**
         * Returns male/female with uniform probability.
         *
         * @return      Organism::Male or Organism::Female
         */
        Organism::Status get_random_sex(float female_threshold=0.5) const {
            if (this->rng_.uniform_01() < female_threshold) {
                return Organism::Male;
            } else {
                return Organism::Female;
            }
        }

        /**
         * Returns new organism, created de novo.
         *
         * @return  a default Organism object.
         */
        Organism * new_organism(CellIndexType cell_index) {
            Organism * organism_ptr = this->organism_memory_manager_.allocate();
            organism_ptr->init(this,
                    this->default_fitness_trait_genotypes_,
                    this->get_random_sex());
            organism_ptr->set_cell_index(cell_index);
            return organism_ptr;
        }

        /**
         * Returns an offspring produced by two organisms.
         *
         * Fitness traits are inherited as the mid-parent value for a trait plus
         * a random error. Genealogies are inherited by selecting an allele from
         * each parent for each locus.
         *
         * @param female    female parent
         * @param male      male parent
         * @return          offspring
         */
        Organism * new_organism(const Organism* female_ptr, const Organism* male_ptr, CellIndexType cell_index, bool evolve_fitness_components) {
            Organism * organism_ptr = this->organism_memory_manager_.allocate();
            organism_ptr->init(this, this->get_random_sex());
            if (evolve_fitness_components) {
                organism_ptr->inherit_fitness_trait_genotypes(
                    *female_ptr,
                    *male_ptr,
                     this->num_fitness_traits_,
                     this->fitness_trait_inheritance_sd_,
                     this->rng_
                );
            } else {
                organism_ptr->set_fitness_trait_genotypes(this->default_fitness_trait_genotypes_);
            }
            organism_ptr->inherit_genealogies(*female_ptr, *male_ptr, this->rng_, cell_index);
            return organism_ptr;
        }

    private:
        /** Copy constructor, private declaration with no definition to disabled */
        Species(const Species&);
        /** Assignment, private declaration with no definition to disabled */
        const Species& operator=(const Species&);

    private:
        /** unique identifier for this species */
        std::string                         label_;
        /** number of active fitness factors */
        unsigned                            num_fitness_traits_;
        /** weighting of the least-squares distance between organism trait value and environmental optimum */
        std::vector<float>                  selection_weights_;
        /** variance in the inheritance equation */
        std::vector<float>                  fitness_trait_inheritance_sd_;
        /** mean number of offspring per female */
        PopulationCountType                 mean_reproductive_rate_;
        /** landscape migration potential for this species */
        std::vector<MovementCountType>      movement_costs_;
        /** is movement capacity fixed or variable? */
        bool                                is_fixed_movement_capacity_;
        /** movement potential of each organism at the start of each round */
        MovementCountType                   movement_capacity_;
        /** variable movement-capacity assignment: vector of probabilities associated with each movement capacity */
        std::vector<float>                  movement_capacity_probabilities_;
        /** cell- and lineage-specific probability of movement of each organism */
        std::vector<float>                  movement_probabilities_;
        /** genotype for organisms created de novo */
        FitnessTraits                       default_fitness_trait_genotypes_;
        /** global selection strength */
        float                               global_selection_strength_;
        /** source of random numbers of various distributions */
        RandomNumberGenerator&              rng_;
        /** tracks the number of labels assigned to organisms of this species */
        unsigned long                       organism_label_index_;
        /** tracks the organism to label assignment */
        std::map<const Organism *, std::string>   organism_labels_;
        /** memory pool */
        OrganismMemoryManager&              organism_memory_manager_;

};
// Species
///////////////////////////////////////////////////////////////////////////////


} // namespace ginkgo

#endif
