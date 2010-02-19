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

#if !defined(GINKGO_ORGANISM_H)
#define GINKGO_ORGANISM_H

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

#if defined(MEMCHECK)
#include "memcheck.hpp"
#endif

namespace ginkgo {

class Organism;
class Species;


///////////////////////////////////////////////////////////////////////////////
// GenealogyNode
/**
 * A single node of a genealogical tree.
 *
 * Represents a single node of a genealogy of a neutral marker, with pointers
 * to its parent node. Uses reference counting to tracking references to self,
 * and deletes self when no other node points to self as parent.
 */
class GenealogyNode {

    public:

        /** Constructs a node with no antecedents. */
        GenealogyNode()
        : parent_(NULL),
          reference_count_(1),
          cell_index_(0)
          { }

        /**
         * Ensures all pointers from and to this object are nulled out before
         * winking out of existence.
         */
        ~GenealogyNode() {
            assert(this->parent_ != this);
            if (this->parent_) {
                this->parent_->decrement_count();
                this->parent_ = NULL;
            }
            assert(this->reference_count_ == 0 || this->reference_count_ == 1);
#if defined(MEMCHECK)
    if (UNRELEASED_NODES_LOG.is_open()) {
        UNRELEASED_NODES_LOG << "GenealogyNode " << this << " destroyed." << std::endl;
    }
#endif
        }

    private:

        /** Copy constructor (disabled by private scoping). */
        GenealogyNode(const GenealogyNode& );

        /** Assignment constructor (disabled by private scoping). */
        const GenealogyNode& operator=(const GenealogyNode&);

    public:

        /**
         * Inserts this node into a genealogy tree as child of given parent.
         *
         * Unlinks this node from its current position in the genealogy, and
         * then requests that the parent add this node as a child (which takes
         * care of all the new linking).
         *
         * @param parent    the immediate ancestor of this node in the tree
         */
        void link(GenealogyNode * parent) {
            assert(parent != this);
            if (this->parent_ != NULL)
                this->parent_->decrement_count();
            this->parent_ = parent;
        	if (this->parent_ != NULL) {
                this->parent_->increment_count();
            }
        }

        /**
         * Registers one less reference to this node.
         *
         * Decrements the count of references to this node (number of objects
         * that point to this node). If this goes to 1 (which means that the
         * only object that references this node object is this node object
         * itself), then this deletes itself.
         */
        void decrement_count() {
            if (this->reference_count_ == 1) {
                delete this;
            } else {
                this->reference_count_ -= 1;
            }
        }

        /**
         * Registers one more reference to this node.
         *
         * Notes that one more object points to or otherwise references this
         * node object.
         */
        void increment_count() {
            this->reference_count_ += 1;
        }

        /**
         * Returns pointer to parent node.
         *
         * @return      pointer to parent node
         */
        GenealogyNode * get_parent() {
            return this->parent_;
        }

        /**
         * Sets pointer to parent node.
         *
         * @param parent    pointer to parent node
         */
        void set_parent(GenealogyNode * parent) {
            this->parent_ = parent;
        }

        /**
         * Returns the cell index of this node.
         */
        CellIndexType get_cell_index() {
            return this->cell_index_;
        }

        /**
         * Sets the cell index.
         * @param cell_index    index of cell occupied by current organism.
         */
        void set_cell_index(CellIndexType cell_index) {
            this->cell_index_ = cell_index;
        }

    private:

        /**
         * Pointer to the parent of this node (<code>NULL</code> if no
         * children).
         */
        GenealogyNode *     parent_;

        /**
         * Number of objects that point to or reference this object (including
         * this object itself).
         */
        unsigned            reference_count_;

        /**
         * Position of containing organism: for reconstruction.
         */
        CellIndexType       cell_index_;

};
// GenealogyNode
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// HaploidMarker
/**
 * The allele of a haploid marker in a particular organism.
 *
 * Wraps a single GenealogyNode pointer, representing the terminal node
 * of a genealogy of a haploid locus.
 */
class HaploidMarker {

    public:

        /** Default constructor. */
        HaploidMarker() : allele_(NULL) {}

        /**
         * Destructor ensures that dropped reference to allele (a
         * GenealogyNode object) is registered.
         */
        ~HaploidMarker() {
            if (this->allele_) {
                this->allele_->decrement_count();
            }
        }

    private:
        /** Copy constructor (disabled by private scoping). */
        HaploidMarker(const HaploidMarker& h) {
            *this = h;
        }

    public:

        /**
         * Assignment operator.
         * Calls on assignment operators of member objects.
         * @param g     the HaploidMarker being copied
         * @return      constant reference to self
         */
        const HaploidMarker& operator=(const HaploidMarker& g) {
            if (this->allele_ != NULL) {
                this->allele_->decrement_count();
            }
            this->allele_ = g.allele_;
            if (this->allele_ != NULL) {
                this->allele_->increment_count();
            }
            return *this;
        }

        /**
         * Connects the allele (a GenealogyNode) at this locus into a genealogy
         * by setting it as the child of the allele wrapped by
         * <code>parent</code>.
         * @param   parent  the HaploidMarker with the allele being inherited
         */
        void inherit(const HaploidMarker& parent) {
            assert(this->allele_ == NULL);
            assert(this != &parent);
            this->allele_ = new GenealogyNode();
            this->allele_->link(parent.allele_);
        }

        /**
         * Returns pointer to GenealogyNode representing allele at this locus.
         * @return GenealogyNode representing allele at this locus
         */
        GenealogyNode* node() const {
            return this->allele_;
        }

        /**
         * Swap allele pointers for rapid assignment.
         */
        void swap(HaploidMarker& h) {
            GenealogyNode * t = this->allele_;
            this->allele_ = h.allele_;
            h.allele_ = t;
        }

        /**
         * Georeferences the current alleles.
         * @param cell_index    index of cell occupied by current organism.
         */
        void set_cell_index(CellIndexType cell_index) {
            if (this->allele_ != NULL) {
                this->allele_->set_cell_index(cell_index);
            }
        }

    private:

        /**
         * A node in the genealogy of this locus representing the
         * alleles at this locus in a particular organism.
         */
        GenealogyNode *      allele_;
};
// HaploidMarker
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Manages the genealogies of a diploid locus.
/**
 * The alleles of a diploid marker in a particualr organism.
 *
 * Wraps two GenealogyNode pointers, representing the terminal nodes
 * of the genealogies of a diploid locus.
 */
class DiploidMarker {

    public:

        /** Default constructor. Initializes genealogy nodes/alleles. */
        DiploidMarker()
            : allele1_(NULL),
              allele2_(NULL) { }

        /**
         * Destructor ensures that dropped reference to alleles
         * (GenealogyNode objects) are registered.
         */
        ~DiploidMarker() {
			if (this->allele1_ != NULL) {
                this->allele1_->decrement_count();
            }
            if (this->allele2_ != NULL) {
                this->allele2_->decrement_count();
            }
        }

    private:
        /** Copy constructor (disabled by private scoping) */
        DiploidMarker(const DiploidMarker&);

    public:

        /**
         * Assignment operator.
         * Calls on assignment operators of member objects.
         * @param g     the HaploidMarker being copied
         * @return      constant reference to self
         */
        const DiploidMarker& operator=(const DiploidMarker& g) {
            if (this->allele1_ != NULL) {
                this->allele1_->decrement_count();
            }
            this->allele1_ = g.allele1_;
            if (this->allele1_ != NULL) {
                this->allele1_->increment_count();
            }
            if (this->allele2_ != NULL) {
                this->allele2_->decrement_count();
            }
            this->allele2_ = g.allele2_;
            if (this->allele2_ != NULL) {
                this->allele2_->increment_count();
            }
            return *this;
        }

        /**
         * Connects the alleles (the GenealogyNode objects) at this locus into
         * existing genealogies represented by <code>female</code> and <code>
         * male</code> markers.
         * @param female a DiploidMarker from which one allele at this locus
         *               will be selected
         * @param male   a DiploidMarker from which another allele at this
         *               locus will be selected
         */
        void inherit(const DiploidMarker& female,
                     const DiploidMarker& male,
                     RandomNumberGenerator& rng) {
            assert(this->allele1_ == NULL);
            this->allele1_ = new GenealogyNode();
            this->allele1_->link(rng.select(female.allele1_, female.allele2_));
            assert(this->allele2_ == NULL);
            this->allele2_ = new GenealogyNode();
            this->allele2_->link(rng.select(male.allele1_, male.allele2_));
        }

        /**
         * Returns pointer to GenealogyNode representing allele 1 at this locus.
         * @return GenealogyNode representing allele 1 at this locus
         */
        GenealogyNode * node1() const {
            return this->allele1_;
        }

        /**
         * Returns pointer to GenealogyNode representing allele 2 at this locus.
         * @return GenealogyNode representing allele 2 at this locus
         */
        GenealogyNode * node2() const {
            return this->allele2_;
        }

        /**
         * Returns random allele of the diploid complement.
         * @param rng   RandomNumberGenerator object
         * @return      random allele from pair in locus
         */
        GenealogyNode * random_node(RandomNumberGenerator& rng) const {
            return rng.select(this->allele1_, this->allele2_);
        }

        /**
         * Swap diploid allele pointers for rapid assignment.
         */
        void swap(DiploidMarker& d) {
            GenealogyNode * t1 = this->allele1_;
            GenealogyNode * t2 = this->allele2_;
            this->allele1_ = d.allele1_;
            this->allele2_ = d.allele2_;
            d.allele1_ = t1;
            d.allele2_ = t2;
        }

        /**
         * Georeferences the current alleles.
         * @param cell_index    index of cell occupied by current organism.
         */
        void set_cell_index(CellIndexType cell_index) {
            if (this->allele1_ != NULL) {
                this->allele1_->set_cell_index(cell_index);
            }
            if (this->allele2_ != NULL) {
                this->allele2_->set_cell_index(cell_index);
            }
        }

    private:

        /**
         * A node in the genealogy of this locus representing one of the
         * alleles at this locus in a particular organism.
         */
        GenealogyNode *      allele1_;

        /**
         * A node in the genealogy of this locus representing the other
         * allele at this locus in a particular organism.
         */
        GenealogyNode *      allele2_;

};
// DiploidMarker
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Organism
/**
 * A single organism of a population of a particular species.
 * Responsible for tracking (non-neutral) genotype and neutral marker
 * histories. Very lightweight, with most functionality delegated to other
 * classes.
 */
class Organism {
    public:

        /** Organism status */
        enum Status {
            // bits:
            // 0: 0/1 = alive/expired
            // 1: 0/1 = female/male
            Female = 0,     // 00 = alive, female
            Expired = 1,    // 01 = expired, female
            Male = 2,       // 10 = alive, male
            ExpiredMale = 3 // 11 = expired, male
        };

        // --- lifecycle and assignment ---

        /**
         * Constructor instantiates a new organism with species, identity,
         * (fitness) genotype, and sex.
         *
         * @param species           pointer to Species of this Organism
         * @param fitness_traits    array of fitness factor values representing
         *                          the inheritable portion of fitness of this
         *                          individual organism
         * @param sex               gender of this organism
         */
        Organism(Species * species,
                 const FitnessTraits& fitness_traits,
                 Organism::Status status)
                : species_(species),
                  fitness_(-1),
                  state_flags_(status) {
            memcpy(this->fitness_trait_genotypes_, fitness_traits, MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));
        }

        /**
         * Constructor instantiates a new organism with species, identity,
         * and sex.
         * @param species_index     index of Species object in species pool
         *                          of the World object to which this organism
         *                          belongs
         * @param sex               gender of this organism
         */
        Organism(Species * species,
                 Organism::Status status)
                : species_(species),
                  fitness_(-1),
                  state_flags_(status) {
        }

        /**
         * Copy constructor, delegate work to assignment operator.
         */
        Organism(const Organism& ind) {
            *this = ind;
        }

        /** Destructor. */
        ~Organism() { }

        /**
         * Assignment.
         * @param   ind     the organism being copied
         * @return          const reference to this object
         */
        const Organism& operator=(const Organism& ind) {
            if (this == &ind) {
                return *this;
            }
            this->species_ = ind.species_;
            memcpy(this->fitness_trait_genotypes_, ind.fitness_trait_genotypes_, MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));
            this->state_flags_ = ind.state_flags_;
            this->fitness_ = ind.fitness_;
            this->neutral_haploid_marker_ = ind.neutral_haploid_marker_;
            for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_loci; ++i) {
                this->neutral_diploid_markers_[i] = ind.neutral_diploid_markers_[i];
            }
            return *this;
        }

        /**
         * Swap one Organism's "contents" with another.
         */
        void swap(Organism& o2) {
            // species
            Species * t_species = this->species_;
            this->species_ = o2.species_;
            o2.species_ = t_species;

            // fitness_trait_genotypes
            FitnessTraits t_fitness_trait_genotypes;
            memcpy(t_fitness_trait_genotypes, this->fitness_trait_genotypes_, MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));
            memcpy(this->fitness_trait_genotypes_, o2.fitness_trait_genotypes_, MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));
            memcpy(o2.fitness_trait_genotypes_, t_fitness_trait_genotypes, MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));

            // sex
            int t_state_flags = this->state_flags_;
            this->state_flags_ = o2.state_flags_;
            o2.state_flags_ = t_state_flags;

            // fitness
            float t_fitness = this->fitness_;
            this->fitness_ = o2.fitness_;
            o2.fitness_ = t_fitness;

            // neutral haploid markers
            this->neutral_haploid_marker_.swap(o2.neutral_haploid_marker_);

            // neutral diploid markers
            for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_loci; ++i) {
                this->neutral_diploid_markers_[i].swap(o2.neutral_diploid_markers_[i]);
            }
        }

        /**
         * Less-than operator to rank organisms according to fitness (which
         * must have been precalculated for both this and the other organism).
         * @param   other   the organism to which this one is being compared
         * @return          <code>true</code> if this organism has a lower
         *                  fitness than the other, <code>false</code>
         *                  otherwise
         */
        bool operator<(const Organism& other) const {
            assert(this->fitness_ >= 0);
            assert(other.fitness_ >= 0);
            return this->fitness_ < other.fitness_;
        }

        /** Returns a reference to the haploid marker of this organism. */
        const HaploidMarker& haploid_marker() const {
            return this->neutral_haploid_marker_;
        }

        /** Returns a reference to the diploid marker of this organism. */
        const DiploidMarker& diploid_marker(unsigned idx) const {
            assert(idx < NUM_NEUTRAL_DIPLOID_loci);
            return this->neutral_diploid_markers_[idx];
        }

        /** Returns a reference to the haploid marker of this organism. */
        GenealogyNode * get_haploid_node() const {
            return this->neutral_haploid_marker_.node();
        }

        /** Returns a reference to allele 1 of the diploid marker of this organism. */
        GenealogyNode * get_diploid_node1(unsigned idx) const {
            assert(idx < NUM_NEUTRAL_DIPLOID_loci);
            return this->neutral_diploid_markers_[idx].node1();
        }

        /** Returns a reference to allele 1 of the diploid marker of this organism. */
        GenealogyNode * get_diploid_node2(unsigned idx) const {
            assert(idx < NUM_NEUTRAL_DIPLOID_loci);
            return this->neutral_diploid_markers_[idx].node2();
        }

        /** Returns a reference to random allele of the diploid marker of this organism. */
        GenealogyNode * get_diploid_random_node(unsigned idx, RandomNumberGenerator& rng) const {
            assert(idx < NUM_NEUTRAL_DIPLOID_loci);
            return this->neutral_diploid_markers_[idx].random_node(rng);
        }

        // --- fitness & survival ---

        /**
         * Returns the i-th fitness factor value.
         * @return  a floating point value corresponding to the i-th fitness
         *          factor.
         */
        FitnessTraitType get_fitness_trait_genotype(unsigned idx) const {
            assert( idx < MAX_FITNESS_TRAITS );
            return this->fitness_trait_genotypes_[idx];
        }

        /** Returns reference to the array of inheritable fitness factors. */
        const FitnessTraits& get_fitness_trait_genotypes() const {
            return this->fitness_trait_genotypes_;
        }

        void set_fitness_trait_genotypes(const FitnessTraits& fitness_traits) {
            memcpy(this->fitness_trait_genotypes_, fitness_traits, MAX_FITNESS_TRAITS*sizeof(FitnessTraitType));
        }

        /**
         * Returns the fitness score of this organism.
         * @return  a floating-point value [0, 1] representing the
         *          pre-calculated fitness of this organism given the its
         *          fitness_trait and environment
         */
        float get_fitness() const {
            return this->fitness_;
        }

        /**
         * Sets the fitness score of this organism.
         * @param fitness a floating-point value [0, 1] representing the
         *                pre-calculated fitness of this organism given the its
         *                fitness_trait and environment
         */
        void set_fitness(float fitness) {
            this->fitness_ = fitness;
        }

        /**
         * Returns <code>true</code> if the organism needs to be removed from
         * the simulation by a higher-power.
         * @return  <code>true</code> if organism should be removed, </code>
         *          false</code> otherwise
         */
        bool is_expired() const {
            return this->state_flags_ & Organism::Expired;
        }

        /**
         * Specifies that the organism needs to be removed from
         * the simulation by a higher-power.
         */
        void set_expired() {
            this->state_flags_ |= Organism::Expired;
        }

        // --- meta info ---

        /**
         * Returns reference to the Species object representing
         * the species of this object in the species pool of the World.
         *
         * @return  reference to Species of this organism
         */
        Species& species() const;

        /**
         * Returns <code>true</code> if this organism is male.
         *
         * @return  <code>true</code> if this organism is male
         */
        bool is_male() const {
            return this->state_flags_ & Organism::Male;
        }

        /**
         * Returns <code>true</code> if this organism is female.
         *
         * @return  <code>true</code> if this organism is female
         */
        bool is_female() const {
            return !this->is_male();
        }

        // --- inheritance ---

        /**
         * Composes the genotypic fitness component of this organism.
         *
         * Composes the genotypic (inheritable) fitness factors of this
         * organism by copying elements with uniform random probability from
         * either the male or female parent, with random mutation.
         *
         * @param female                the female parent
         * @param male                  the other parent
         * @param num_fitness_traits   the number of fitness factors (which
         *                              should be less than or equal to
         *                              MAX_FITNESS_TRAITS)
         * @param rng                   source of random numbers
         */
        void inherit_fitness_trait_genotypes(const Organism& female,
                                               const Organism& male,
                                               unsigned num_fitness_traits,
                                               const std::vector<float>& fitness_trait_inheritance_sd,
                                               RandomNumberGenerator& rng) {
            assert(num_fitness_traits <= MAX_FITNESS_TRAITS);
            for (unsigned i = 0; i < num_fitness_traits; ++i) {
                FitnessTraitType ff_value = static_cast<FitnessTraitType>(female.fitness_trait_genotypes_[i] + male.fitness_trait_genotypes_[i])/2
                            + rng.normal(0, fitness_trait_inheritance_sd.at(i));
//                            + rng.normal(0, 0.00001);
                            // OK:      0.0000001
                            //          0.000001
                            //          0.00001
                            // EXTINCT: 0.0001
//                FitnessTraitType ff_value = rng.select(female.fitness_trait_genotypes_[i], male.fitness_trait_genotypes_[i]);
//                std::cout << "F=" << female.fitness_trait_genotypes_[i];
//                std::cout << ", M=" << male.fitness_trait_genotypes_[i];
//                std::cout << ", sd=" << fitness_trait_inheritance_sd.at(i);
//                std::cout << ", Offspring=" << ff_value << std::endl;
                this->fitness_trait_genotypes_[i] = ff_value;
            }
        }

        /**
         * Inserts the alleles/genealogical nodes of this organism's neutral
         * markers into the genealogy of its parents.
         *
         * @param female                the female parent
         * @param male                  the other parent
         * @param rng                   source of random numbers
         */
        void inherit_genealogies(const Organism& female,
                                 const Organism& male,
                                 RandomNumberGenerator& rng,
                                 CellIndexType cell_index) {
            this->neutral_haploid_marker_.inherit(female.neutral_haploid_marker_);
            this->neutral_haploid_marker_.set_cell_index(cell_index);
            for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_loci; ++i) {
                this->neutral_diploid_markers_[i].inherit(female.neutral_diploid_markers_[i], male.neutral_diploid_markers_[i], rng);
                this->neutral_diploid_markers_[i].set_cell_index(cell_index);
            }
        }

        void set_cell_index(CellIndexType cell_index) {
            this->neutral_haploid_marker_.set_cell_index(cell_index);
            for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_loci; ++i) {
                 this->neutral_diploid_markers_[i].set_cell_index(cell_index);
            }
        }

        Species * species_ptr() const {
            return this->species_;
        }

    private:

        /**
         * pointer to Species object in the World species pool of the
         * species of this organism
         */
        Species *               species_;

        /** values of inheritable component of fitness */
        FitnessTraits           fitness_trait_genotypes_;

        /** genealogy inherited through female alone */
        HaploidMarker           neutral_haploid_marker_;

        /** diploid genealogies */
        DiploidMarker           neutral_diploid_markers_[NUM_NEUTRAL_DIPLOID_loci];

        /** to cache pre-calculated fitness values */
        float                   fitness_;

        /** alive/expired, female/male */
        int                     state_flags_;
};
// Organism
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// compare_organism_fitness
/**
 * Serves as key for Organism sorting / set insertion etc.; if fitnesses are
 * equal, then returns random order.
 */
bool compare_organism_fitness(const Organism * o1, const Organism * o2);
typedef bool(*CompareOrganismFitnessFuncPtrType)(const Organism *, const Organism *);

// compare_organism_fitness
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// collection of Organism objects
typedef std::vector<Organism>   OrganismVector;
//
///////////////////////////////////////////////////////////////////////////////

} // ginkgo namespace


///////////////////////////////////////////////////////////////////////////////
// Specialization of std::swap when dealing with Organisms (for efficiency)

namespace std
{
    /**
     * Template specialization of std::swap() for Organism objects.
     */
    template<>
    void swap(ginkgo::Organism& o1, ginkgo::Organism& o2);
}

// Specialization of std::swap
///////////////////////////////////////////////////////////////////////////////

#endif
