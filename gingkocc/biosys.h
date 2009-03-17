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

#if !defined(GINGKO_BIOSYS_H)
#define GINGKO_BIOSYS_H

/** @file biosys.h */

#include <cassert>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "gingko_defs.h"
#include "randgen.h"

namespace gingko {

class Organism;
class Species;

/** collection of pointers to Species objects */
typedef std::vector<Species *>  SpeciesPointerVector;

/** collection of Organism objects */
typedef std::vector<Organism>   OrganismVector;

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
          first_child_(NULL),
          next_sib_(NULL),
          reference_count_(1) { }
                  
        /** 
         * Ensures all pointers from and to this object are nulled out before 
         * winking out of existence. 
         */        
        ~GenealogyNode() {
            assert(this->parent_ != this);
            if (this->parent_) {                
                this->unlink();                         
            }                
            assert(this->first_child_ == NULL);                
            assert(this->reference_count_ == 0 || this->reference_count_ == 1);
        }
          
    private:		
		
        /** Copy constructor (disabled by private scoping). */   		
        GenealogyNode(const GenealogyNode& );

        /** Assignment constructor (disabled by private scoping). */ 
        const GenealogyNode& operator=(const GenealogyNode&);          
        
    public:        
          
        /**          
         * Nullifies all pointers to the specified node in this node and 
         * its children.
         *
         * Ensures that nothing points to <code>c</code>, either this node 
         * itself (through <code>first_child_</code>, or any of this node's 
         * children (through this node's children's <code>next_sib_</code>).
         * Also sets the child node's parent to NULL, and then decrements the 
         * count of references to self.
         *
         * @param child     pointer to the child node to be removed
         */
		void remove_child(GenealogyNode * child) {
			assert(this);
			GenealogyNode * g = this->first_child_;
			if (g == child) {
				this->first_child_ = g->next_sib_;
			} else { // parent's first child is not self ...
				// search for self amongst parent's children
				while (g->next_sib_ != child) {
					g = g->next_sib_;
					assert(g);
				}
				assert(g->next_sib_ == child);
				// set previous sib to point to self's next sib
				g->next_sib_ = child->next_sib_;
			}
			child->parent_ = NULL;
			this->decrement_count();
		}
		
        /**          
         * Adds a node as a child of this node.
         *
         * If this node has no children, then this node's 
         * <code>first_child_</code> pointer is set to point to the new child
         * node, otherwise the first null <code>next_sib_</code> is set to 
         * point to child. Also sets the child node's parent to this node, and 
         * increments count of references to self.
         *
         * @param child     pointer to the child node to be added
         */
		void add_child(GenealogyNode * child) {
			assert(this != child);
			child->parent_ = this;
			this->increment_count();
			if (this->first_child_ == NULL) {
				// parent has no children: self is first child.
				this->first_child_ = child;
			} else {
				// parent has children
				GenealogyNode * g = this->first_child_;
				// search for the first child with no sibling
				while (g->next_sib_ != NULL) {
					g = g->next_sib_;       
				}
				// insert self as sibling                    
				g->next_sib_ = child;
			}
		}		
		
        /**          
         * Removes this node from a genealogy tree.
         *
         * Disconnects this node from the pointer graph of a genealogy 
         * tree by asking this node's parent to remove self, and then nulling 
         * out this node's parent pointers.
         */
        void unlink() {
            if (this->parent_ == NULL)
            	return;
            this->parent_->remove_child(this);
			this->parent_ = NULL;
			this->next_sib_ = NULL;
        }
          
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
        	this->unlink();
        	if (parent == NULL) {
        	    return;
        	}
			parent->add_child(this);
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
            	assert(this->first_child_ == 0L);
                delete this;
            }
            this->reference_count_ -= 1;
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
         * Returns pointer to first child.
         * 
         * @return      pointer to first child
         */          
        GenealogyNode * get_first_child() {
            return this->first_child_;
        }
        
        /**
         * Sets pointer to first child.
         * 
         * @param first_child    pointer to first child
         */         
        void set_first_child(GenealogyNode * first_child) {
            this->first_child_ = first_child;
        }
        
        /**
         * Returns pointer to this node's (right) sibling.
         * 
         * @return      pointer to this node's right sibling
         */          
        GenealogyNode * get_next_sib() {
            return this->next_sib_;
        }
        
        /**
         * Sets pointer to this node's (right) sibling.
         * 
         * @param      pointer to this node's right sibling
         */         
        void set_next_sib(GenealogyNode * next_sib) {
            this->next_sib_ = next_sib;
        }
        
        /**
         * Returns the label of this node.
         *
         * The label serves to identify this node as an OTU on a tree 
         * representing the genealogy of this node.
         * 
         * @return      the label of the OTU represented by this node
         */          
        std::string get_label() const {
            return this->label_;
        }
                        
        /**
         * Sets the label of this node.
         *
         * The label serves to identify this node as an OTU on a tree 
         * representing the genealogy of this node.
         * 
         * @label      the label of the OTU represented by this node
         */          
        void set_label(const std::string& label) {
            this->label_ = label;                       
        }      
                
    private:
    
        /** 
         * Pointer to the parent of this node (<code>NULL</code> if no 
         * children).
         */
        GenealogyNode *     parent_;
        
        /** 
         * Pointer to the first child of this node (<code>NULL</code> if no 
         * children).
         */
        GenealogyNode *     first_child_;
        
        /** 
         * Pointer to next sibling of this node (<code>NULL</code> if this node
         * is the last child).
         */        
        GenealogyNode *     next_sib_;
        
        /** 
         * Number of objects that point to or reference this object (including
         * this object itself).
         */        
        unsigned            reference_count_;
        
        /** 
         * Identifier or label that will represent this node as an OTU on 
         * a tree.
         */            
		std::string         label_;
               
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
         * Sets the string representation or identification of the allele at 
         * this locus.
         * @param label     string representation of allele as an OTU on a tree
         */
        void set_label(const std::string& label) {       
            if (this->allele_) {
                this->allele_->set_label(label);
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
         * Sets the OTU label prefix for the nodes of the genealogies of this
         * locus.
         * @param label     the string for the OTU label prefix (generally, the
         *                  label of the organism to which this locus belongs)
         */
        void set_label(const std::string& label) {
            if (this->allele1_ != NULL) {
                this->allele1_->set_label(label + "_a");
            }              
            if (this->allele2_ != NULL) {
                this->allele2_->set_label(label + "_b");
            }             
        }
        
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
    
        /** Flags indicating gender for sexual reproduction */
        enum Sex {
            Male,
            Female
        };

        // --- lifecycle and assignment ---
        
        /** 
         * Constructor instantiates a new organism with species, identity,
         * (fitness) genotype, and sex.
         *
         * @param species_index     index of Species object in species pool
         *                          of the World object to which this organism 
         *                          belongs
         * @param label             string label that will identify alleles
         *                          belonging to this organism in genealogies
         * @param genotype          array of fitness factor values representing
         *                          the inheritable portion of fitness of this
         *                          individual organism
         * @param sex               gender of this organism
         */ 
        Organism(unsigned species_index,
                 const std::string& label,
                 const FitnessFactors& new_genotype, 
                 Organism::Sex new_sex) 
                : species_index_(species_index),
                  label_(label),
                  sex_(new_sex),
                  fitness_(-1),
                  expired_(false) {
            memcpy(this->genotypic_fitness_factors_, new_genotype, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));
        }
        
        /** 
         * Constructor instantiates a new organism with species, identity,
         * and sex.
         * @param species_index     index of Species object in species pool
         *                          of the World object to which this organism 
         *                          belongs
         * @param label             string label that will identify alleles
         *                          belonging to this organism in genealogies
         * @param sex               gender of this organism
         */        
        Organism(unsigned species_index, 
                 const std::string& label,        
                 Organism::Sex new_sex) 
                : species_index_(species_index), 
                  label_(label),                
                  sex_(new_sex),
                  fitness_(-1),
                  expired_(false) { }
        
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
            this->species_index_ = ind.species_index_;
            memcpy(this->genotypic_fitness_factors_, ind.genotypic_fitness_factors_, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));                    
            this->sex_ = ind.sex_;
            this->fitness_ = ind.fitness_;
            this->expired_ = ind.expired_;
            this->label_ = ind.label_;
            this->neutral_haploid_marker_ = ind.neutral_haploid_marker_;
            for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_LOCII; ++i) {
                this->neutral_diploid_markers_[i] = ind.neutral_diploid_markers_[i];
            }             
            return *this;
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
                   
        /** Returns reference to the array of inheritable fitness factors. */
        const FitnessFactors& genotype() const {
            return this->genotypic_fitness_factors_;
        }
        
        /** Returns a reference to the haploid marker of this organism. */
        const HaploidMarker& haploid_marker() const {
            return this->neutral_haploid_marker_;
        }
        
        /** Returns a reference to the diploid marker of this organism. */
        const DiploidMarker& diploid_marker(unsigned idx) const {
            assert(idx < NUM_NEUTRAL_DIPLOID_LOCII);
            return this->neutral_diploid_markers_[idx];
        }                
        
        // --- fitness & survival ---
        
        /** 
         * Returns the fitness score of this organism. 
         * @return  a floating-point value [0, 1] representing the 
         *          pre-calculated fitness of this organism given the its 
         *          genotype and environment
         */  
        float get_fitness() const {
            return this->fitness_;
        }
        
        /** 
         * Sets the fitness score of this organism. 
         * @param fitness a floating-point value [0, 1] representing the 
         *                pre-calculated fitness of this organism given the its 
         *                genotype and environment
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
            return this->expired_;
        }
        
        /** 
         * Specifies whether or not the organism needs to be removed from
         * the simulation by a higher-power. 
         * @param val  <code>true</code> if organism should be removed, </code>
         *             false</code> otherwise
         */           
        void set_expired(bool val) {
            this->expired_ = val;
        }          
        
        /** 
         * Returns the index of the pointer to the Species object representing
         * the species of this object in the species pool of the World.
         *
         * @return  index of pointer to Species object in species pool
         */   
        unsigned species_index() const {
            return this->species_index_;
        }
        
        /** 
         * Returns <code>true</code> if this organism is male.         
         *
         * @return  <code>true</code> if this organism is male
         */   
        bool is_male() const {
            return this->sex_ == Organism::Male;
        }
        
        /** 
         * Returns <code>true</code> if this organism is female.         
         *
         * @return  <code>true</code> if this organism is female
         */   
        bool is_female() const {
            return this->sex_ == Organism::Female;
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
         * @param num_fitness_factors   the number of fitness factors (which 
         *                              should be less than or equal to 
         *                              MAX_FITNESS_FACTORS)
         * @param mutation_rate         the probability of mutation per 
         *                              factor inheritance
         * @param max_mutation_size     the window (+/-) that a mutation event
         *                              can perturb a value
         * @param rng                   source of random numbers
         */
        void inherit_genotypic_fitness_factors(const Organism& female, 
                                               const Organism& male, 
                                               unsigned num_fitness_factors,
                                               float mutation_rate,
                                               int max_mutation_size,
                                               RandomNumberGenerator& rng) {
            assert(num_fitness_factors <= MAX_FITNESS_FACTORS);                                               
            for (unsigned i = 0; i < num_fitness_factors; ++i) {
                FitnessFactorType ff_value = rng.select(female.genotypic_fitness_factors_[i], male.genotypic_fitness_factors_[i]);  
                if (rng.uniform_real() < mutation_rate) {
                    ff_value += rng.uniform_int(-max_mutation_size, max_mutation_size);
                }
                this->genotypic_fitness_factors_[i] = ff_value;
            }
        }        

        void inherit_genealogies(const Organism& female, 
                                 const Organism& male,
                                 RandomNumberGenerator& rng) {
            // // std::cout << "INHERITING: " << this << " = " << &female << " + " << &male << std::endl;
            this->neutral_haploid_marker_.inherit(female.neutral_haploid_marker_);
            this->neutral_haploid_marker_.set_label(this->label_);        
            for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_LOCII; ++i) {
                this->neutral_diploid_markers_[i].inherit(female.neutral_diploid_markers_[i], male.neutral_diploid_markers_[i], rng);
                this->neutral_diploid_markers_[i].set_label(this->label_);
            }                
        }
        
    private:
        unsigned            species_index_;                                 // species 
        std::string         label_;                                         // label        
        FitnessFactors      genotypic_fitness_factors_;                     // non-neutral genotype: maps to fitness phenotype
        HaploidMarker        neutral_haploid_marker_;                        // track the genealogy of neutral genes in this organism    
        DiploidMarker        neutral_diploid_markers_[NUM_NEUTRAL_DIPLOID_LOCII];  // track the genealogy of neutral genes in this organism        
        Organism::Sex       sex_;                                           // male or female
        float               fitness_;                                       // cache this organism's fitness
        bool                expired_;                                       // flag an organism to be removed allowing for use of std::remove_if() and std::resize() or v.erase()
};
// Organism

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//! A collection of processes and properties that determine the ecologies of
//! organisms.
class Species {

    public:
    
        // --- lifecycle and assignment ---        
        Species(unsigned index,
                const char* label, 
                unsigned num_fitness_factors,
                RandomNumberGenerator& rng);              
        ~Species() {}        
                
        // --- access and mutation ---
        unsigned get_index() const {
            return this->index_;
        }
        void set_num_fitness_factors(unsigned i) {
            this->num_fitness_factors_ = i;
        }
        unsigned get_num_fitness_factors() const {
            return this->num_fitness_factors_;
        }
        void set_index(unsigned i) {
            this->index_ = i;
        }            
        float get_mutation_rate() const {
            return this->mutation_rate_;
        }
        void set_mutation_rate(float i) {
            this->mutation_rate_ = i;
        }
        FitnessFactorType get_max_mutation_size() const {
            return this->max_mutation_size_;
        }
        void set_max_mutation_size(FitnessFactorType i) {
            this->max_mutation_size_ = i;
        }
        unsigned get_mean_reproductive_rate() const {
            return this->mean_reproductive_rate_;
        }
        void set_mean_reproductive_rate(unsigned i) {
            this->mean_reproductive_rate_ = i;
        }     
        unsigned get_reproductive_rate_mutation_size() const {
            return this->reproductive_rate_mutation_size_;
        }
        void set_reproductive_rate_mutation_size(unsigned i) {
            this->reproductive_rate_mutation_size_ = i;
        }
        int get_movement_capacity() const {
            return this->movement_capacity_;
        }
        void set_movement_capacity(int i) {
            this->movement_capacity_ = i;
        }   
        
        void set_movement_costs(const std::vector<int>& costs) {
            this->movement_costs_ = costs;
        }
        int movement_cost(CellIndexType i) {
            assert( (i >= 0 ) and (static_cast<unsigned>(i) < this->movement_costs_.size()) );
            return this->movement_costs_[i];
        }
        
        void set_selection_strengths(const std::vector<float>& strengths) {
            this->selection_strengths_ = strengths;
        }
        void set_default_genotype(const FitnessFactors& genotype) {
            memcpy(this->default_genotypic_fitness_factors_, genotype, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));
        }         
        
        // --- fitness ---
        float calc_fitness(const Organism& organism, const FitnessFactors environment) const {
            const FitnessFactors& genotype = organism.genotype();        
            const FitnessFactorType * g = genotype;
            const FitnessFactorType * e = environment;
            std::vector<float>::const_iterator s = this->selection_strengths_.begin();
            float weighted_distance = 0.0;
            for (unsigned i = 0; i < this->num_fitness_factors_; ++i, ++g, ++e, ++s) {
                weighted_distance += pow((*e - *g), 2) * *s; // each distance weighted by selection strength
            }
            return exp(-weighted_distance);
        }                    
        
        // --- organism generation and reproduction ---
        void new_generation();
        
        std::string new_organism_label() {
            std::ostringstream label_ostr;
            label_ostr << this->label_ << "_" << this->organism_counter_++;
            return label_ostr.str();
        }
        
        Organism::Sex get_random_sex(float female_threshold=0.5) const {
            if (this->rng_.uniform_real() < female_threshold) {
                return Organism::Male;
            } else {
                return Organism::Female;
            }
        }  
        
        Organism new_organism() {
            return Organism(this->index_, 
                            this->new_organism_label(), 
                            this->default_genotypic_fitness_factors_, 
                            this->get_random_sex());
        }
        
        Organism new_organism(const Organism& female, const Organism& male) {
            Organism organism(this->index_, this->new_organism_label(), this->get_random_sex());
            // // std::cout << "\n\nCreating new organism: " << &organism << std::endl;
            organism.inherit_genotypic_fitness_factors(female, 
                                                       male, 
                                                       this->num_fitness_factors_, 
                                                       this->mutation_rate_, 
                                                       this->max_mutation_size_, 
                                                       this->rng_);
            organism.inherit_genealogies(female, male, this->rng_);
            return organism;
        }        
        
    private:
        // declared as private (and undefined) to prevent copying/assignment
        Species(const Species& species);    
        const Species& operator=(const Species& species);
        
    private:        
        unsigned                    index_;                     // "slot" in cell's pop vector    
        std::string                 label_;                     // arbitrary identifier
        unsigned                    num_fitness_factors_;       // so genotypes of appropriate length can be composed                
        std::vector<float>          selection_strengths_;       // weighted_distance = distance / (sel. strength)
        float                       mutation_rate_;             // rate of mutations
        FitnessFactorType           max_mutation_size_;         // window "size" of mutations
        unsigned                    mean_reproductive_rate_;    // "base" reproductive rate
        unsigned                    reproductive_rate_mutation_size_;  // if reprod. rate evolves, size of step
        std::vector<int>            movement_costs_;            // the movement surface: the "cost" to enter into every cell on the landscape
        int                         movement_capacity_;         // number of cells per round an individual can move
        FitnessFactors              default_genotypic_fitness_factors_;          // genotype of individuals generated de novo        
        RandomNumberGenerator&      rng_;                       // rng to use
        unsigned long               organism_counter_;

};
// Species
///////////////////////////////////////////////////////////////////////////////

} // gingko namespace

#endif
