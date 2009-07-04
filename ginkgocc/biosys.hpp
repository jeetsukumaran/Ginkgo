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

#if !defined(GINKGO_BIOSYS_H)
#define GINKGO_BIOSYS_H

/** @file biosys.h */

#include <algorithm>
#include <cassert>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>

#include "ginkgo_defs.hpp"
#include "randgen.hpp"

#if defined(MEMCHECK)
#include "memcheck.hpp"
#endif

namespace ginkgo {

class Organism;
class Species;

/** collection of Species Pointers */
typedef std::map<std::string, Species *> SpeciesByLabel;

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
        GenealogyNode * random_node(RandomNumberGenerator rng) const {
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
         * @param species           pointer to Species of this Organism
         * @param genotype          array of fitness factor values representing
         *                          the inheritable portion of fitness of this
         *                          individual organism
         * @param sex               gender of this organism
         */ 
        Organism(Species * species,
                 const FitnessFactors& new_genotype, 
                 Organism::Sex new_sex) 
                : species_(species),
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
         * @param sex               gender of this organism
         */        
        Organism(Species * species,      
                 Organism::Sex new_sex) 
                : species_(species),              
                  sex_(new_sex),
                  fitness_(-1),
                  expired_(false) {                 
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
            memcpy(this->genotypic_fitness_factors_, ind.genotypic_fitness_factors_, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));                    
            this->sex_ = ind.sex_;
            this->fitness_ = ind.fitness_;
            this->expired_ = ind.expired_;
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
            
            // genotypic_fitness_factors            
            FitnessFactors t_genotypic_fitness_factors;
            memcpy(t_genotypic_fitness_factors, this->genotypic_fitness_factors_, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));
            memcpy(this->genotypic_fitness_factors_, o2.genotypic_fitness_factors_, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));
            memcpy(o2.genotypic_fitness_factors_, t_genotypic_fitness_factors, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));            

            // sex
            Organism::Sex t_sex = this->sex_;
            this->sex_ = o2.sex_;
            o2.sex_ = t_sex;            
            
            // fitness
            float t_fitness = this->fitness_;
            this->fitness_ = o2.fitness_;
            o2.fitness_ = t_fitness;
            
            // expired
            bool t_expired = this->expired_;
            this->expired_ = o2.expired_;
            o2.expired_ = t_expired;            
            
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
                if (rng.uniform_01() < mutation_rate) {
                    ff_value += rng.uniform_int(-max_mutation_size, max_mutation_size);
                }
                this->genotypic_fitness_factors_[i] = ff_value;
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
        
    private:
    
        /** 
         * pointer to Species object in the World species pool of the
         * species of this organism
         */
        Species *               species_;
        
        /** values of inheritable component of fitness */
        FitnessFactors          genotypic_fitness_factors_;
        
        /** genealogy inherited through female alone */
        HaploidMarker           neutral_haploid_marker_;
        
        /** diploid genealogies */
        DiploidMarker           neutral_diploid_markers_[NUM_NEUTRAL_DIPLOID_loci];
        
        /** for reproduction */
        Organism::Sex           sex_;
        
        /** to cache pre-calculated fitness values */       
        float                   fitness_;
        
        /** flag organism to be removed */        
        bool                    expired_;
};
// Organism
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// CompareOrganismFitness
/**
 * Serves as function object in Organism sorting / set insertion etc.
 */
class CompareOrganismFitness {

    public:
        CompareOrganismFitness() {}
        bool operator()(Organism * o1, Organism * o2) {
            assert(o1->get_fitness() >= 0);
            assert(o2->get_fitness() >= 0);
            return (o1->get_fitness() < o2->get_fitness());
        }
};
// CompareOrganismFitness
///////////////////////////////////////////////////////////////////////////////

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
         * @param label                 identifier string for this species
         * @param num_fitness_factors   number of factors in fitness function
         * @param rng                   random number generator
         */
        Species(const std::string& label, 
                unsigned num_fitness_factors,
                unsigned fitness_factor_grain,
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
        unsigned get_num_fitness_factors() const {
            return this->num_fitness_factors_;
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
        void set_num_fitness_factors(unsigned num_fitness_factors) {
            this->num_fitness_factors_ = num_fitness_factors;
        }
        
        /**
         * Returns fitness factor grain.
         *
         * Returns the scaling factor for fitness calculation.
         *
         * @return number of fitness factor grain
         */        
        unsigned get_fitness_factor_grain() const {
            return this->fitness_factor_grain_;
        }        

        /**
         * Returns the probability of mutation when inheriting an organism 
         * inherits a genotypic fitness factor from its parent.
         *
         * @return rate of mutation per factor inheritance 
         */           
        float get_mutation_rate() const {
            return this->mutation_rate_;
        }
        
        /**
         * Sets the probability of mutation when inheriting an organism 
         * inherits a genotypic fitness factor from its parent.
         *
         * @param rate of mutation per factor inheritance 
         */           
        void set_mutation_rate(float mutation_rate) {
            this->mutation_rate_ = mutation_rate;
        }
        
        /**
         * Returns the window size that a genotypic fitness factor can be 
         * changed when inheriting a genotypic fitness factor from its parent.
         *
         * @return size of mutation (+/-)
         */           
        FitnessFactorType get_max_mutation_size() const {
            return this->max_mutation_size_;
        }
        
        /**
         * Sets the window size that a genotypic fitness factor can be 
         * changed when inheriting a genotypic fitness factor from its parent.
         *
         * @param i size of mutation (+/-)
         */          
        void set_max_mutation_size(FitnessFactorType mutation_size) {
            this->max_mutation_size_ = mutation_size;
        }
        
        FitnessFactorType get_max_mutation_size() {
            return this->max_mutation_size_;
        }
        
        /**
         * Returns the mean number of offspring per female.
         *
         * @return mean number of offspring per female
         */         
        unsigned get_mean_reproductive_rate() const {
            return this->mean_reproductive_rate_;
        }
        
        /**
         * Sets the mean number of offspring per female.
         *
         * @param mean number of offspring per female
         */           
        void set_mean_reproductive_rate(unsigned rate) {
            this->mean_reproductive_rate_ = rate;
        }     
        
        /**
         * Returns the movement currency or potential for an organism.
         *
         * @return the movement potential of an organism
         */          
        int get_movement_capacity() const {
            return this->movement_capacity_;
        }
        
        /**
         * Sets the movement currency or potential for an organism.
         *
         * @param moves the movement potential of an organism
         */          
        void set_movement_capacity(int moves) {
            this->movement_capacity_ = moves;
        }  
        
        /**
         * Sets the probability of movement of an organism.
         *
         * @param movement_probability probability of movement
         */
        void set_movement_probability(float movement_probability) {
            this->movement_probability_ = movement_probability;
        }
                
        /**
         * Gets the probability of movement of an organism.
         *
         * @param movement_probability probability of movement
         */
        float get_movement_probability() {
            return this->movement_probability_;
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
        void set_movement_costs(const std::vector<long>& costs) {
            this->movement_costs_ = costs;
        }
        
        /**
         * Returns the cost to be paid by an organism to enter a particular 
         * cell.
         * @param  i    index of the cell to be entered
         * @return      value to be deducted from organism's movement capacity
         */
        int movement_cost(CellIndexType i) {
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
            this->selection_weights_ = strengths;
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
        void set_default_genotypic_fitness_factors(const FitnessFactors& genotype) {
            memcpy(this->default_genotypic_fitness_factors_, genotype, MAX_FITNESS_FACTORS*sizeof(FitnessFactorType));
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
        void set_default_genotypic_fitness_factors(std::vector<FitnessFactorType> genotype) {
            assert(genotype.size() <= MAX_FITNESS_FACTORS);
            for (unsigned i = 0; i != genotype.size(); ++i) {
                this->default_genotypic_fitness_factors_[i] = genotype[i];
            }
        }
        
        std::vector<FitnessFactorType> get_default_genotypic_fitness_factors() {
            std::vector<FitnessFactorType> v;
            v.reserve(this->get_num_fitness_factors());
            for (unsigned i = 0; i < this->get_num_fitness_factors(); ++i) {
                v.push_back(this->default_genotypic_fitness_factors_[i]);    
            }
            return v;
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
        float calc_fitness(const Organism& organism, const FitnessFactors environment) const {
            const FitnessFactors& genotype = organism.genotype();        
            const FitnessFactorType * g = genotype;
            const FitnessFactorType * e = environment;
            std::vector<float>::const_iterator s = this->selection_weights_.begin();
            float weighted_distance = 0.0;
            float diff = 0.0;
            for (unsigned i = 0; i < this->num_fitness_factors_; ++i, ++g, ++e, ++s) {
//                 weighted_distance += pow((*e - *g), 2) * *s; // each distance weighted by selection strength
                diff = static_cast<float>(*e)/this->fitness_factor_grain_ - static_cast<float>(*g)/this->fitness_factor_grain_;
                weighted_distance += (diff * diff) * *s; // each distance weighted by selection strength
            }
            return exp(-weighted_distance);
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
        std::string compose_organism_label(CellIndexType x, CellIndexType y) {
            std::ostringstream label_ostr;
            label_ostr << this->label_ << "_x" << x << "_y" << y << "_" << this->organism_label_index_++;
            return label_ostr.str();
        }        
        
        /**
         * Returns a label for the given organism.
         *
         * @param   organism    organism to be labelled
         * @return              unique label for organism
         */
         const std::string& get_organism_label(const Organism& organism) {
            std::map<const Organism *, std::string>::const_iterator ol = this->organism_labels_.find(&organism);
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
         * @param   organism    organism to be labelled
         * @param   x           x-coordinate of organism's location
         * @param   y           y-coordinate of organism's location
         * @return              unique label for organism
         */
         const std::string& set_organism_label(const Organism& organism, CellIndexType x, CellIndexType y) {
            this->organism_labels_[&organism] = this->compose_organism_label(x, y);
            std::map<const Organism *, std::string>::const_iterator ol = this->organism_labels_.find(&organism);
            return ol->second;
         }         
         
        /**
         * Erases a label assigned to an organism.
         *
         * @param       organism label mapping to be erased
         */
         void erase_organism_label(const Organism& organism) {
            this->organism_labels_.erase(&organism);
         }         
        
        /**
         * Returns male/female with uniform probability.
         *
         * @return      Organism::Male or Organism::Female
         */
        Organism::Sex get_random_sex(float female_threshold=0.5) const {
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
        Organism new_organism(CellIndexType cell_index) {
            Organism    org(this, 
                            this->default_genotypic_fitness_factors_, 
                            this->get_random_sex());
            org.set_cell_index(cell_index);     
            return org;
        }
        
        /**
         * Returns an offspring produced by two organisms.
         *
         * Genotypic fitness factors of the offspring are random mosaic 
         * composed of factors selected with uniform random probability from
         * either male or female genotype of corresponding elements, with 
         * random mutation of values. Genealogies are linked in with parental
         * nodes as ancestors: haploid node of offspring takes female's haploid
         * node as ancestor, whereas in diploid loci, each of the two nodes at
         * each loci selects one node at random from the male, and one from 
         * the female parent. Sex is assigned randomly.
         *
         * @param female    female parent
         * @param male      male parent
         * @return          offspring 
         */
        Organism new_organism(const Organism& female, const Organism& male, CellIndexType cell_index) {
            Organism organism(this, this->get_random_sex());
            organism.inherit_genotypic_fitness_factors(female, 
                                                       male, 
                                                       this->num_fitness_factors_, 
                                                       this->mutation_rate_, 
                                                       this->max_mutation_size_, 
                                                       this->rng_);
            organism.inherit_genealogies(female, male, this->rng_, cell_index);
            return organism;
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
        unsigned                            num_fitness_factors_;
        /** scaling factor for fitness coefficients, for finer resolution of fitness levels */
        unsigned int                        fitness_factor_grain_;
        /** species-level weighting coefficients for the fitness functions */
        std::vector<float>                  selection_weights_;
        /** rate of mutation for the genotypic fitness factors */
        float                               mutation_rate_;
        /** window for perturbations of fitness factor values */
        FitnessFactorType                   max_mutation_size_;
        /** mean number of offspring per female */
        unsigned                            mean_reproductive_rate_;
        /** allowing for evolution in fecundity */
        unsigned                            reproductive_rate_mutation_size_;
        /** landscape migration potential for this species */
        std::vector<long>                   movement_costs_;
        /** movement potential of each organism at the start of each round */
        int                                 movement_capacity_;
        /** probability of movement of each organism */
        float                               movement_probability_;
        /** genotype for organisms created de novo */
        FitnessFactors                      default_genotypic_fitness_factors_;
        /** source of random numbers of various distributions */
        RandomNumberGenerator&              rng_;
        /** tracks the number of labels assigned to organisms of this species */
        unsigned long                       organism_label_index_;
        /** tracks the organism to label assignment */
        std::map<const Organism *, std::string>   organism_labels_;

};
// Species
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