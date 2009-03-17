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

#if !defined(GINGKO_TREE_H)
#define GINGKO_TREE_H

#include "biosys.h"
#include <functional>
#include <iterator>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

namespace gingko {

class Tree {

    typedef std::map< GenealogyNode*, long >        NodeIndexMap;
    typedef std::vector<GenealogyNode*>             NodeVector;
    typedef std::vector<long>                       IndexVector;

    public:
        Tree() 
            : coalesce_multiple_roots_(true) {
        }
    
        Tree(std::vector<Organism>& organisms) 
                : coalesce_multiple_roots_(true) {
            std::vector<Organism*> organism_ptrs;
            organism_ptrs.reserve(organisms.size());
            for (std::vector<Organism>::iterator optr = organisms.begin();
                    optr != organisms.end();
                    ++optr) {
                organism_ptrs.push_back(&(*optr));                    
            }
            this->init(organism_ptrs);
        }

        Tree(std::vector<Organism*>& organism_ptrs) 
                : coalesce_multiple_roots_(true) { 
            this->init(organism_ptrs);
        }

        void init(std::vector<Organism*>& organism_ptrs) {
            assert(this->nodes_to_coalesce_.size() == 0);
            this->nodes_to_coalesce_.reserve(organism_ptrs.size());
            for (std::vector<Organism*>::iterator optr_iter = organism_ptrs.begin();
                    optr_iter != organism_ptrs.end();
                    ++optr_iter) {
                Organism& organism = *(*optr_iter);
                // for now: just focus on the haploid markers
                GenealogyNode* gnode = organism.haploid_marker().node();
                this->process_node(gnode);
            }
        }
        
        // adds node if it is note already in the ancestor array
        // returns index in array
        unsigned long process_node(GenealogyNode* node) {
            if (node == NULL) {
                return -1;
            }
//             node->suppress_outdegree1();
            NodeIndexMap::iterator nidx = this->node_indexes_.find(node);
            if (nidx != this->node_indexes_.end()) {
                return nidx->second;
            }            
            this->tree_nodes_.push_back(this->process_node(node->get_parent())); 
            unsigned long idx = this->tree_nodes_.size() - 1;
            this->node_indexes_.insert(std::make_pair(node, idx));
            if (not node->has_label()) {
                std::ostringstream label_os;
                label_os << "K" << std::setw(6) << std::setfill('0') << idx;
                this->labels_.push_back(label_os.str());
            } else {                
                this->labels_.push_back(node->get_label());
            }                
//             this->edge_lens_.push_back(node->get_edge_len());
            return idx;
        }
        
        std::vector<long> get_children(long parent) {
            std::vector<long> children;
//             std::remove_copy_if(this->tree_nodes_.begin(), 
//                                 this->tree_nodes_.end(), 
//                                 std::back_inserter(children), 
//                                 std::bind2nd(std::not_equal_to<long>(), parent));
            for (unsigned long i = 0; i < this->tree_nodes_.size(); ++i) {
                if (this->tree_nodes_[i] == parent) {
                    children.push_back(i);
                }
            }
            return children;
        }                
        
        void write_newick_tree(std::ostream& out) {
            int num_roots = std::count(this->tree_nodes_.begin(), this->tree_nodes_.end(), -1);         
            assert(num_roots > 0);
            if (this->coalesce_multiple_roots_ and num_roots > 1) {
                IndexVector::iterator root = std::find(this->tree_nodes_.begin(), this->tree_nodes_.end(), -1);
                out << "(";
                this->write_newick_node(root-this->tree_nodes_.begin(), out);                
                while (root != this->tree_nodes_.end()) {
                    root = std::find(root+1, this->tree_nodes_.end(), -1);
                    if (root != this->tree_nodes_.end()) {
                        out << ",";
                        this->write_newick_node(root-this->tree_nodes_.begin(), out);
                    }                        
                }
                out << ");"; // add infinite branch length?                
            } else {
                assert(num_roots < 2); 
                IndexVector::iterator root = std::find(this->tree_nodes_.begin(),
                        this->tree_nodes_.end(), -1);
                assert(root != this->tree_nodes_.end());                    
                this->write_newick_node(root-this->tree_nodes_.begin(), out);
                out << ";";
            }                
        }
        
        void write_newick_node(long node_idx, std::ostream& out) {
            unsigned edge_length = 0; 
            std::vector<long> children = this->get_children(node_idx);
            while (children.size() == 1) {
                // this deals with nodes of outdegree 1 still in the structure
                ++edge_length; 
                children = this->get_children(children[0]);
            }
            if (children.size() > 0) {
                out << "(";
                for (std::vector<long>::iterator child_iter = children.begin();
                     child_iter != children.end();
                     ++child_iter) {
                     if (child_iter != children.begin()) {
                        out << ", ";
                     }
                    this->write_newick_node(*child_iter, out);    
                }
                out << ")";
            } else {
                out << this->labels_.at(node_idx);
            }
            out << ":" << edge_length;
        }
        
        bool get_coalesce_multiple_roots() const {
            return this->coalesce_multiple_roots_;
        }
        
        void set_coalesce_multiple_roots(bool val) {
            this->coalesce_multiple_roots_ = val;
        }
                
        void dump(std::ostream& out) {
            out << std::setw(10) << "idx" << "   ";
            out << std::setw(10) << "parent" << "   ";
            out << std::setw(10) << "label" << "\n";        
            for (unsigned i = 0; i < this->tree_nodes_.size(); ++i) {
                out << std::setw(10) << i << "   ";
                out << std::setw(10) << this->tree_nodes_[i] << "   ";
                out << std::setw(10) << this->labels_[i] << "\n";
            }
        }

    private:
        NodeIndexMap                        node_indexes_;
        NodeVector                          nodes_to_coalesce_;
        IndexVector                         tree_nodes_;
        std::vector<std::string>            labels_;
        bool                                coalesce_multiple_roots_;
};


} // namespace gingko

#endif 