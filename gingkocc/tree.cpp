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

#include <algorithm>
#include <map>
#include <utility>
#include <vector>
#include "biosys.hpp"
#include "tree.hpp"

namespace gingko {

///////////////////////////////////////////////////////////////////////////////
// Path

Path::Path(unsigned long index,
           std::vector<Path>* paths_ptr,
           std::map<GenealogyNode *, std::string>* node_labels_map_ptr,
           std::map<GenealogyNode *, unsigned long>* node_path_index_map_ptr,
           std::map<GenealogyNode *, CellIndexType>* node_cell_indexes_map_ptr,
           Landscape * landscape_ptr)
    : index_(index),
      paths_ptr_(paths_ptr),
      node_labels_map_ptr_(node_labels_map_ptr),
      node_path_index_map_ptr_(node_path_index_map_ptr),
      node_cell_indexes_map_ptr_(node_cell_indexes_map_ptr),
      landscape_ptr_(landscape_ptr) { }
      
// Path::Path(const Path& p) 
//     : path_nodes_(p.path_nodes_),
//       child_paths_(p.child_paths_),
//       node_labels_map_ptr_(p.node_labels_map_ptr_),
//       node_paths_map_ptr_(p.node_paths_map_ptr_),
//       node_cell_indexes_map_ptr_(p.node_cell_indexes_map_ptr_),
//       landscape_ptr_(p.landscape_ptr_) {
// 
// }
//       
// Path::~Path() { }
// 
// const Path& Path::operator=(const Path& p) {
//     this->path_nodes_ = p.path_nodes_;
//     this->child_paths_ = p.child_paths_;
//     this->node_labels_map_ptr_ = p.node_labels_map_ptr_;
//     this->node_paths_map_ptr_ = p.node_paths_map_ptr_;
//     this->node_cell_indexes_map_ptr_ = p.node_cell_indexes_map_ptr_;
//     this->landscape_ptr_ = p.landscape_ptr_;
//     return *this;
// }

void Path::add_node(GenealogyNode * node) {
    this->path_nodes_.push_back(node);
    if (node != NULL) {
        this->node_cell_indexes_map_ptr_->insert(std::make_pair(node, node->get_cell_index()));
    }
    (*(this->node_path_index_map_ptr_))[node] = this->index_;
}

void Path::split_after_node(GenealogyNode * node, 
        unsigned long new_child_path_idx, 
        unsigned long other_child_path_idx) {
    std::vector<GenealogyNode *>::iterator node_loc = std::find(this->path_nodes_.begin(), this->path_nodes_.end(), node);    
    assert(node_loc !=  this->path_nodes_.end());
    for (std::vector<GenealogyNode *>::iterator i = this->path_nodes_.begin(); i < node_loc; ++i) {
        this->paths_ptr_->at(new_child_path_idx).add_node(*i);
        assert((*(this->node_path_index_map_ptr_))[node] == new_child_path_idx);
    }
    this->path_nodes_.erase(this->path_nodes_.begin(), node_loc);   
    this->paths_ptr_->at(new_child_path_idx).child_path_indexes_.swap(this->child_path_indexes_);
    this->child_path_indexes_.push_back(other_child_path_idx);    
    this->child_path_indexes_.push_back(new_child_path_idx);    
}

void Path::write_newick(std::ostream& out) {
    if (this->child_path_indexes_.size() == 0) {
        std::map<GenealogyNode *, std::string>::iterator pni = this->node_labels_map_ptr_->find(this->path_nodes_.front());
        assert(pni != this->node_labels_map_ptr_->end());
        out << pni->second << ":" << this->size();            
    } else {
        out << "(";
        for (std::vector<unsigned long>::iterator pi = this->child_path_indexes_.begin(); pi != this->child_path_indexes_.end(); ++pi) {
            if (pi != this->child_path_indexes_.begin()) {
                out << ", ";
            }
            (this->paths_ptr_->at(*pi)).write_newick(out);
        }
        out << ")";
        if (this->landscape_ptr_ != NULL) {
            std::map<GenealogyNode *, CellIndexType>::iterator ci = this->node_cell_indexes_map_ptr_->find(this->path_nodes_.front());
            if (ci != this->node_cell_indexes_map_ptr_->end()) {
                CellIndexType cell_index = ci->second;
                CellIndexType x = this->landscape_ptr_->index_to_x(cell_index);        
                CellIndexType y = this->landscape_ptr_->index_to_y(cell_index);
                out << "x" << x << "_" << "y" << y;
            }
        }
        out << ":" << this->size();
    }
}

///////////////////////////////////////////////////////////////////////////////
// Tree

Tree::Tree(Landscape * landscape_ptr) 
    : landscape_ptr_(landscape_ptr) {
}

void Tree::add_leaf(GenealogyNode* node, const std::string * label) {
    if (label != NULL) {
        this->node_labels_map_.insert(std::make_pair(node, *label));
    }
    if (this->paths_.size() == 0) {
        Path * first_path = this->add_new_path();   
        while (node != NULL) {
            first_path->add_node(node);
            node = node->get_parent();
        }
    } else {
        Path * new_node_path = this->add_new_path();
        bool coalesced = false;
        while (node != NULL) {
            std::map<GenealogyNode *, unsigned long>::iterator npi = this->node_path_index_map_.find(node);            
            if (npi != this->node_path_index_map_.end()) {
                Path * new_split_child_path = this->add_new_path();
                this->paths_.at(npi->second).split_after_node(node, 
                                                                new_split_child_path->get_index(), 
                                                                new_node_path->get_index());
                coalesced = true;
                return;
            } else {
                new_node_path->add_node(node);
                node = node->get_parent();
            }            
        }
        assert(coalesced); // TODO: if not coalesced here, multiple roots
    }
}

void Tree::write_newick_tree(std::ostream& out) {
    this->paths_.at(0).write_newick(out);
}  

} // namespace gingko
