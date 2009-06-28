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

Path::Path(Tree * tree_ptr)
    : tree_ptr_(tree_ptr) { }
      
void Path::add_node(GenealogyNode * node) {
    this->path_nodes_.push_back(node);
    this->tree_ptr_->set_node_path(node, this);
}

void Path::add_split_below_node(GenealogyNode * node, Path * other_child_path) {
    Path * new_child_path = this->tree_ptr_->add_new_path();
    std::vector<GenealogyNode *>::iterator node_loc = std::find(this->path_nodes_.begin(), this->path_nodes_.end(), node);    
    assert(node_loc !=  this->path_nodes_.end());
    new_child_path->path_nodes_.reserve( new_child_path->path_nodes_.size() + (node_loc - this->path_nodes_.begin()) );
    for (std::vector<GenealogyNode *>::iterator i = this->path_nodes_.begin(); i < node_loc; ++i) {
        new_child_path->add_node(*i);
        assert(this->tree_ptr_->get_node_path(*i) == new_child_path);
    }
    this->path_nodes_.erase(this->path_nodes_.begin(), node_loc);   
    new_child_path->child_paths_.swap(this->child_paths_);
    this->child_paths_.push_back(new_child_path);    
    this->child_paths_.push_back(other_child_path);    
}

void Path::write_newick(std::ostream& out) {
    if (this->child_paths_.size() == 0) {
        out << this->tree_ptr_->get_node_label(this->path_nodes_.front()) << ":" << this->size();            
    } else {
        out << "(";
        for (std::vector<Path *>::iterator pi = this->child_paths_.begin(); pi != this->child_paths_.end(); ++pi) {
            if (pi != this->child_paths_.begin()) {
                out << ", ";
            }
            (*pi)->write_newick(out);
        }
        out << ")";
        this->tree_ptr_->write_node_cell_xy(this->path_nodes_.front(), out);
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
    if (this->paths_list_.size() == 0) {
        Path * first_path = this->add_new_path();   
        while (node != NULL) {
            first_path->add_node(node);
            this->store_node_cell_index(node);
            node = node->get_parent();
        }
    } else {
        Path * new_node_path = this->add_new_path();
        bool coalesced = false;
        while (node != NULL) {
            std::map<GenealogyNode *, Path *>::iterator npi = this->node_path_map_.find(node);            
            if (npi != this->node_path_map_.end()) {
                npi->second->add_split_below_node(node, new_node_path);
                coalesced = true;
                return;
            } else {
                new_node_path->add_node(node);
                this->store_node_cell_index(node);
                node = node->get_parent();
            }            
        }
        if (!coalesced) {
            throw TreeStructureMultipleRootError("failure to coalesce to single ancestor in sample");
        }
    }
}

void Tree::write_newick_tree(std::ostream& out) {
    this->paths_list_.front().write_newick(out);
}

Path * Tree::get_node_path(GenealogyNode * node) {
    std::map<GenealogyNode *, Path *>::iterator pi = this->node_path_map_.find(node);
    if (pi == this->node_path_map_.end()) {
        return NULL;
    } else {
        return pi->second;
    }
}

void Tree::set_node_path(GenealogyNode * node, Path * path) {
    this->node_path_map_[node] = path;
}

void Tree::write_node_cell_xy(GenealogyNode * node, std::ostream& out) {
    std::map<GenealogyNode *, CellIndexType>::iterator ci = this->node_cell_indexes_map_.find(node);
    if (ci != this->node_cell_indexes_map_.end()) {
        return;
    } else {
        if (this->landscape_ptr_ != NULL) {
            CellIndexType cell_index = ci->second;
            CellIndexType x = this->landscape_ptr_->index_to_x(cell_index);        
            CellIndexType y = this->landscape_ptr_->index_to_y(cell_index);
            out << "x" << x << "_" << "y" << y;
        } else {
            return;
        }
    }        
}

void Tree::store_node_cell_index(GenealogyNode * node) {
    if (node != NULL) {
        this->node_cell_indexes_map_.insert(std::make_pair(node, node->get_cell_index()));
    }        
}

Path * Tree::add_new_path() {
    this->paths_list_.push_back(Path(this));
    return &this->paths_list_.back();
}

std::string& Tree::get_node_label(GenealogyNode * node) {
    std::map<GenealogyNode *, std::string>::iterator pni = this->node_labels_map_.find(node);
    assert(pni != this->node_labels_map_.end());
    return pni->second;
}        
   

} // namespace gingko
