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

#include <map>
#include <utility>
#include "biosys.hpp"
#include "tree.hpp"

namespace gingko {

///////////////////////////////////////////////////////////////////////////////
// Path

Path::Path()
    : first_child_(NULL),
      next_sib_(NULL) { }
      
Path::~Path() { }

void Path::add_node(GenealogyNode * node) {
    this->path_nodes_.push_back(node);
    this->node_to_indexes_.insert(std::make_pair(node, this->path_nodes_.size()-1));
}

long Path::find_node(GenealogyNode * node) {
    std::map<GenealogyNode *, long>::iterator i = this->node_to_indexes_.find(node);
    if ( i == this->node_to_indexes_.end() ) {
        return -1;
    } else {
        return i->second;
    }
}

Path Path::split_on_index(long idx) {
    assert(idx < this->path_nodes_.size() && idx >= 0);
    Path p;
    for (long i = this->path_nodes_.size(); i>= idx; --i) {
        p.add_node(this->path_nodes_[i]);
        this->node_to_indexes_.erase(this->path_nodes_[i]);        
    }
    this->path_nodes_.erase(this->path_nodes_.begin()+idx, this->path_nodes_.end());    
    return p;
}

///////////////////////////////////////////////////////////////////////////////
// Tree

Tree::Tree(Landscape * landscape_ptr, bool coalesce_multiple_roots) 
    : landscape_ptr_(landscape_ptr),
      coalesce_multiple_roots_(coalesce_multiple_roots) {
}

void Tree::process_node(GenealogyNode* node, const std::string * label) {

    bool coalesced = false;
    this->paths_.push_back(Path());
    Path& new_node_path = this->paths_.back();
    while (node != NULL and not coalesced) {
        for (std::vector<Path>::iterator pi = this->paths_.begin(); pi != this->paths_.end(); ++pi) {
            long idx = (*pi).find_node(node);
            if (idx >= 0) {
                this->paths_.push_back((*pi).split_on_index(idx));
                coalesced = true;
                break;
            }
        }
        if (not coalesced) {
            new_node_path.add_node(node);
            node = node->get_parent();
        }            
    }
    if (node == NULL) {
        this->start_path_ = &new_node_path;
    }
    if (label != NULL) {
        this->labels_.insert(std::make_pair(node, *label));
    } 
}


void Tree::write_newick_tree(std::ostream& out) {
               
}

bool Tree::get_coalesce_multiple_roots() const {
    return this->coalesce_multiple_roots_;
}

void Tree::set_coalesce_multiple_roots(bool val) {
    this->coalesce_multiple_roots_ = val;
}
       

} // namespace gingko
