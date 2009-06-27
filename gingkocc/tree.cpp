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

Path::Path() { }
      
Path::~Path() { }

void Path::add_node(GenealogyNode * node) {
    this->path_nodes_.push_back(node);
    this->node_to_indexes_.insert(std::make_pair(node, this->path_nodes_.size()-1));
}

long Path::get_node_index(GenealogyNode * node) {
//     std::map<GenealogyNode *, long>::iterator i = this->node_to_indexes_.find(node);
//     if ( i == this->node_to_indexes_.end() ) {
//         return -1;
//     } else {
//         return i->second;
//     }
    for (std::vector<GenealogyNode *>::iterator i = this->path_nodes_.begin(); i != this->path_nodes_.end(); ++i) {
        if (*i == node) {
            return i-this->path_nodes_.begin();
        }
    }
    return -1;
}

Path * Path::find_node(GenealogyNode * node, long& idx) {
    idx = this->get_node_index(node);
    if (idx >= 0) {
        return this;
    } else {
        for (std::vector<Path>::iterator pi = this->child_paths_.begin(); pi != this->child_paths_.end(); ++pi) {
            Path * p = pi->find_node(node, idx);
            if (p != NULL) {
                return p;
            }
        }
    }
    return NULL;
}

Path * Path::split_on_index(long idx, Path& new_child) {
    assert(static_cast<unsigned long>(idx) < this->path_nodes_.size());
    assert(idx >= 0);
    Path split_path;
    for (std::vector<GenealogyNode *>::iterator i = this->path_nodes_.begin(); i < (this->path_nodes_.begin() + idx); ++i) {
        split_path.add_node(*i);
        this->node_to_indexes_.erase(*i);        
    }
    this->path_nodes_.erase(this->path_nodes_.begin(), this->path_nodes_.begin()+idx);
    assert(this->node_to_indexes_.size() == this->path_nodes_.size());    
    split_path.child_paths_.swap(this->child_paths_);
    this->child_paths_.push_back(new_child);    
    this->child_paths_.push_back(split_path);    
    return &this->child_paths_.back();
}

void Path::write_newick(std::ostream& out, std::map<GenealogyNode *, std::string>& labels) {
    if (this->child_paths_.size() == 0) {
        std::map<GenealogyNode *, std::string>::iterator pni = labels.find(this->path_nodes_.front());
        //assert(pni != labels.end());
        if (pni == labels.end()) {
            out << "x:" << this->size();
        } else {
            out << pni->second << ":" << this->size();            
        }
    } else {
        out << "(";
        for (std::vector<Path>::iterator pi = this->child_paths_.begin(); pi != this->child_paths_.end(); ++pi) {
            if (pi != this->child_paths_.begin()) {
                out << ",";
            }
            pi->write_newick(out, labels);
        }
        out << "):" << this->size();
    }
}

///////////////////////////////////////////////////////////////////////////////
// Tree

Tree::Tree(Landscape * landscape_ptr, bool coalesce_multiple_roots) 
    : landscape_ptr_(landscape_ptr),
      coalesce_multiple_roots_(coalesce_multiple_roots) {
}

void Tree::process_node(GenealogyNode* node, const std::string * label) {

    if (label != NULL) {
        this->labels_.insert(std::make_pair(node, *label));
    } 

    if (this->start_path_.size() == 0) {
        this->paths_.push_back(&this->start_path_);
        while (node != NULL) {
            this->start_path_.add_node(node);
            node = node->get_parent();
        }
    } else {
        Path new_node_path;
        bool coalesced = false;
        while (node != NULL) {
            long idx = -1;
            Path * p = this->start_path_.find_node(node, idx);
            if (p != NULL && idx >= 0) {
                p->split_on_index(idx, new_node_path);
                coalesced = true;
                return;
            } else {
                new_node_path.add_node(node);
                node = node->get_parent();
            }            
        }
        assert(coalesced); // TODO: if not coalesced here, multiple roots
    }
}

void Tree::write_newick_tree(std::ostream& out) {
//     Path * p = this->start_path_.find_leaf_path();
//     if (p != NULL) {
//         GenealogyNode * n = p->get_tip_node();
//         std::map<GenealogyNode *, std::string>::iterator pni = this->labels_.find(n);
//         if (pni == this->labels_.end()) {
//             std::cout << "**** hssssere ****" << std::endl << std::endl;
//         }
//     }
    this->start_path_.write_newick(out, this->labels_);
}

bool Tree::get_coalesce_multiple_roots() const {
    return this->coalesce_multiple_roots_;
}

void Tree::set_coalesce_multiple_roots(bool val) {
    this->coalesce_multiple_roots_ = val;
}
       

} // namespace gingko
