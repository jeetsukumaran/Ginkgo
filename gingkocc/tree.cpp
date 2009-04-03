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

Tree::Tree(bool coalesce_multiple_roots) 
    : coalesce_multiple_roots_(coalesce_multiple_roots) {
}

long Tree::process_node(GenealogyNode* node, const std::string * label) {
    if (node == NULL) {
        // Before a node can be inserted, it needs the index of its 
        // parent; 
        // Thus a node typically calls <code>process_node</code>
        // on its parent pointer to obtain this index before inserting
        // itself. 
        // If a node has no parent, then the pointer passed here will 
        // be NULL, and so a flag is returned indicating that the node
        // is a root.
        return -1;
    }
    
    // search for the node in the node to parent array index map
    NodeToIndexMap::iterator nidx = this->node_indexes_.find(node);                
    if (nidx != this->node_indexes_.end()) {
        // if found, then return its index
        return nidx->second;
    }
    
    // This node is not already in the parent array, and needs to be
    // added.
    // We need to get the index of this node's parent first so we
    // know what value to store as this node's parent.
    // So we call <code>process_node</code> again, passing it the 
    // pointer to this node's parent: if this node's parent is not 
    // already in the array, it will then be added and its new index 
    // returned, otherwise the index of the existing parent node in the 
    // array will be returned. 
    // This index will be inserted into the array in the position 
    // corresponding to this node.
    // If this node's parent pointer is NULL, then -1 will be returned
    // by the <code>process_node</code> and stored in this node's
    // location in lieu of an index, indicating that this node is a 
    // root.
    this->tree_nodes_.push_back(this->process_node(node->get_parent())); 
    
    // record this node's index in the node to array index map
    unsigned long idx = this->tree_nodes_.size() - 1;
    this->node_indexes_.insert(std::make_pair(node, idx));
    
    // and, if a label was passed, store it in the node index to label
    // map
    if (label != NULL) {
        this->labels_.insert(std::make_pair(idx, *label));
    } 
//     else {
//         assert(node->get_first_child() != NULL);
//     }
    return idx;
}

std::vector<long> Tree::get_children(long parent) {
    std::vector<long> children;
    for (unsigned long i = parent+1; i < this->tree_nodes_.size(); ++i) {
        if (this->tree_nodes_[i] == parent) {
            children.push_back(i);
        }
    }
    return children;
}

const std::string& Tree::get_label_for_node(long node_idx) {
    NodeIndexToLabelMap::iterator node_label = this->labels_.find(node_idx);
    assert(node_label != this->labels_.end());
    return node_label->second;
}

// std::vector<std::string> Tree::compose_newick_tree() {
//     std::vector<std::string> trees_as_newick;
//     std::map<long, std::string> nodes_as_newick;
//     for (long node_idx = this->tree_nodes_.size()-1; node_idx >= 0; --node_idx) {
//         long parent_idx = this->tree_nodes_[node_idx];
//         if (parent_idx == -1) {
//             trees_as_newick.push_back("(" + nodes_as_newick[node_idx] + ")");
//         } else {
//             std::string& parent_newick_string = nodes_as_newick[parent_idx];
//             if (parent_newick_string.size() > 0) {
//                 parent_newick_string += ",";
//             }        
//             if (nodes_as_newick[node_idx].size() == 0) {
//                 // leaf node: add label to parent
//                 parent_newick_string += this->get_label_for_node(node_idx);
//             } else {
//                 // internal node: wrap in parentheses and then write to parent
//                 // PROBLEM: Nodes of outdegree 1!
//                 parent_newick_string += "(" + nodes_as_newick[node_idx] +")";
//             }        
//         }
//         nodes_as_newick.erase(node_idx);
//     }
// }

void Tree::write_newick_tree(std::ostream& out) {
    int num_roots = std::count(this->tree_nodes_.begin(), this->tree_nodes_.end(), -1);         
    if (num_roots == 0) {
        throw TreeStructureMissingRootError("no root nodes found (possibly because node list was empty)");
    }
    if (this->coalesce_multiple_roots_ and num_roots > 1) {
        ParentIndexVector::iterator root = std::find(this->tree_nodes_.begin(), this->tree_nodes_.end(), -1);
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
        if (num_roots >= 2)  {
            throw TreeStructureMultipleRootError("multiple roots found");
        }
        ParentIndexVector::iterator root = std::find(this->tree_nodes_.begin(),
                this->tree_nodes_.end(), -1);
        if (root == this->tree_nodes_.end())  {
            throw TreeStructureMissingRootError("no root nodes found (possibly because node list was empty)");
        }
        this->write_newick_node(root-this->tree_nodes_.begin(), out);
        out << ";";
    }                
}

void Tree::write_newick_node(long node_idx, std::ostream& out) {
    unsigned edge_length = 1; 
    std::vector<long> children = this->get_children(node_idx);
    while (children.size() == 1) {
        // this deals with nodes of outdegree 1 still in the structure
        ++edge_length;
        long only_child = children[0];
        children = this->get_children(only_child);
        if (children.size() == 0) {
            // node has chain of outdegree 1 descendents all the way
            // to single terminal: collapse to child
            node_idx = only_child;
        }
    }
    if (children.size() > 0) {
        out << "(";
        for (std::vector<long>::iterator child_iter = children.begin();
             child_iter != children.end();
             ++child_iter) {
             if (child_iter != children.begin()) {
                out << ",";
             }
            this->write_newick_node(*child_iter, out);    
        }
        out << ")";
    } else {
        NodeIndexToLabelMap::iterator node_label = this->labels_.find(node_idx);
        assert(node_label != this->labels_.end());
        out << node_label->second;
    }
    out << ":" << edge_length;
}

bool Tree::get_coalesce_multiple_roots() const {
    return this->coalesce_multiple_roots_;
}

void Tree::set_coalesce_multiple_roots(bool val) {
    this->coalesce_multiple_roots_ = val;
}

void Tree::add_indexed_node(long parent_index, const char * label) {
    this->tree_nodes_.push_back(parent_index);
    if (label != NULL) {
        std::string label_str(label);
        unsigned long idx = this->tree_nodes_.size() - 1;
        this->labels_.insert(std::make_pair(idx, label_str));
    }
}

void Tree::dump(std::ostream& out) {
    out << std::setw(10) << "idx" << "   ";
    out << std::setw(10) << "parent" << "   ";
    out << std::setw(10) << "label" << "\n";        
    for (unsigned i = 0; i < this->tree_nodes_.size(); ++i) {
        out << std::setw(10) << i << "   ";
        out << std::setw(10) << this->tree_nodes_[i] << "   ";
        out << std::setw(10) << this->labels_[i] << "\n";
    }
}        

} // namespace gingko
