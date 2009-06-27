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

#include "biosys.hpp"
#include "landscape.hpp"
#include <functional>
#include <iterator>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <stdexcept>

namespace gingko {

/**
 * General i/o error.
 */
class TreeIOError : public std::runtime_error {
    public:
        TreeIOError(const char * msg) : std::runtime_error(msg) {}
        TreeIOError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

/**
 * General format error.
 */
class TreeStructureError : public std::runtime_error {
    public:
        TreeStructureError(const char * msg) : std::runtime_error(msg) {}
        TreeStructureError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

/**
 * No root on tree.
 */
class TreeStructureMissingRootError : public TreeStructureError {
    public:
        TreeStructureMissingRootError(const char * msg) : TreeStructureError(msg) {}
        TreeStructureMissingRootError(const std::string& msg) : TreeStructureError(msg.c_str()) {}
};

/**
 * No root on tree.
 */
class TreeStructureMultipleRootError : public TreeStructureError {
    public:
        TreeStructureMultipleRootError(const char * msg) : TreeStructureError(msg) {}
        TreeStructureMultipleRootError(const std::string& msg) : TreeStructureError(msg.c_str()) {}
};

/**
 * Tracks a subset of paths on a tree.
 */
class Path {

    public:
        Path(unsigned long index,
             std::vector<Path>* paths_ptr,
             std::map<GenealogyNode *, std::string>* node_labels_map_ptr,
             std::map<GenealogyNode *, unsigned long>* node_path_index_map_ptr,
             std::map<GenealogyNode *, CellIndexType>* node_cell_indexes_map_ptr,
             Landscape * landscape_ptr);
        void add_node(GenealogyNode * node);
        void split_after_node(GenealogyNode * node, 
                              unsigned long new_child_path_idx, 
                              unsigned long other_child_path_idx);
        unsigned long size() {
            return this->path_nodes_.size();
        }
        GenealogyNode * get_tip_node() {
            return this->path_nodes_.front();
        }
        unsigned long get_index() {
            return this->index_;
        }
        void write_newick(std::ostream& out);

    private:
        /** Index of this path in master path pool. */
        unsigned long                               index_;
        /** Pointer to master path pool. */
        std::vector<Path>*                          paths_ptr_;
        /** Tracks nodes within this path. */
        std::vector<GenealogyNode *>                path_nodes_;
        /** Tracks paths branching off this one. */
        std::vector<unsigned long>                  child_path_indexes_;
        /** Maps node indexes to their corresponding label. */
        std::map<GenealogyNode *, std::string>*     node_labels_map_ptr_;
        /** Tracks path of nodes. */
        std::map<GenealogyNode *, unsigned long>*   node_path_index_map_ptr_;        
        /** Tracks indexes of cells. */
        std::map<GenealogyNode *, CellIndexType>*   node_cell_indexes_map_ptr_;          
        /** Landscape from which the nodes are derived. */
        Landscape *                                 landscape_ptr_;      
        
};

/**
 * Encapsulates the building of trees from a collection of GenealogyNode 
 * objects.
 */
class Tree {

    public:
    
        /**
         * Constructor.
         * @param landscape               needed to convert cell indices to x, 
         *                                y coordinates.
         */
        Tree(Landscape* landscape_ptr=NULL);
        
        /**
         * Adds a node (and its lineage to the original ancestor) to the tree 
         * if is not already there, and returns index of the node in the 
         * parent array structure representing the tree.
         *
         * @param  node     pointer to GenealogyNode to be added to the tree
         * @param  label    pointer to std::string reprensting the label for 
         *                  this node (required if it is a leaf node)
         * @return          index of the node in the parent array (or -1 if node 
         *                  was a null node, as would be the case if a node 
         *                  without a parent passed its parent pointer to be 
         *                  inserted into the array)
         */
        void add_leaf(GenealogyNode* node, const std::string * label=NULL);

        /**
         * Writes newick string representing the tree structure to the given
         * output stream.
         * 
         * @param out   output stream to which to write the tree
         */
        void write_newick_tree(std::ostream& out);
        
        Path * add_new_path() {
            this->paths_.push_back(Path(this->paths_.size(), 
                                        &this->paths_, 
                                        &this->node_labels_map_, 
                                        &this->node_path_index_map_, 
                                        &this->node_cell_indexes_map_, 
                                        this->landscape_ptr_));
            return &this->paths_.back();
        }

    private:    
        /** Maps node indexes to their corresponding label. */
        std::map<GenealogyNode *, std::string>      node_labels_map_;
        /** Tracks path of nodes. */
        std::map<GenealogyNode *, unsigned long>    node_path_index_map_;
        /** Tracks indexes of nodes. */
        std::map<GenealogyNode *, CellIndexType>    node_cell_indexes_map_;        
        /** Landscape from which the nodes are derived. */
        Landscape *                                 landscape_ptr_;       
        /** Paths, with [0] being main path */        
        std::vector<Path>                           paths_;
        
};


} // namespace gingko

#endif 
