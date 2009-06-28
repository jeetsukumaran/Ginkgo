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
#include <list>
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

class Tree;
class Path;

/**
 * Tracks a subset of paths on a tree.
 */
class Path {

    public:
        Path(Tree * tree_ptr);
        void add_node(GenealogyNode * node);
        void add_split_below_node(GenealogyNode * node, Path * other_child_path);
        unsigned long size() {
            return this->path_nodes_.size();
        }
        GenealogyNode * get_tip_node() {
            return this->path_nodes_.front();
        }
        void write_newick(std::ostream& out);

    private:
        /** Tracks nodes within this path. */
        std::vector<GenealogyNode *>                path_nodes_;
        /** Tracks paths branching off this one. */
        std::vector<Path *>                         child_paths_;
        /** Pointer to parent tree. */
        Tree *                                      tree_ptr_;     
        
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
        
        Path * get_node_path(GenealogyNode * node);
        void set_node_path(GenealogyNode * node, Path * path);
        void write_node_cell_xy(GenealogyNode * node, std::ostream& out);
        void store_node_cell_index(GenealogyNode * node);
        Path * add_new_path();
        std::string& get_node_label(GenealogyNode * node);
        
    private:    
        /** Maps node indexes to their corresponding label. */
        std::map<GenealogyNode *, std::string>      node_labels_map_;
        /** Tracks path of nodes. */
        std::map<GenealogyNode *, Path *>           node_path_map_;
        /** Tracks indexes of nodes. */
        std::map<GenealogyNode *, CellIndexType>    node_cell_indexes_map_;        
        /** Landscape from which the nodes are derived. */
        Landscape *                                 landscape_ptr_;       
        /** Paths, with [0] being main path */        
        std::list<Path>                             paths_list_;
        
};


} // namespace gingko

#endif 
