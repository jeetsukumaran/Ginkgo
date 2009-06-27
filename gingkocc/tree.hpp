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
        Path();
        ~Path();
        void add_node(GenealogyNode * node);
        long find_node(GenealogyNode * node);
        Path split_on_index(long idx);

    private:
        std::vector<GenealogyNode *>                path_nodes_;
        std::map<GenealogyNode *, long>             node_to_indexes_;
        Path *                                      first_child_;
        Path *                                      next_sib_;
};
 

/**
 * Encapsulates the building of trees from a collection of GenealogyNode 
 * objects.
 */
class Tree {

    public:
    
        /**
         * Constructor.
         *
         * @param coalesce_multiple_roots if <code>false</code> throws
         *                                exception if nodes do not coalesce
         *                                into single ancestor
         * @param landscape               needed to convert cell indices to x, 
         *                                y coordinates.
         */
        Tree(Landscape* landscape_ptr=NULL, bool coalesce_multiple_roots=true);
        
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
        void process_node(GenealogyNode* node, const std::string * label=NULL);

        /**
         * Writes newick string representing the tree structure to the given
         * output stream.
         * 
         * @param out   output stream to which to write the tree
         */
        void write_newick_tree(std::ostream& out);
        
        /**
         * Returns the current mode of treating multiple roots (as errors
         * or as coalescence in infinity).
         *
         * @return  <code>true</code> if multiple roots get automatically 
         *          coalesced
         */
        bool get_coalesce_multiple_roots() const;
        
        /**
         * Sets the current mode of treating multiple roots (as errors
         * or as coalescence in infinity).
         *
         * @param val   <code>true</code> if multiple roots get automatically 
         *              coalesced
         */        
        void set_coalesce_multiple_roots(bool val);


    private:
        /** List of paths. */        
        std::vector<Path>                           paths_;        
        /** Primary path. */        
        Path *                                      start_path_;  
        /** Maps node indexes to their corresponding label. */
        std::map<GenealogyNode *, std::string>   labels_;
        /** True if multiple roots are to be coalesced into a single node. */
        bool                                        coalesce_multiple_roots_;
        /** Landscape from which the nodes are derived. */
        Landscape *                                 landscape_ptr_;
        

};


} // namespace gingko

#endif 
