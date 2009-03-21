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

/**
 * Encapsulates the building of trees from a collection of GenealogyNode 
 * objects.
 */
class Tree {

    typedef std::map< long, std::string >           NodeIndexToLabelMap;
    typedef std::map< GenealogyNode*, long >        NodeToIndexMap;
    typedef std::vector<GenealogyNode *>            NodeVector;
    typedef std::vector<long>                       ParentIndexVector;

    public:
    
        /**
         * Constructor.
         *
         * @param coalesce_multiple_roots if <code>false</code> throws
         *                                exception if nodes do not coalesce
         *                                into single ancestor
         */
        Tree(bool coalesce_multiple_roots=true);
        
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
        long process_node(GenealogyNode* node, const std::string * label=NULL);
        
        /**
         * Given an index of a node in the parent array, returns the indexes
         * of all its children.
         *
         * @param   parent  index of a node in the parent array
         * @returns         vector of indexes of all nodes with this node as
         *                  parent
         */
        std::vector<long> get_children(long parent);
        
        /**
         * Returns label for given node.
         *
         * @param   node_idx    index of node
         * @return              label for node
         */
         const std::string& get_label_for_node(long node_idx);
         
        /**
         * Composes newick string representation of node relationships given
         * in the parent array. Multiple trees are returned if the nodes have 
         * not coalesced.
         * 
         * @param out   output stream to which to write the tree
         */
//         std::vector<std::string> Tree::compose_newick_tree()         
        
        /**
         * Writes newick string representing the tree structure to the given
         * output stream.
         * 
         * @param out   output stream to which to write the tree
         */
        void write_newick_tree(std::ostream& out);
        
        /**
         * Writes the newick representation of a single node specified by its
         * index the parent array structure to the given output stream.
         *
         * @param node_idx  index of node in the parent array structure
         * @param out       output stream to which to write the newick string
         */
        void write_newick_node(long node_idx, std::ostream& out);
        
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

        /**
         * Writes out tree structure in human readable format for debugging.
         * @param out   output stream to write to
         */
        void dump(std::ostream& out);
        
        /**
         * Directly adds an element into the parent index array.
         *
         * @param parent_idx    index of parent of the node
         * @param label         label of the node
         */
        void add_indexed_node(long parent_index, const char * label = NULL);

    private:
        /** Maps node pointers to indexes of the corresponding node in the parent array */
        NodeToIndexMap                      node_indexes_;
        /** 
         * Parent array structure tracking nodes in a tree, where the value
         * at each location in the array is the index of the parent of that
         * node.
         */
        ParentIndexVector                   tree_nodes_;
        /** Maps node indexes to their corresponding label. */
        NodeIndexToLabelMap                 labels_;
        /** True if multiple roots are to be coalesced into a single node. */
        bool                                coalesce_multiple_roots_;
};


} // namespace gingko

#endif 
