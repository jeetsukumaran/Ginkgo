#! /usr/bin/env python

###############################################################################
##
## GINGKO Biogeographical Evolution Simulator.
##
## Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
## This program is free software; you can redistribute it and#or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License along
## with this program. If not, see <http:##www.gnu.org#licenses#>.
##
###############################################################################

import unittest
import subprocess

from gingko_tests import get_logger
from gingko_tests import get_gingko_program_path
from gingko_tests import run_program
from gingko_tests import run_external_tests
import sys
from dendropy import datasets
from dendropy import splits
from dendropy import treedists

_LOG = get_logger("test_tree")

def to_parent_array(tree, include_labels=False, include_edge_lens=False):
    pa = []
    edge_lens = []
    node_pa = {}
    for node in tree.preorder_node_iter():
        node_pa[node] = len(pa)
        if node.parent_node is not None:
            if node.taxon is not None:
                pa.append("%s:'%s'" % (node_pa[node.parent_node], node.taxon.label))
            else:
                pa.append(str(node_pa[node.parent_node]))
        else:
            pa.append(str(-1))
        if include_edge_lens:
            assert node.edge is not None
            edge_lens.append(node.edge.length)
    return pa, edge_lens
    

class ParentArrayTreeRoundTripTest(unittest.TestCase):
    def setUp(self):
        self.prog_path = get_gingko_program_path("test_tree_from_parent_indexes")
        self.src_trees = [
            "(ahli:0.264213,(((garmani:0.106838,grahami:0.086367)1.00:0.069511,valencienni:0.164263)0.87:0.020752,lineatopus:0.195726)1.00:0.077682,(((((((aliniger:0.160001,coelestinus:0.193231)1.00:0.071920,bahorucoensis:0.226688)0.68:0.023043,(equestris:0.022702,luteogularis:0.030641)1.00:0.198165,occultus:0.423120)0.89:0.056277,(barahonae:0.211489,cuvieri:0.168670)1.00:0.084190,(insolitus:0.243882,olssoni:0.256877)1.00:0.050618)0.86:0.031679,((brevirostris:0.180130,distichus:0.115136)1.00:0.123136,((cristatellus:0.214436,krugi:0.157330)0.93:0.036788,stratulus:0.197347)1.00:0.081037)1.00:0.056582)0.77:0.021826,((alutaceus:0.161906,vanidicus:0.205996)1.00:0.118216,((angusticeps:0.085710,paternus:0.059511)1.00:0.153413,loysiana:0.183628)1.00:0.042858)1.00:0.057139,(marcanoi:0.235912,strahmi:0.197766)1.00:0.141032,darwinii:0.636493)1.00:0.067869,(ophiolepis:0.094501,sagrei:0.096758)1.00:0.179398)0.96:0.044895)",
        ]
                
    def check_tree(self, tree_str):          
        d1 = datasets.Dataset()
        tree1 = d1.trees_from_string(tree_str, format="newick")[0]
        pa, edge_lens = to_parent_array(tree1, True, False)
        _LOG.info('Original tree: %s' % tree_str)
        cmd = self.prog_path + " " + " ".join(pa)
        stdout, stderr, returncode = run_program(cmd)
        assert returncode == 0, "Program exited with error:\n%s" % stderr
        _LOG.info('Returned tree: %s' % stdout)
        tree2 = d1.trees_from_string(stdout, format="newick")[0]
        splits.encode_splits(tree1)
        splits.encode_splits(tree2)
        d = treedists.symmetric_difference(tree1, tree2)
        assert d == 0, "Symmetric distance = %d:\n%s;\n%s;" % (d, tree_str, stdout)

    def testRoundTripParentArray(self):
        _LOG.info('Testing round-tripping of trees into GINGKO parent array and back')
        for t in self.src_trees:
            self.check_tree(t)
            
class GenealogyTreeTest(unittest.TestCase):
    def setUp(self):
        self.prog_path = get_gingko_program_path("test_tree_from_genealogies")

    def testTreesFromGenealogies(self):
        run_external_tests(self.prog_path, _LOG, "testing trees from genealogies") 
            
if __name__ == "__main__":
    unittest.main()
    