#! /usr/bin/env python

###############################################################################
##
## GINGKO Biogeographical Evolution Simulator Post-Processing Library.
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

import re
import math

try:
    from rpy2 import robjects
except:
    pass

from dendropy import treedists

def georeference_tree_taxa(tree):
    """
    Decorates taxa with x and y cell coordinates as parsed from the taxon 
    labels.
    """
    for taxon in tree.taxa_block:
        taxon.x = int(re.match(".* x([\d]+) .*", taxon.label).groups(1)[0])
        taxon.y = int(re.match(".* y([\d]+) .*", taxon.label).groups(1)[0])        
    
def euclidean_distance(p1, p2):
    """
    Returns the euclidean distance between two points, `p1` and `p2`, which
    each is tuple of two integers, (x,y).
    """
    return math.sqrt( pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) )
            
def geo_distance(node1, node2):
    """
    Returns the spatial euclidean distance between two nodes, assuming that the
    node taxa have been marked up with "x" and "y" attributes.
    """
    return euclidean_distance((node1.taxon.x, node1.taxon.y), (node2.taxon.x, node2.taxon.y))

def calc_pearsons_r(tree):
    """
    Calculates Pearson's R on a georeferenced tree for geographic vs. patristic 
    distances. Taxa on trees must be georeferenced, and splits encoded.
    """
    gd = []
    pd = []
    leaves = [x for x in tree.leaf_iter()]
    for idx1, leaf1 in enumerate(leaves):
        for idx2, leaf2 in enumerate(leaves[idx1+1:]):
            gd.append(geo_distance(leaf1, leaf2))
            pd.append(treedists.patristic_distance(leaf1, leaf2))    
    k = robjects.r.cor(robjects.IntVector(gd), robjects.IntVector(pd), method='pearson')
    tree.pearson_r_val = k[0]
    return tree.pearson_r_val
    