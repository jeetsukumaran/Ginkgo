#! /usr/bin/env python

###############################################################################
##
## GINKGO Biogeographical Evolution Simulator Post-Processing Library.
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

"""
Helps in performing calculations/analysis on/of results.
"""

import re
import math

try:
    from rpy2 import robjects
except:
    pass

import dendropy
from dendropy import treecalc

class GeoreferencingError(Exception):
    def __init__(self, err_type, text):
        Exception.__init__(self, "%s: '%s'" % (err_type, text))

class GeoreferencingParseError(GeoreferencingError):
    def __init__(self, text):
        GeoreferencingError.__init__(self, "Failed to parse coordinates", text)

class GeoreferencingValueError(GeoreferencingError):
    def __init__(self, text):
        GeoreferencingError.__init__(self, "Invalid coordinate value(s)", text)

class GeoreferencingIndexError(GeoreferencingError):
    def __init__(self, text):
        GeoreferencingError.__init__(self, "Insufficient or incomplete coordinates", text)

_NODE_XY_PAT = re.compile(r"x([\d]+) y([\d]+)")
_TAXON_XY_PAT = re.compile(r".* x([\d]+) y([\d]+) .*")
_TAXON_CELL_INDEX_PAT = re.compile(".* i([\d]+) .*")

def _process_xy_pattern(pattern, text):
    xy_match = pattern.match(text)
    if xy_match is None:
        raise GeoreferencingParseError(text)
    try:
        return int(xy_match.group(1)), int(xy_match.group(2))
    except ValueError:
        raise GeoreferencingValueError(text)
    except IndexError:
        raise GeoreferencingIndexError(text)

def georeference_taxa(taxon_set):
    """
    Decorates taxa with x and y cell coordinates as parsed from the taxon
    labels.
    """
    for taxon in taxon_set:
        taxon.cell_index = int(_TAXON_CELL_INDEX_PAT.match(taxon.label).group(1))
        taxon.x, taxon.y = _process_xy_pattern(_TAXON_XY_PAT, taxon.label)

def georeference_nodes(tree, ignore_errors=False):
    """
    Decorates nodes with x and y cell coordinates, as parsed from the node
    labels.
    """
    for nd in tree.postorder_node_iter():
        if nd.label is not None:
            label = nd.label
            xypat = _NODE_XY_PAT
        elif nd.taxon is not None:
            label = nd.taxon.label
            xypat = _TAXON_XY_PAT
        else:
            nd.x = None
            nd.y = None
            continue
        try:
            nd.x, nd.y = _process_xy_pattern(xypat, label)
        except GeoreferencingError:
            if not ignore_errors:
                raise
            else:
                nd.x = None
                nd.y = None

class GeographicDistanceMatrix(object):
    """
    Calculates and maintains geographic distance information of taxa on a tree.
    """

    def __init__(self, taxon_set=None):
        self.taxon_set = None
        self._geo_dists = {}
        if taxon_set is not None:
            self.calc(taxon_set)

    def __call__(self, taxon1, taxon2):
        """
        Returns patristic distance between two taxon objects.
        """
        try:
            return self._geo_dists[taxon1][taxon2]
        except KeyError, e:
            return self._geo_dists[taxon2][taxon1]

    def calc(self, taxon_set):
        """
        Calculates the distances.
        """
        self.taxon_set = taxon_set
        georeference_taxa(self.taxon_set)
        self._geo_dists = {}
        for i1, t1 in enumerate(self.taxon_set):
            self._geo_dists[t1] = {}
            for i2, t2 in enumerate(self.taxon_set[i1+1:]):
                self._geo_dists[t1][t2] = euclidean_distance( (t1.x, t1.y), (t2.x, t2.y) )

    def distances(self):
        """
        Returns list of patristic distances.
        """
        dists = []
        for dt in self._geo_dists.values():
            for d in dt.values():
                dists.append(d)
        return dists

    def sum_of_distances(self):
        """
        Returns sum of patristic distances on tree.
        """
        return sum(self.distances())

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
    pdm = PatristicDistanceMatrix(tree)
    leaves = [x for x in tree.leaf_iter()]
    for idx1, leaf1 in enumerate(leaves):
        for idx2, leaf2 in enumerate(leaves[idx1+1:]):
            gd.append(geo_distance(leaf1, leaf2))
            pd.append(pdm(leaf1.taxon, leaf2.taxon))
    k = robjects.r.cor(robjects.IntVector(gd), robjects.IntVector(pd), method='pearson')
    tree.pearson_r_val = k[0]
    return tree.pearson_r_val

