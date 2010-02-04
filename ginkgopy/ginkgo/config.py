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
Helps in building configuration files.
"""

from cStringIO import StringIO

INDENT_SIZE = 4

class EnvironmentGrid(object):

    def __init__(self, type_name, attribute_dict, grid_filepath):
        self.type_name = type_name
        self.attribute_dict = attribute_dict
        self.grid_filepath = grid_filepath

    def __str__(self):
        s = StringIO()
        s.write("<%s" % self.type_name)
        if self.attribute_dict:
            parts = []
            for k, v in self.attribute_dict.items():
                parts.append('%s="%s"' % (k,v))
            s.write(" %s" % (" ".join(parts)))
        s.write(">%s</%s>" % (self.grid_filepath, self.type_name))
        return s.getvalue()

class Environment(object):

    def __init__(self, gen, indent_level=0):
        self.gen = gen
        self.grids = []
        self.indent_level = indent_level

    def add_grid(self, type_name, attribute_dict, grid_filepath):
        self.grids.append(EnvironmentGrid(type_name, attribute_dict, grid_filepath))
        return self.grids[-1]

    def set_carrying_capacity_grid(self, grid_filepath):
        return self.add_grid("carrying_capacity",
                attribute_dict={},
                grid_filepath=grid_filepath)

    def set_movement_probabilities_grid(self, lineage, grid_filepath):
        return self.add_grid("movementProbabilities",
                attribute_dict={'lineage': lineage},
                grid_filepath=grid_filepath)

    def set_movement_costs_grid(self, lineage, grid_filepath):
        return self.add_grid("movementCosts",
                attribute_dict={'lineage': lineage},
                grid_filepath=grid_filepath)

    def set_fitness_trait_optima_grid(self, trait, grid_filepath):
        return self.add_grid("fitnessTraitOptima",
                attribute_dict={"trait": trait},
                grid_filepath=grid_filepath)

    def __str__(self):
        if len(self.grids) == 0:
            return ""
        parts = []
        top_indent = (self.indent_level * INDENT_SIZE) * ' '
        parts.append('%s<environment gen="%d">' % (top_indent, self.gen))
        sub_indent = ((self.indent_level + 1) * INDENT_SIZE) * ' '
        for grid in self.grids:
            parts.append('%s%s' % (sub_indent, str(grid)))
        parts.append('%s</environment>' % top_indent)
        return "\n".join(parts)

def Sample(object):

    def __init__(self, lineage, gen, **kwargs):

        self.lineage = lineage
        self.gen = gen
        self.label = kwargs.get('label', None)
        self.

e = Environment(0)
e.set_carrying_capacity_grid('cc.grd')
e.set_movement_probabilities_grid('Zx', 'movp_zx.grd')
e.set_movement_costs_grid('Zx', 'movp_zx.grd')
e.set_fitness_trait_optima_grid('1', 'movp_zx.grd')
print(e)
