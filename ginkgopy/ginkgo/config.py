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

class CarryingCapacityGrid(EnvironmentGrid):

    def __init__(self, grid_filepath):
        EnvironmentGrid.__init__(self,
                "carrying_capacity",
                attribute_dict={},
                grid_filepath=grid_filepath)

class MovementProbabilitiesGrid(EnvironmentGrid):

    def __init__(self, lineage, grid_filepath):
        EnvironmentGrid.__init__(self,
                "movementProbabilities",
                attribute_dict={'lineage': lineage},
                grid_filepath=grid_filepath)

class MovementCostsGrid(EnvironmentGrid):

    def __init__(self, lineage, grid_filepath):
        EnvironmentGrid.__init__(self,
                "movementCosts",
                attribute_dict={'lineage': lineage},
                grid_filepath=grid_filepath)

class FitnessTraitOptimaGrid(EnvironmentGrid):

    def __init__(self, trait, grid_filepath):
        EnvironmentGrid.__init__(self,
                "fitnessTraitOptima",
                attribute_dict={'trait': trait},
                grid_filepath=grid_filepath)

class Environment(object):

    def __init__(self, gen, indent_level=0):
        self.gen = gen
        self.grids = []
        self.indent_level = indent_level

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
        self.individuals_per_cell = kwargs.get('individuals_per_cell', None)
        self.cells = kwargs.get('cells', None)
        self.indent_level = kwargs.get('indent_level', 0)

    def __str__(self):
        top_indent = (self.indent_level * INDENT_SIZE) * ' '
        if self.label:
            label = ' label="%s"' % label
        else:
            label = ''
        parts.append('%s<sample gen="%d" lineage=%s%s>' % (top_indent, self.gen, self.lineage, label))
        sub_indent = ((self.indent_level + 1) * INDENT_SIZE) * ' '
        if self.individuals_per_cells:
            parts.append('%s<individualsPerCell>%s</individualsPerCell>' % (sub_indent, self.individuals_per_cells))
        if self.cells:
            parts.append('%s<cells>' % sub_indent)
            sub_sub_indent = ((self.indent_level + 2) * INDENT_SIZE) * ' '
            for c in self.cells:
                parts.append('%s%s,%s' % (sub_sub_indent, c[0], c[1]))
            parts.append('%s</cells>' % sub_indent)
        parts.append('%s</sample>' % top_indent)
        return "\n".join(parts)

e = Environment(0)
e.grids.append(CarryingCapacityGrid('cc.grd'))
e.grids.append(MovementProbabilitiesGrid('Zx', 'movp_zx.grd'))
e.grids.append(MovementCostsGrid('Zx', 'movp_zx.grd'))
e.grids.append(FitnessTraitOptimaGrid('1', 'movp_zx.grd'))
print(e)
