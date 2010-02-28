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

class GinkgoConfiguration(object):

    def __init__(self, **kwargs):
        self.title = kwargs.get("title", "results")
        self.log_frequency = kwargs.get("log_frequency", None)
        self.multifurcating_trees = kwargs.get("multifurcating_trees", None)
        self.final_output = kwargs.get("final_output", None)
        self.full_complement_diploid_trees = kwargs.get("full_complement_diploid_trees", None)
        self.system = kwargs.get("system", None)
        self.landscape = kwargs.get("landscape", None)
        self.lineages = kwargs.get("lineages", [])
        self.initialization_regime = kwargs.get("initialization_regime", None)
        self.environments = kwargs.get("environments", [])
        self.samples = kwargs.get("samples", [])

    def __str__(self):
        s = StringIO()
        s.write('<?xml version="1.0"?>\n')
        s.write('<ginkgo')
        indent = INDENT_SIZE * ' '
        if self.title is not None:
            s.write('\n%stitle="%s"' % (indent, self.title))
        if self.log_frequency is not None:
            s.write('\n%slog_frequency="%s"' % (indent, self.log_frequency))
        if self.multifurcating_trees is not None:
            s.write('\n%smultifurcating_trees="%s"' % (indent, self.multifurcating_trees))
        if self.final_output is not None:
            s.write('\n%sfinal_output="%s"' % (indent, self.final_output))
        if self.full_complement_diploid_trees is not None:
            s.write('\n%sfull_complement_diploid_trees="%s"' % (indent, self.full_complement_diploid_trees))
        s.write('>\n')
        sub_indent = indent * 1
        if self.system is not None:
            s.write(str(self.system))
        if self.landscape is not None:
            s.write(str(self.landscape))
        if self.lineages:
            s.write('%s<lineages>\n' % sub_indent)
            for lineage in self.lineages:
                s.write(str(lineage))
            s.write('%s</lineages>\n' % sub_indent)
        if self.initialization_regime is not None:
            s.write(str(self.initialization_regime))
        if self.environments:
            s.write('%s<environments>\n' % sub_indent)
            for environment in self.environments:
                s.write(str(environment))
            s.write('%s</environments>\n' % sub_indent)
        if self.samples:
            s.write('%s<samples>\n' % sub_indent)
            for sample in self.samples:
                s.write(str(sample))
            s.write('%s</samples>\n' % sub_indent)
        s.write('</ginkgo>\n')
        return s.getvalue()

class System(object):

    def __init__(self, **kwargs):
        self.random_seed = kwargs.get("random_seed", None)
        self.fitness_dimensions = kwargs.get("fitness_dimensions", None)
        self.global_selection_strength = kwargs.get("global_selection_strength", None)
        self.ngens = kwargs.get("ngens", None)

    def __str__(self):
        s = StringIO()
        indent = (INDENT_SIZE * ' ')
        sub_indent = indent * 2
        s.write('%s<system>' % indent)
        if self.random_seed is not None:
            s.write('\n%s<random_seed>%s</random_seed>' % (sub_indent, self.random_seed))
        if self.fitness_dimensions is not None:
            s.write('\n%s<fitness_dimensions>%s</fitness_dimensions>' % (sub_indent, self.fitness_dimensions))
        if self.global_selection_strength is not None:
            s.write('\n%s<global_selection_strength>%s</global_selection_strength>' % (sub_indent, self.global_selection_strength))
        if self.ngens is not None:
            s.write('\n%s<ngens>%s</ngens>' % (sub_indent, self.ngens))
        s.write('\n%s</system>\n' % indent)
        return s.getvalue()

class Landscape(object):

    def __init__(self, ncols, nrows, **kwargs):
        self.ncols = ncols
        self.nrows = nrows
        self.default_cell_carrying_capacity = kwargs.get("default_cell_carrying_capacity", None)

    def __str__(self):
        s = StringIO()
        indent = (INDENT_SIZE * ' ')
        sub_indent = indent * 2
        s.write('%s<landscape ncols="%s" nrows="%s">' % (indent, self.ncols, self.nrows))
        if self.default_cell_carrying_capacity is not None:
            s.write('\n%s<default_cell_carrying_capacity>%s</default_cell_carrying_capacity>' % (sub_indent, self.default_cell_carrying_capacity))
        s.write('\n%s</landscape>\n' % indent)
        return s.getvalue()

class Lineage(object):

    def __init__(self, lineage_id, **kwargs):
        self.lineage_id = lineage_id
        self.fitness_trait_relative_selection_weights = kwargs.get("fitness_trait_relative_selection_weights", None)
        self.fitness_trait_default_genotypes = kwargs.get("fitness_trait_default_genotypes", None)
        self.fecundity = kwargs.get("fecundity", None)
        self.__movement_capacity = None
        self.movement_capacity = kwargs.get("movement_capacity", None)
        self.indent_level = kwargs.get('indent_level', 2)

    def _set_movement_capacity(self, val):
        if not (isinstance(val, list) or isinstance(val, tuple))\
                or len(val) != 2 \
                or not isinstance(val[0], str):
            raise ValueError("Expecting tuple of (<dist>, <val>) for movement capacity, " \
                           + "where <dist> is string specifying the distribution ('constant' or 'poisson'), " \
                           + "and <val> is either the fixed movement capacity value or the mean of the Poisson")
        self.__movement_capacity = val

    def _get_movement_capacity(self):
        return self.__movement_capacity

    movement_capacity = property(_get_movement_capacity, _set_movement_capacity)

    def __str__(self):
        top_indent = (self.indent_level * INDENT_SIZE) * ' '
        s = StringIO()
        s.write('%s<lineage id="%s">\n' % (top_indent, self.lineage_id))
        sub_indent = ((self.indent_level + 1) * INDENT_SIZE) * ' '
        if self.fitness_trait_relative_selection_weights:
            s.write('%s<fitness_trait_relative_selection_weights>%s</fitness_trait_relative_selection_weights>\n' %
                    (sub_indent, (" ".join([str(x) for x in self.fitness_trait_relative_selection_weights]))))
        if self.fitness_trait_default_genotypes:
            s.write('%s<fitness_trait_default_genotypes>%s</fitness_trait_default_genotypes>\n' %
                    (sub_indent, (" ".join([str(x) for x in self.fitness_trait_default_genotypes]))))
        if self.fecundity:
            s.write('%s<fecundity>%s</fecundity>\n' % (sub_indent, self.fecundity))
        if self.movement_capacity:
            s.write('%s<movement_capacity distribution="%s">%s</movement_capacity>\n' \
                    % (sub_indent, self.movement_capacity[0], self.movement_capacity[1]))
        s.write('%s</lineage>\n' % (top_indent))
        return s.getvalue()

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

#class MovementProbabilitiesGrid(EnvironmentGrid):
#
#    def __init__(self, lineage, grid_filepath):
#        EnvironmentGrid.__init__(self,
#                "movement_probabilities",
#                attribute_dict={'lineage': lineage},
#                grid_filepath=grid_filepath)

class MovementCostsGrid(EnvironmentGrid):

    def __init__(self, lineage, grid_filepath):
        EnvironmentGrid.__init__(self,
                "movement_costs",
                attribute_dict={'lineage': lineage},
                grid_filepath=grid_filepath)

class FitnessTraitOptimaGrid(EnvironmentGrid):

    def __init__(self, trait, grid_filepath):
        EnvironmentGrid.__init__(self,
                "fitness_trait_optima",
                attribute_dict={'trait': trait},
                grid_filepath=grid_filepath)

class Environment(object):

    def __init__(self, gen, indent_level=2):
        self.gen = gen
        self.grids = []
        self.indent_level = indent_level

    def __str__(self):
        if len(self.grids) == 0:
            return ""
        parts = []
        top_indent = (self.indent_level * INDENT_SIZE) * ' '
        if self.gen is not None:
            gen = ' gen="%d"' % self.gen
        else:
            gen = ""
        parts.append('%s<environment%s>' % (top_indent, gen))
        sub_indent = ((self.indent_level + 1) * INDENT_SIZE) * ' '
        for grid in self.grids:
            parts.append('%s%s' % (sub_indent, str(grid)))
        parts.append('%s</environment>\n' % top_indent)
        return "\n".join(parts)

class InitializationCellPopulations(object):

    def __init__(self, **kwargs):
        self.x = kwargs.get("x", None)
        self.y = kwargs.get("y", None)
        self.index = kwargs.get("index", None)
        self.pops = kwargs.get("pops", {})
        self.indent_level = kwargs.get("indent_level", 3)

    def __str__(self):
        parts = []
        top_indent = (self.indent_level * INDENT_SIZE) * ' '
        if self.index is not None:
            parts.append('%s<cell index="%s">' % (top_indent, self.index))
        elif self.x is not None and self.y is not None:
            parts.append('%s<cell x="%s" y="%s">' % (top_indent, self.x, self.y))
        else:
            raise TypeError("Either 'index' or both 'x' and 'y' must be specified")
        sub_indent = ((self.indent_level+1) * INDENT_SIZE) * ' '
        for sp, sz in self.pops.items():
            parts.append('%s<population lineage="%s" size="%s" />' % (sub_indent, sp, sz))
        parts.append("%s</cell>" % top_indent)
        return "\n".join(parts)

class InitializationRegime(object):

    def __init__(self, pops, environment, max_cycles=None, indent_level=1):
        self.pops = pops
        self.environment = environment
        self.max_cycles = max_cycles
        self.indent_level = indent_level

    def __str__(self):
        parts = []
        top_indent = (self.indent_level * INDENT_SIZE) * ' '
        sub_indent = ((self.indent_level+1) * INDENT_SIZE) * ' '
        if self.max_cycles is not None:
            max_cycles = ' max_cycles="%s"' % self.max_cycles
        else:
            max_cycles = ""
        parts.append('%s<initialization%s>' % (top_indent, max_cycles))
        e = str(self.environment)
        if e.endswith('\n'):
            e = e[:-1]
        parts.append(e)
        parts.append('%s<populations>' % (sub_indent))
        for cp in self.pops:
            parts.append(str(cp))
        parts.append('%s</populations>' % (sub_indent))
        parts.append("%s</initialization>\n" % top_indent)
        return "\n".join(parts)

class Sample(object):

    def __init__(self, lineage_id, gen, **kwargs):
        self.lineage_id = lineage_id
        self.gen = gen
        self.label = kwargs.get('label', None)
        self.individuals_per_cell = kwargs.get('individuals_per_cell', None)
        self.cell_coordinates = kwargs.get('cell_coordinates', None)
        self.cell_indexes = kwargs.get('cell_indexes', None)
        self.indent_level = kwargs.get('indent_level', 2)
        self.build_trees = kwargs.get('build_trees', True)

    def __str__(self):
        top_indent = (self.indent_level * INDENT_SIZE) * ' '
        if self.label:
            label = ' label="%s"' % self.label
        else:
            label = ''
        parts = []
        if self.build_trees:
            close = '>'
        else:
            close = ' trees="False" />\n'
        parts.append('%s<sample gen="%d" lineage="%s"%s%s' % (top_indent, self.gen, self.lineage_id, label, close))
        if self.build_trees:
            sub_indent = ((self.indent_level + 1) * INDENT_SIZE) * ' '
            if self.individuals_per_cell:
                parts.append('%s<individuals_per_cell>%s</individuals_per_cell>' % (sub_indent, self.individuals_per_cell))
            if self.cell_indexes or self.cell_coordinates:
                parts.append('%s<cells>' % sub_indent)
                sub_sub_indent = ((self.indent_level + 2) * INDENT_SIZE) * ' '
                if self.cell_indexes:
                    for c in self.indexes:
                        parts.append('%s<cell index="%s" />' % (sub_sub_indent, c))
                if self.cell_coordinates:
                    for c in self.cell_coordinates:
                        parts.append('%s<cell x="%s" y="%s"/>' % (sub_sub_indent, c[0], c[1]))
                parts.append('%s<cells>' % sub_indent)
            parts.append('%s</sample>\n' % top_indent)
        return "\n".join(parts)

if __name__ == "__main__":
    num_fitness_traits = 1

    g = GinkgoConfiguration(
            title="results1",
            log_frequency=100,
            multifurcating_trees=True,
            final_output=False,
            full_complement_diploid_trees=False,
            system=System(random_seed=3124,
                         fitness_dimensions=num_fitness_traits,
                         global_selection_strength=1.0,
                         ngens=1000),
            landscape=Landscape(ncols=10, nrows=10, default_cell_carrying_capacity=100)
        )
    cells = [(1,1),(2,1),(3,1),(4,1),
             (1,2),(2,2),(3,2),(4,2)]
    s1 = Lineage("Zx",
            fitness_trait_relative_selection_weights=[1]*num_fitness_traits,
            fitness_trait_default_genotypes=[0]*num_fitness_traits,
            fecundity=16,
            movement_capacity=("poisson", 10))
    g.lineages.append(s1)

    e0 = Environment(None)
    e0.grids.append(CarryingCapacityGrid('cc.grd'))
    e0.grids.append(MovementCostsGrid('Zx', 'movp_zx.grd'))
    e0.grids.append(FitnessTraitOptimaGrid('1', 'movp_zx.grd'))
    cell_pops = []
    for cell in cells:
        cell_pops.append(InitializationCellPopulations(x=cell[0], y=cell[1], pops={"Zx": 100}))
    init_reg = InitializationRegime(pops=cell_pops, environment=e0)
    g.initialization_regime = init_reg

    e1 = Environment(1)
    e1.grids.append(CarryingCapacityGrid('cc.grd'))
    e1.grids.append(MovementCostsGrid('Zx', 'movp_zx.grd'))
    e1.grids.append(FitnessTraitOptimaGrid('1', 'movp_zx.grd'))
    g.environments.append(e1)

    for t in xrange(10):
        gen = (t + 10) * 1000
        s = Sample(lineage_id="Zx", gen=gen, individuals_per_cell=10, cell_coordinates=cells)
        g.samples.append(s)

    print(str(g))
