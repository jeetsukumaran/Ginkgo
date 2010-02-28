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

class World(object):

    def __init__(self, x_range, y_range, num_gens, **kwargs):
        self.x_range = x_range
        self.y_range = y_range
        self.num_gens = num_gens
        self.label = kwargs.get("label", None)
        self.num_fitness_traits = kwargs.get("num_fitness_traits", None)
        self.global_selection_strength = kwargs.get("global_selection_strength", None)
        self.default_cell_carrying_capacity = kwargs.get("default_cell_carrying_capacity", None)
        self.random_seed = kwargs.get("random_seed", None)
        self.log_frequency = kwargs.get("log_frequency", None)
        self.multifurcating_trees = kwargs.get("multifurcating_trees", None)
        self.final_output = kwargs.get("final_output", None)
        self.full_complement_diploid_trees = kwargs.get("full_complement_diploid_trees", None)
        self.lineages = []
        self.environments = []
        self.initialization_regime = kwargs.get("initialization_regime", None)
        self.samples = []
        self._child_elements = [self.lineages, self.environments, self.samples]

    def __str__(self):
        s = StringIO()
        s.write('<?xml version="1.0"?>\n')
        s.write('<ginkgo>\n')
        indent = (INDENT_SIZE * ' ')
        parts = []
        parts.append('world')
        parts.append('x_range="%s"' % self.x_range)
        parts.append('y_range="%s"' % self.y_range)
        parts.append('num_gens="%s"' % self.num_gens)
        if self.label:
            parts.append('label="%s"' % self.label)
        if self.num_fitness_traits:
            parts.append('num_fitness_traits="%s"' % self.num_fitness_traits)
        if self.global_selection_strength:
            parts.append('global_selection_strength="%s"' % self.num_fitness_traits)
        if self.default_cell_carrying_capacity:
            parts.append('default_cell_carrying_capacity="%s"' % self.default_cell_carrying_capacity)
        if self.random_seed:
            parts.append('random_seed="%s"' % self.num_fitness_traits)
        if self.log_frequency:
            parts.append('log_frequency="%s"' % self.num_fitness_traits)
        if self.multifurcating_trees:
            parts.append('multifurcating_trees="%s"' % self.num_fitness_traits)
        if self.final_output:
            parts.append('final_output="%s"' % self.num_fitness_traits)
        if self.full_complement_diploid_trees:
            parts.append('full_complement_diploid_trees="%s"' % self.num_fitness_traits)
        sep = "\n%s" % (indent * 2)
        world_elem = sep.join(parts)
        s.write('%s<%s>\n' % (indent, world_elem))
        sub_indent = indent * 2
        if self.lineages:
            s.write('%s<lineages>\n' % sub_indent)
            for lineage in self.lineages:
                s.write(str(lineage))
            s.write('%s</lineages>\n' % sub_indent)
        if self.initialization_regime:
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
        s.write('%s</world>\n' % indent)
        s.write('</ginkgo>\n')
        return s.getvalue()

class Lineage(object):

    def __init__(self, lineage_id, **kwargs):
        self.lineage_id = lineage_id
        self.fitness_trait_relative_selection_weights = kwargs.get("fitness_trait_relative_selection_weights", None)
        self.fitness_trait_default_genotypes = kwargs.get("fitness_trait_default_genotypes", None)
        self.fecundity = kwargs.get("fecundity", None)
        self.__movement_capacity = None
        self.movement_capacity = kwargs.get("movement_capacity", None)
        self.indent_level = kwargs.get('indent_level', 3)

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
            s.write('%s<movementCapacity distribution="%s">%s</movementCapacity>\n' \
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

    def __init__(self, gen, indent_level=3):
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
        self.indent_level = kwargs.get("indent_level", 4)

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

    def __init__(self, pops, environment, max_cycles=None, indent_level=2):
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
        self.indent_level = kwargs.get('indent_level', 3)
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
        parts.append('%s<sample gen="%d" lineage=%s%s%s' % (top_indent, self.gen, self.lineage_id, label, close))
        if self.build_trees:
            sub_indent = ((self.indent_level + 1) * INDENT_SIZE) * ' '
            if self.individuals_per_cell:
                parts.append('%s<individualsPerCell>%s</individualsPerCell>' % (sub_indent, self.individuals_per_cell))
            if self.cell_indexes:
                parts.append('%s<cellIndexes>' % sub_indent)
                sub_sub_indent = ((self.indent_level + 2) * INDENT_SIZE) * ' '
                if isinstance(self.cell_indexes, str):
                    rows = self.cell_indexes.split('\n')
                    for r in rows:
                        parts.append('%s%s' % (sub_sub_indent, r))
                else:
                    parts.append('%s%s' % (sub_sub_indent, " ".join([str(s) for s in self.cell_indexes])))
                parts.append('%s</cellIndexes>' % sub_indent)
            elif self.cell_coordinates:
                parts.append('%s<cellCoordinates>' % sub_indent)
                sub_sub_indent = ((self.indent_level + 2) * INDENT_SIZE) * ' '
                if isinstance(self.cell_coordinates, str):
                    rows = self.cell_coordinates.split('\n')
                    for r in rows:
                        parts.append('%s%s' % (sub_sub_indent, r))
                else:
                    for c in self.cell_coordinates:
                        parts.append('%s%s,%s' % (sub_sub_indent, c[0], c[1]))
                parts.append('%s</cellCoordinates>' % sub_indent)
            parts.append('%s</sample>\n' % top_indent)
        return "\n".join(parts)

if __name__ == "__main__":
    num_fitness_traits = 1
    w = World(
            x_range=10,
            y_range=10,
            num_gens=20000,
            label="ginkgo_run",
            num_fitness_traits=num_fitness_traits,
            global_selection_strength=1.0,
            default_cell_carrying_capacity=100,
            log_frequency=1000,
            multifurcating_trees=True,
            final_output=False,
            full_complement_diploid_trees=False)
    cells = [(1,1),(2,1),(3,1),(4,1),
             (1,2),(2,2),(3,2),(4,2)]
    s1 = Lineage("Zx",
            fitness_trait_relative_selection_weights=[1]*num_fitness_traits,
            fitness_trait_default_genotypes=[0]*num_fitness_traits,
            fecundity=16,
            movement_capacity=("poisson", 10))
    w.lineages.append(s1)

    e0 = Environment(None)
    e0.grids.append(CarryingCapacityGrid('cc.grd'))
    e0.grids.append(MovementProbabilitiesGrid('Zx', 'movp_zx.grd'))
    e0.grids.append(MovementCostsGrid('Zx', 'movp_zx.grd'))
    e0.grids.append(FitnessTraitOptimaGrid('1', 'movp_zx.grd'))
    cell_pops = []
    for cell in cells:
        cell_pops.append(InitializationCellPopulations(x=cell[0], y=cell[1], pops={"Zx": 100}))
    init_reg = InitializationRegime(pops=cell_pops, environment=e0)
    w.initialization_regime = init_reg

    e1 = Environment(1)
    e1.grids.append(CarryingCapacityGrid('cc.grd'))
    e1.grids.append(MovementProbabilitiesGrid('Zx', 'movp_zx.grd'))
    e1.grids.append(MovementCostsGrid('Zx', 'movp_zx.grd'))
    e1.grids.append(FitnessTraitOptimaGrid('1', 'movp_zx.grd'))
    w.environments.append(e1)

    for t in xrange(10):
        gen = (t + 10) * 1000
        s = Sample(lineage_id="Zx", gen=gen, individuals_per_cell=10, cells=cells)
        w.samples.append(s)

    print(str(w))
