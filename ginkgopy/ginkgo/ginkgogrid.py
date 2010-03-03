#! /usr/bin/env python

import random
import sys
import os
from ginkgo import argparse

###############################################################################\\
# Basic grids

class Grid(object):

    def __init__(self, ncols, nrows, **kwargs):
        self.ncols = ncols
        self.nrows = nrows
        if 'values' in kwargs:
            self.values = kwargs['values']
        elif 'func' in kwargs:
            self.populate(kwargs['func'])
        else:
            self.values = {}
        self.is_int = kwargs.get('is_int', False)

    def format_value(self, value):
        return str(value)

    def populate(self, func):
        self.values = {}
        for x in range(self.ncols):
            self.values[x] = {}
            for y in range(self.nrows):
                self.values[x][y] = func(x, y)

    def grid_header(self):
        return ("""ncols         {0}
nrows         {1}
xllcorner     0.0
yllcorner     0.0
cellsize      50.0""").format(self.ncols, self.nrows)

    def __str__(self):
        rows = []
        rows.append(self.grid_header())
        max_field_len = 8
        for y in range(self.nrows):
            if y % 5 == 0:
                rows.append("")
            row = []
            for x in range(self.ncols):
                leader = " " if (x and (x % 5 == 0)) else ""
                v = self.values[x][y]
                if isinstance(v, float):
                    row.append("{2}{0:>{1}.4}".format(v, max_field_len, leader))
                elif isinstance(v):
                    row.append("{2}{0:>{1}}".format(v, max_field_len, leader))
            rows.append("  ".join(row))
        return "\n".join(rows)

class RealGrid(Grid):

    def __init__(self, ncols, nrows, **kwargs):
        Grid.__init__(self, ncols, nrows, **kwargs)
        self.is_int = False

class IntGrid(Grid):

    def __init__(self, ncols, nrows, **kwargs):
        Grid.__init__(self, ncols, nrows, **kwargs)
        self.is_int = True

###############################################################################\\
# Build and return the grids

def random_gaussian_grid(ncols, nrows, mean=0, sd=1, output=sys.stdout):
    return RealGrid(ncols, nrows, func = lambda x, y: random.gauss(mean, sd))

def random_uniform_real_grid(ncols, nrows, a, b, output=sys.stdout):
    return RealGrid(ncols, nrows, func = lambda x, y: random.gauss(a, b))

def random_uniform_int_grid(ncols, nrows, a, b, output=sys.stdout):
    return IntGrid(ncols, nrows, func = lambda x, y: random.randint(a, b))

def fixed_value_grid(ncols, nrows, val, output=sys.stdout):
    return IntGrid(ncols, nrows, func = lambda x, y: val)

