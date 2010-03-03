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
        self._max_formatted_value_len = None

    def formatted_values(self):
        fv = {}
        fv_lens = []
        for x in range(self.ncols):
            fv[x] = {}
            for y in range(self.nrows):
                v = self.values[x][y]
                if isinstance(v, float):
                    fv[x][y] = "{0:>.4}".format(v)
                else:
                    fv[x][y] = "{0:>}".format(v)
                fv_lens.append(len(fv[x][y]))
        self._max_formatted_value_len = max(fv_lens)
        return fv

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
        fv = self.formatted_values()
        for y in range(self.nrows):
            if y % 5 == 0:
                rows.append("")
            row = []
            for x in range(self.ncols):
#                leader = ("{0:{1}}".format(" ", self._max_formatted_value_len)) if (x and (x % 5 == 0)) else ""
                leader = "   " if (x and (x % 5 == 0)) else ""
                v = fv[x][y]
                row.append("{2}{0:>{1}}".format(v, self._max_formatted_value_len, leader))
            rows.append("  ".join(row))
        return "\n".join(rows)

###############################################################################\\
# Build and return the grids

def random_gaussian_grid(ncols, nrows, mean=0, sd=1):
    return Grid(ncols, nrows, func = lambda x, y: random.gauss(mean, sd))

def random_uniform_real_grid(ncols, nrows, a, b):
    return Grid(ncols, nrows, func = lambda x, y: random.uniform(a, b))

def random_uniform_int_grid(ncols, nrows, a, b):
    return Grid(ncols, nrows, func = lambda x, y: random.randint(a, b))

def fixed_value_grid(ncols, nrows, val):
    return Grid(ncols, nrows, func = lambda x, y: val)

