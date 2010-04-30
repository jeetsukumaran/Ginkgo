#! /usr/bin/env python

import random
import sys
import os
from ginkgo import argparse

##############################################################################\\
# Grid

class Grid(object):

    def __init__(self, **kwargs):
        self.ncols = kwargs.get("ncols", None)
        self.nrows = kwargs.get("nrows", None)
        self.value_type = kwargs.get("value_type", int)
        self.values = None
        self.matrix = None
        if 'values' in kwargs:
            self.values = kwargs['values']
        elif 'pop_func' in kwargs:
            self.populate(kwargs['pop_func'])
        elif 'filepath' in kwargs:
            self.read(open(kwargs['filepath'], "rU"))
        elif 'stream' in kwargs:
            self.read(kwargs['stream'])
        else:
            self.values = {}
        self._max_formatted_value_len = None

    def __str__(self):
        return self.as_string(include_header=True)

    def populate(self, func):
        self.values = {}
        for x in range(self.ncols):
            self.values[x] = {}
            for y in range(self.nrows):
                self.values[x][y] = func(x, y)

    def read(self, src):
        self.values = []
        for line in src:
            line = line.replace('\n', '').strip()
            parts = line.split(' ',1)
            kw = parts[0].lower()
            if kw == 'ncols':
                assert len(parts) == 2
                self.ncols = int(parts[1])
                continue
            elif kw == 'nrows':
                assert len(parts) == 2
                self.nrows = int(parts[1])
                continue
            elif kw in ['xllcorner', 'yllcorner', 'cellsize', 'nodata_value']:
                continue
            else:
                parts = line.split(' ')
                self.values.extend([self.value_type(i) for i in parts])
                break
        assert self.ncols > 0
        assert self.nrows > 0

        for line in src:
            line = line.replace('\n', '').strip()
            parts = line.split(' ')
            self.values.extend([self.value_type(i) for i in parts])

        return self.matrix_from_values()

    def matrix_from_values(self):
        assert len(self.values) == self.ncols * self.nrows
        self.matrix = []
        for r in range(self.nrows):
            self.matrix.append([])
            for c in range(self.ncols):
                self.matrix[r].append(self.values[(r * self.ncols) + c])
            assert len(self.matrix[r]) == self.ncols
        assert len(self.matrix) == self.nrows
        return self.matrix

    def formatted_value_matrix(self, cell_width=None):
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
        if cell_width is None:
            self._max_formatted_value_len = max(fv_lens)
        else:
            self._max_formatted_value_len = cell_width
        return fv

    def ascii_grid_header(self):
        return ("""ncols         {0}
nrows         {1}
xllcorner     0.0
yllcorner     0.0
cellsize      50.0
NODATA_value  -9999""").format(self.ncols, self.nrows)

    def as_string(self, include_header=True, cell_width=None):
        rows = []
        if include_header:
            rows.append(self.grid_header())
        fv = self.formatted_value_matrix(cell_width=cell_width)
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
            #rows.append("")
        return "\n".join(rows)

###############################################################################\\
# Occurrences

class Occurrences(Grid):

    def __init__(self, filepath=None):
        Grid.__init__(self)
        self.filepath = None
        if filepath is not None:
            self.read(open(filepath, "rU"))

    def __str__(self):
        s = []
        for r in range(self.nrows):
            s.append(" ".join(["{0:>3}".format(self.matrix[r][c]) for c in range(self.ncols)]))
        return "\n".join(s)

###############################################################################\\
# Input Grid Generation

def random_gaussian_grid(ncols, nrows, mean=0, sd=1):
    return Grid(ncols=ncols, nrows=nrows, pop_func=lambda x, y: random.gauss(mean, sd))

def random_uniform_real_grid(ncols, nrows, a, b):
    return Grid(ncols=ncols, nrows=nrows, pop_func=lambda x, y: random.uniform(a, b))

def random_uniform_int_grid(ncols, nrows, a, b):
    return Grid(ncols=ncols, nrows=nrows, pop_func=lambda x, y: random.randint(a, b))

def fixed_value_grid(ncols, nrows, val):
    return Grid(ncols=ncols, nrows=nrows, pop_func=lambda x, y: val)

