#! /usr/bin/env python

import random
import sys
import os
import baker

_prog_usage = '%prog [options] NCOLS NROWS VAL1 VAL2'
_prog_version = 'GINKGOGRID Version 1.0'
_prog_description = 'Generates grid files used as input for the Ginkgo simulator.'
_prog_author = 'Jeet Sukumaran'
_prog_copyright = 'Copyright (C) 2010 Jeet Sukumaran.'

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
                if self.is_int:
                    row.append("{2}{0:>{1}}".format(v, max_field_len, leader))
                else:
                    row.append("{2}{0:>{1}.4}".format(v, max_field_len, leader))
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

def random_gaussian_grid(ncols, nrows, mean=0, sd=1, output=sys.stdout):
    return RealGrid(ncols, nrows, func = lambda x, y: random.gauss(mean, sd))

def random_uniform_real_grid(ncols, nrows, a, b, output=sys.stdout):
    return RealGrid(ncols, nrows, func = lambda x, y: random.gauss(a, b))

def random_uniform_int_grid(ncols, nrows, a, b, output=sys.stdout):
    return IntGrid(ncols, nrows, func = lambda x, y: random.randint(a, b))

def fixed_value_grid(ncols, nrows, val, output=sys.stdout):
    return IntGrid(ncols, nrows, func = lambda x, y: val)

if __name__ == '__main__':
    print random_uniform_int_grid(4, 5, 1, 10)
