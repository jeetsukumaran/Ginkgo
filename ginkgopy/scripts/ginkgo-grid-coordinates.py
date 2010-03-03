#! /usr/bin/env python

import random
import sys
import os
import math
from ginkgo import ginkgogrid
from ginkgo import argparse


_prog_usage = '%prog [options] NCOLS NROWS'
_prog_version = 'GINKGO-GRID-COORDINATES Version 1.0'
_prog_description = 'Generates grid with display of cell coordinates (or indexes).'
_prog_author = 'Jeet Sukumaran'
_prog_copyright = 'Copyright (C) 2010 Jeet Sukumaran.'

def grid_coords(ncols=10, nrows=10, origin_upper_left=True, as_indexes=False):
    rows = []
    if as_indexes:
        max_field_len = int(math.log(ncols*nrows, 10)) + 1
    else:
        max_field_len = int(math.log(ncols, 10)) + int(math.log(nrows, 10)) + 3
    idx = 0
    for y in xrange(nrows):
        if y % 5 == 0:
            rows.append("")
        row = []
        for x in xrange(ncols):
            leader = " " if (x and (x % 5 == 0)) else ""
            if as_indexes:
                row.append("{2}{0:>{1}}".format("{0}".format(idx), max_field_len, leader))
            else:
                if origin_upper_left:
                    yc = y
                else:
                    yc = nrows - y - 1
                row.append("{2}{0:>{1}}".format("{0},{1}".format(x, yc), max_field_len, leader))
            idx += 1
        rows.append("  ".join(row))
    rows.append("")
    return "\n".join(rows)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=_prog_description)
    parser.add_argument('-o', '--output-file', default=None)
    parser.add_argument('-O', '--origin',
                        default='ul',
                        choices=['ul', 'll'],
                        help="grid origin ('ul': upper-left [default], 'll': lower-left)")
    parser.add_argument('--ul',
                        action='store_const',
                        const='ul',
                        dest='origin',
                        help="use upper left-hand corner as grid origin [default]")
    parser.add_argument('--ll',
                        action='store_const',
                        const='ll',
                        dest='origin',
                        help="use lower left-hand corner as grid origin")
    parser.add_argument('-i', '--indexes',
                        action='store_true',
                        default=False,
                        dest='show_indexes',
                        help='show grid cell indexes instead of coordinates')
    parser.add_argument('ncols',
                        type=int,
                        help="number of columns in the grid (i.e., the x-dimension)")
    parser.add_argument('nrows',
                        type=int,
                        help="number of rows in the grid (i.e., the y-dimension)")
    args = parser.parse_args()
    if args.output_file is None:
        out = sys.stdout
    else:
        out = open(os.expanduser(os.expandpaths(args.output_file)))
    out.write(grid_coords(args.ncols, args.nrows, args.origin=='ul', as_indexes=args.show_indexes))




