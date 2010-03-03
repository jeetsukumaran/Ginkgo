#! /usr/bin/env python

import random
import sys
import os
from ginkgo import ginkgogrid
from ginkgo import argparse


_prog_usage = '%prog [options] NCOLS NROWS VAL1 VAL2'
_prog_version = 'GINKGOGRID Version 1.0'
_prog_description = 'Generates ESRI ASCII format grid files (for use in the Ginkgo simulator).'
_prog_author = 'Jeet Sukumaran'
_prog_copyright = 'Copyright (C) 2010 Jeet Sukumaran.'

def gen_rand_gauss(args):
    g = ginkgogrid.random_gaussian_grid(args.ncols, args.nrows, args.mean, args.sd)
    args.output.write(str(g))

def gen_rand_unif_real(args):
    g = ginkgogrid.random_uniform_real_grid(args.ncols, args.nrows, args.a, args.b)
    args.output.write(str(g))

def gen_rand_unif_int(args):
    g = ginkgogrid.random_uniform_int_grid(args.ncols, args.nrows, args.a, args.b)
    args.output.write(str(g))

def gen_fixed_val(args):
    g = ginkgogrid.fixed_value_grid(args.ncols, args.nrows, args.val)
    args.output.write(str(g))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=_prog_description)
    parser.add_argument('-o', '--output', default=sys.stdout)
    subparsers = parser.add_subparsers(title='grid types')
#                                       description='valid subcommands',
#                                       help='additional help')

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument('ncols', type=int, help="number of columns in the grid (i.e., the x-dimension)")
    common_parser.add_argument('nrows', type=int, help="number of rows in the grid (i.e., the y-dimension)")

    parse_rand_gauss = subparsers.add_parser('gauss',
            parents=[common_parser],
            help="generate grid with Gaussian (normally) distributed values")
    parse_rand_gauss.add_argument('mean', type=float, help="the mean of the Gaussian distribution")
    parse_rand_gauss.add_argument('sd', type=float, help="the standard deviation of the Gaussian distribution")
    parse_rand_gauss.set_defaults(func=gen_rand_gauss)

    parse_rand_unif_real = subparsers.add_parser('uniform-real',
            parents=[common_parser],
            help="generate grid with uniformly distributed real values")
    parse_rand_unif_real.add_argument('a', type=float, help="the lower bound of the distribution range")
    parse_rand_unif_real.add_argument('b', type=float, help="the upper bound of the distribution range")
    parse_rand_unif_real.set_defaults(func=gen_rand_unif_real)

    parse_rand_unif_int = subparsers.add_parser('uniform-int',
            parents=[common_parser],
            help="generate grid with uniformly distributed integer values")
    parse_rand_unif_int.add_argument('a', type=int, help="the lower bound of the distribution range")
    parse_rand_unif_int.add_argument('b', type=int, help="the upper bound of the distribution range")
    parse_rand_unif_int.set_defaults(func=gen_rand_unif_int)

    parse_fixed_val = subparsers.add_parser('constant',
        parents=[common_parser],
        help="generate grid with constant (fixed) value")
    parse_fixed_val.add_argument('val', help="the value for all the cells of the grid")
    parse_fixed_val.set_defaults(func=gen_fixed_val)

    args = parser.parse_args()
    args.func(args)
