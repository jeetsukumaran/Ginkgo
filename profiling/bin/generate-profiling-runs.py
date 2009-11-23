#! /usr/bin/env python

from optparse import OptionGroup
from optparse import OptionParser
import os

GINKGO_XML = """
<?xml version="1.0"?>
<ginkgo>
    <world label="run_%(run_title)s"
           x_range = "%(size_x)d"
           y_range = "%(size_y)d"
           num_gens = "%(ngens)d"
           fitness_dimensions = "1"
           fitness_grain ="1"
           suppress_final_output = "True"
           default_cell_carrying_capacity = "%(cell_carrying_capacity)d">
        <biota>
            <lineage id="Pseudozoon">
                <selectionWeights>1</selectionWeights>
                <genotypicFitness>0</genotypicFitness>
                <genotypicFitnessMutationRate>0.0</genotypicFitnessMutationRate>
                <genotypicFitnessMutationSize>0</genotypicFitnessMutationSize>
                <fecundity>20</fecundity>
                <seedPopulations>
                    <seedPopulation x="1" y="5" size="10">
                        <ancestralPopulationSize>50</ancestralPopulationSize>
                        <ancestralGenerations>1000</ancestralGenerations>
                    </seedPopulation>
                </seedPopulations>
                <movementCapacity>1</movementCapacity>
                <movementProbability>1.0</movementProbability>
            </lineage>
        </biota>
        <samples>
            <occurrence lineage="Pseudozoon" gen="100" />
            <occurrence lineage="Pseudozoon" gen="1000" />
            <occurrence lineage="Pseudozoon" gen="10000" />
            <occurrence lineage="Pseudozoon" gen="100000" />
            <occurrence lineage="Pseudozoon" gen="1000000" />
            <occurrence lineage="Pseudozoon" gen="10000000" />
        </samples>
    </world>
</ginkgo>
"""

JOB_SGE = """
# /bin/sh
#$ -cwd
#$ -V
#$ -q "general.q"
#$ -l h_vmem=15G
#$ -l vf=15G
#$ -N %(run_title)s
syrupy.py -r -t %(run_title)s_profile -B ginkgo %(run_title)s.xml
"""

def setup_ginkgo_xml(output_dir, run_title, x, y, cc, ngens=1000000):
    params = {
        'run_title' : run_title,
        'size_x' : x,
        'size_y' : y,
        'cell_carrying_capacity' : cc,
        'ngens' : ngens
    }
    gingko = GINKGO_XML % params
    f = open(os.path.join(output_dir, run_title + ".xml"), "w")
    f.write(gingko)
    f.close()

def setup_job(output_dir, run_title):
    params = {
        'run_title' : run_title
    }
    job = JOB_SGE % params
    f = open(os.path.join(output_dir, run_title + ".job"), "w")
    f.write(job)
    f.close()

def setup_run(output_dir, run_title, x, y, cc, ngens=1000000):
    setup_ginkgo_xml(output_dir, run_title, x, y, cc, ngens)
    setup_job(output_dir, run_title)

def sweep_over_pop_sizes(output_dir, begin=10000, end=1280000):
    x = 10
    y = 10
    i = begin
    while i <= end:
        cc = i / (x*y)
        setup_run(output_dir, "P%07d" % i, x, y, cc)
        i = i * 2

def main():
    """
    Main CLI handler.
    """

    parser = OptionParser(add_help_option=True)

    parser.add_option('-o', '--output-dir',
        action='store',
        dest='output_dir',
        type='string',
        default='.',
        metavar='PATH',
        help="path to output directory (default='%default')")

    parser.add_option('--sweep-popsizes',
        action='store_true',
        dest='sweep_popsizes',
        default=False,
        help="generate suite of runs to sweep over population sizes")

    parser.add_option('-t', '--run-title',
        action='store',
        dest='run_title',
        type='string',
        default='GINKGO',
        metavar='NAME',
        help="title of run, to be prefixed to all output files (default='%default')")

    parser.add_option('-x', '--size-x',
        action='store',
        dest='size_x',
        type='int',
        default='10',
        metavar='DIM',
        help="x-dimension size (default=%default)")

    parser.add_option('-y', '--size-y',
        action='store',
        dest='size_y',
        type='int',
        default='10',
        metavar='DIM',
        help="y-dimension size (default=%default)")

    parser.add_option('-c', '--carrying-capacity',
        action='store',
        dest='cc',
        type='int',
        default='1000',
        metavar='K',
        help="cell carrying capacity (default=%default)")

    parser.add_option('-n', '--num-generations',
        action='store',
        dest='ngens',
        type='int',
        default='10000',
        metavar='NUM',
        help="number of generations to run (default=%default)")

    (opts, args) = parser.parse_args()

    if not os.path.exists(opts.output_dir):
        os.makedirs(opts.output_dir)

    if opts.sweep_popsizes:
        sweep_over_pop_sizes(opts.output_dir)
    else:
        setup_run(output_dir=opts.output_dir,
                run_title=opts.run_title,
                x=opts.size_x,
                y=opts.size_y,
                cc=opts.cc,
                ngens=opts.ngens)

if __name__ == '__main__':
    main()



