#######################################################
Ginkgo Phylogeographical Evolution Simulator User Guide
#######################################################

Introduction
============

Ginkgo is a agent-based forward-time simulation to produce gene genealogies and incidence (occurrence) data for multiple populations of multiple species in spatially-explicit and environmentally-heterogenous framework.

Features of Ginkgo include:

    - a virtual landscape of arbitrary dimensions
    - large populations of multiple species of diploid sexually-reproducing (male/female) individuals
    - forward-time spatially-explicit tracking of 10 independent diploid loci and one maternally-inherited locus
    - species-specific, spatially- and temporally-heterogenous connectivity between cells of the landscape (i.e., migration rates can be different for different species in different parts of the landscape at different times)
    - species-specific, spatially- and temporally-heterogenous selection regime (i.e., different species can be more or less sensitive to different aspects of the environment at different times)
    - carrying-capacity driven competition between individuals, with relative fitness a function of inheritable phenotypes and the (spatially- and temporally-dynamic) environment

The entire framework is fully configurable by the user, and can be set to accomodate scenarios that approximate completely neutral classic Wright-Fisher population conditions at the one extreme, to hyper-complex conditions with full micro-, meso-, and macro-scale spatial structuring and selection under multi-level hetereogeneous environmental regimes at the other.

Three kinds of output are produced by Ginkgo:

    - genealogies for each of the independent loci tracked (10 diploid and 1 haploid), along with spatial history (i.e., the position on the landscape occupied not only by the current generation of alleles, but all ancestral alleles as well)
    - incidence or occurrence data, i.e., the numbers of individuals of each species occupying each cell of the landscape
    - fitness trait values and fitness scores for each individual

More Information
================

.. toctree::
    :maxdepth: 2

    installing.rst
    running.rst
