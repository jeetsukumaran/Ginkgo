*************************************
Setting Up and Executing a Ginkgo Run
*************************************

.. highlight:: xml

Introduction
============

All Ginkgo simulations are driven by an XML-format configuration file.
As with any valid XML document, the configuration file should begin with an XML version declaration (i.e, ``<?xml version="1.0"?>``).
Following this, the configuration file should have a ``<ginkgo>`` element as its root level element.
This element has various attributes that customize logging, reporting and other output of the simulation, and are described below.
The ``<ginkgo>`` element has the following top-level elements which define in the detail the actual simulation and the sampling of genealogies:

    ``<system>``
        This defines details of the simulation framework, such as the the number of dimensions to the fitness equation, the number of generations for which to run the simulation, the global selection strength, etc.

    ``<landscape>``
        This defines the dimensions of the landscape of the simulation.

    ``<lineages>``
        This defines the different species ecologies of the organisms that populate and interact in the landscape.

    ``<initialization>``
        This defines the "burn-in" or "bootstrapping" phase of the simulation, where the landscape is populated and life-cycles are run until the population structures of the organsms are conditioned under the specified environment.

    ``<environments>``
        This defines the changing geographic and environmental templates of the landscapes at different points in the simulation.

    ``<samples>``
        This defines the schedule of genealogies and occurrence samples that takes place during the course of the simulation.

Thus, a typical configuration file would have a structure line this::

    <?xml version="1.0"?>
    <ginkgo>
        <system> ... </system>
        <landscape> ... </landscape>
        <lineages> ... </lineages>
        <initialization> ... </initialization>
        <environments> ... </environments>
        <samples> ... </samples>
    </ginkgo>

The Ginkgo Configuration File
=============================

The ``<ginkgo>`` Element: Defining the Reporting Parameters
-----------------------------------------------------------

This is the (mandatory) root element of the configuration document, and contains as its children all other elements that, taken together, provide that directives that define and drive the simulation, from start to initialization to the main cycles to reporting and output.

This ``<ginkgo>`` can optionally have the following attributes:

    ``title``
        This defines a common prefix for all output files from a run of this simulation. If the "``-i``" command-line options is passed to the Ginkgo invocation, the identifier given by the "``-i``" flag will be appended to the end of the string given by this attribute. Typically, the ``title`` attribute would be used to provide some sort of description of the scenario being simulated, while the command-line "``-i``" flag describes the replicate number or identity.
        This combination of an optional replicate identifier in the Ginkgo invocation and simulation schema ttile allow for the running multiple realizations of the same simulation regime, but maintain distinct sets of output files for each realization in the same directory, as each output file will have the prefix "``<title><i>``".

    ``log_frequency``
        A non-zero positive integer specifying frequency of log output in terms of numbers of generations. Defaults to 10.

    ``multifurcating_trees``
        If "``True``", then multifurcating trees are preserved. If "``False``", then multifurcating trees are converted into bifurcating trees by arbitrarily resolving the polytomies with new edges of length 0. Defaults to "``True``" (i.e., preserve polytomies).

    ``final_output``
        By default, at the end of the simulation Ginkgo samples and reports the genealogies, traits and occurrences of *every* organism of *every* lineage in *every* cell.  Setting ``final_output`` to "``False``" suppresses the automatic reporting of the full results, thus restricting the output to a own custom fine-tuned sampling regime. Note that if ``final_output`` is set to "``False``", if no ``<sample>`` directives are given, then no output will be produced apart from the log files!

    ``full_complement_diploid_trees``
        By default, the genealogies for the diploid locii are built by randomly sampling one of the alleles of the diploid complement for each organism (and saved with the "``.diploid1.tre``" extension). By setting this attribute to "``False``", you can ask that Ginkgo produce a file (which will be saved with the "``.diploid2.tre``"  extension) where *both* alleles are sampled.

For example::

    <ginkgo
        title="results"
        log_frequency="1"
        multifurcating_trees="True"
        final_output = "True"
        full_complement_diploid_trees = "False">
        :
        :
    </ginkgo>


The ``<system>`` Element: Defining the Simulation Meta-Parameters
-----------------------------------------------------------------

The mandatory ``<system>`` element defines some meta-level aspects of the simulation in its child elements:

    ``<ngens>``
        This element is mandatory. It requires a positive non-zero integer} argument which defines the number of generations of cycles that the simulation will run.

    ``<random_seed>``
        Seed for pseudo-random number generator. The command-line specified seed (using the "``-z``" flag) overrides this, and if neither are specified, then the seed defaults to system time.

    ``<fitness_dimensions>``
        This element takes a single non-zero positive integer as its content, and specifies the number of fitness factors (trait types, dimensions of environment variables, etc.) in the simulation. If not specified, it defaults to 10.

    ``<global_selection_strength>``
        This requires a positive real number argument, and weights the overall multi-dimensional Euclidean distance between the vector of an organism's traits and the corresponding environmental optima. A value of 0.0 means that **no** selection takes place. A value of 1.0 results in the survival probability of the organisms given directly by the exponentiated distance, while higher values increase the strength of selection by lowering the survival probability of an organism for a given trait-optima distance. If not specified, defaults to 1.0.

For example::

    <system>
        <fitness_dimensions>1</fitness_dimensions>
        <global_selection_strength>1.0</global_selection_strength>
        <ngens>1000</ngens>
    </system>

The ``<lineages>`` Element: Defining the Biota
----------------------------------------------

This mandatory element defines the ecologies of all species or lineages in the simulation.
It contains one or more ``<lineage>`` child element, with each ``<lineage>`` child element defining the ecology of a single species in the simulation.

A ``<lineage>`` element in turn has one mandatory attribute, ``id``, which serves as the identifier or label for the species it defines, as well as the following child elements:

    ``fitness_trait_default_genotypes``
        The traits of organisms evolve under the specified selection pressures once the simulation begins, but the initial values for the traits for the first generation of organisms of each lineage are specified by this element, the contents of which should be a space-delimited list of numeric values.
        If not given, this defaults to all zeros.

    ``fitness_trait_relative_selection_weights``
        The contents of this element should be a space-delimited list of numeric values that provide the *relative* weights of the selection pressures for each of traits in the fitness function
        For example, "``1 1 1 1 1``" models a niche in which the organism is equally sensitive to all 5 environmental dimensions, while "``5 1 1 1 0``" models a nich in which the organism is extremely sensitive to the first environmental dimension, and, conversely extremely insensitive to the last environmental dimension, while having a moderate sensitivity to the remaining middle three dimensions.
        As these are relative weights, then "``1 1 1 1 1``" and "``2 2 2 2 2``" etc., all describe the same niche.
        If not provided, then this vector defaults to equal-weighting on all environment distances.

    ``fecundity``
        The number of offspring per mating: specified as positive integer value.

    ``movement_capacity``
        The element determines the number of "movement credits" available to an organism at the beginning of the migration phase, specified as a positive integer.
        This element has an attribute, "``distribution``", which can take on one of two values: "constant" or "poisson". If "constant", then all organisms get the same, fixed, movement capacity, equal to this element's value. If "poisson", then organisms get a random movement capacity, drawn from a Poisson distribution with a mean given by this element's value.

For example::

    <lineages>
        <lineage id="Zu">
            <fitness_trait_default_genotypes>0</fitness_trait_default_genotypes>
            <fecundity>16</fecundity>
            <movement_capacity distribution="constant">1</movement_capacity>
        </lineage>
        <lineage id="Zv">
            <fitness_trait_default_genotypes>0</fitness_trait_default_genotypes>
            <fecundity>1</fecundity>
            <movement_capacity distribution="poisson">9</movement_capacity>
        </lineage>
    </lineages>

The ``<initialization>`` Element: Defining the Starting Conditions
------------------------------------------------------------------

This element defines the initialization regime of the simulation.
In the initialization regime, the landscape is seeded with populations of organisms of various species, and the simulation is run under a particular set of environmental conditions, until all tracked neutral locii have coalesced into their respective common ancestors.
This way, the genetic structure and associated genealogies can be calibrated or conditioned under a known scenario or geo-demographic model before the main cycles of the simulation begin.
The ``<initialization>`` element has a ``<populations>`` element and an ``<environment>`` element as child elements::

    <initialization>
        <populations>
        :
        :
        </populations>
        <environment>
        :
        :
        </environment>
    </initialization>

The ``<populations>`` element describes the initial seeding of the cells of the landscape with organisms of the various species.
The ``<populations>`` element contains multiple ``<cell>`` child elements, each of which has an ``x`` and ``y`` attribute specifying the cell position on the landscape.
Each ``<cell>`` element, in turn, has multiple ``<population>`` child elements, each of which defines the number of organisms of a particular lineage to be introduced into that that cell at the beginning of the initialization phase.
For example::

    <populations>
        <cell x="0" y="0">
            <population lineage="Zu" size="100" />
            <population lineage="Zv" size="100" />
        </cell>
        <cell x="0" y="1">
            <population lineage="Zu" size="100" />
            <population lineage="Zv" size="100" />
        </cell>
        <cell x="1" y="0">
            <population lineage="Zu" size="100" />
            <population lineage="Zv" size="100" />
        </cell>
        <cell x="1" y="1">
            <population lineage="Zu" size="100" />
            <population lineage="Zv" size="100" />
        </cell>
        :
        :
    </populations>

The ``<environment>`` element describes the geographical and evolutionary template that structures or conditions the organisms in the initialization phase: the connectivity between cells, the fitness optima of the various fitness dimensions, etc.
This syntax and semantics of this element is identical to the ``<environment>`` child element of the ``<environments>`` element, described in detail below, except that the ``gen`` attribute is ignored.

A complete example of an initialization section is seen here::

    <initialization>
        <populations>
            <cell x="0" y="0">
                <population lineage="Zu" size="100" />
                <population lineage="Zv" size="100" />
            </cell>
            <cell x="0" y="1">
                <population lineage="Zu" size="100" />
                <population lineage="Zv" size="100" />
            </cell>
            <cell x="1" y="0">
                <population lineage="Zu" size="100" />
                <population lineage="Zv" size="100" />
            </cell>
            <cell x="1" y="1">
                <population lineage="Zu" size="100" />
                <population lineage="Zv" size="100" />
            </cell>
        </populations>
        <environment>
            <carrying_capacity>cc100.asc</carrying_capacity>
            <movement_costs lineage="Zu">g25x25_mv1.asc</movement_costs>
            <movement_costs lineage="Zv">g25x25_mv9.asc</movement_costs>
            <fitness_trait_optima trait="0">trait_unif0.asc</fitness_trait_optima>
            <fitness_trait_optima trait="1">trait_rand1.asc</fitness_trait_optima>
            <fitness_trait_optima trait="2">trait_rand2.asc</fitness_trait_optima>
        </environment>
    </initialization>

The ``<environments>`` Element: Defining the Landscape and Climatic History
---------------------------------------------------------------------------

The ``<environments>`` element controls the geographical template (i.e., the connectivity and movement costs of the landscape), carrying capacities and fitness regimes over the course of the simulation.
This element consists of a list of ``<Environment>`` elements, with each ``<environment>`` element encapsulating the suite of movement costs, carrying capacities, fitness trait optima, etc. to be activate at a particular generation.

The ``<environment>`` element has a single mandatory attribute, ``gen``, which takes a positive integer value specifying the generation number that this element's directive should take affect (with the first generation of the simulation being generation 0). For example::

    <environments>
        <environment gen="100"> ... </environment>
        <environment gen="200"> ... </environment>
    </environments>

The child elements of the ``<environment>`` element can be one of the following:

    ``<carrying_capacity>``
        This element specifies the carrying capacity of each cell of the landscape, i.e. the maximum number of organisms that can occupy the cell at the end of each generation. The content of this element should be an ESRI ASCII format grid file of the same dimensions as the landscape, with the cell values specifying the maximum carrying capacity of the corresponding cell on the landscape. If this element is not included for a particular ``<environment>`` suite, then the previous setting is retained.


    ``<movement_costs>``
        This element specifies the cell entry costs of the landscape for a particular lineage or species.
        It has a single mandatory attribute, ``lineage``, which should be the label or identifier of the species for which these costs should be applied.
        The content of this element should be an ESRI ASCII format grid file of the same dimensions as the landscape, with the cell values specifying the entry costs for the correspondng cells on the landscape. If this element is not included for a particular ``<environment>`` suite, then the previous setting is retained.


    ``<fitness_trait_optima>``
        This element specifies the optimum environmental value for each cell of the landscape for a particular trait.
        It has a single mandatory attribute, ``trait``, which should be the 0-based index of the trait.
        The content of this element should be an ESRI ASCII format grid file of the same dimensions as the landscape, with the cell values specifying the optimum trait value for the correspondng cells on the landscape. If this element is not included for a particular ``<environment>`` suite, then the previous setting is retained.

The following is an example of a more thoroughly specified environmental regime::

    <environments>
        <environment gen="100">
            <movement_costs lineage="Sp1">mov_sp1_gen100.asc</movement_costs>
            <movement_costs lineage="Sp2">mov_sp1_gen100.asc</movement_costs>
            <fitness_trait_optima trait="0">traits_x.asc</fitness_trait_optima>
            <fitness_trait_optima trait="1">traits_y.asc</fitness_trait_optima>
        </environment>
        <environment gen="500">
            <movement_costs lineage="Sp1">mov_sp1_gen500.asc</movement_costs>
            <carrying_capacity>cc_gen500.asc</carrying_capacity>
        </environment>
        <environment gen="1000">
            <carrying_capacity>cc_gen1000.asc</carrying_capacity>
            <movement_costs lineage="Sp2">mov_sp2_gen1000.asc</movement_costs>
            <fitness_trait_optima trait="0">traits_a.asc</fitness_trait_optima>
            <fitness_trait_optima trait="1">traits_b.asc</fitness_trait_optima>
        </environment>
    </environments>

In generation 100 of the simulation, cell entry costs are set for species "Sp1" and "Sp2", as well as the environmental fitness optima for the first and second fitness factors. Until these settings were changed in generation "100", they would have retained the values specified in the initialization phase. In addition, settings not specified in this generation (for example, the carrying capacity) would also retain the values given in the initialization phase. In generation 500, the movement costs for "Sp1" and the overall carrying capacity of landscape are changed. Again, all other settings (e.g., the movement costs for "Sp1", the fitness trait optima) remain unchanged from the values they had in generation "100".

The ``<samples>`` Element: Defining the Sampling Regime
-------------------------------------------------------

The ``<samples>`` element serves as container for one or more ``<sample>`` elements, each of which describes the sampling design for a set of genealogies and occurrence data that will be sampled at a particular point in time during the simulation for a particular lineage.
Each ``<sample>`` element has two mandatory attributes: the ``gen`` attribute is a positive integer that defines the 0-based index of the generation number that this sample should be taken, and the ``lineage`` attribute which is a string giving the identifier of the lineage to be sampled::

    <samples>
        <sample gen="100" lineage="Sp1"> ... </sample>
        <sample gen="100" lineage="Sp2"> ... </sample>
        <sample gen="200" lineage="Sp1"> ... </sample>
        :
    </samples>

Each sample, by default, will result in three files:

    * an occurrence matrix, showing the number of individuals of the specified lineage in each cell of the landscape
    * a genealogy showing the relationship of the alleles of the neutral haploid locus of every organism of the specified lineage from every cell of the landscape
    * a genealogy for a random allele sampled from each of the diploid loci of every organism of the specified lineage from every cell of the landscape
    * a traits file

The occurrence matrix will be given as an ESRI ASCII format integer grid.
The genealogies will be given as NEXUS-format tree files, with the haploid locus genealogy file having a single tree, and the diploid locus genealogy file have multiple trees (one tree per diploid locus).
The traits file will be a NEXUS-format character matrix, with one column per fitness factor or dimension, and a final column showing the current fitness of the sampled individual.

There is an exact correspondence between individuals in the genealogy files and the traits files **for a particular sampling**, and, within the same sample, the same individuals will have the same taxon labels across the files.
The taxon labels will be of the following form:

    ``<LINEAGE_ID>_x<X>_y<Y>_<UNIQUE#>``

So, for example, given a lineage named "Sp1", sampled from a cell with coordinates ``(25,5)``, the corresponding taxon label might be:

    ``Sp1_x25_y5_12001``

In many cases, it is not neccessary to sample every individual from every cell in the landscape. The number of individuals sampled from each cell can be restricted using the ``<individuals_per_cell>`` child element of the ``<sample>`` element, which takes a positive number as a value. For example, to limit the samples to 10 individuals per cell::

    <samples>
        <sample gen="100" lineage="Sp1">
            <individuals_per_cell>10</individuals_per_cell>
        </sample>
        :
    </samples>

It is also possible to restrict the sampling to particular cells, using the ``<cells>`` child element of ``<sample>``, which takes in turn one or more ``<cell>`` elements, whose ``x`` and ``y`` attributes specify the coordinates of the cells to be sampled::

    <samples>
        <sample gen="100" lineage="Sp1">
            <cells>
                <cell x="5" y="5" />
                <cell x="6" y="5" />
                <cell x="7" y="5" />
                <cell x="8" y="5" />
                <cell x="9" y="5" />
                <cell x="10" y="5" />
            </cells>
        </sample>
        :
    </samples>

Of course, you can can combined the ``<individuals_per_cell>`` and the ``<cell>`` directives to restrict both the numbers of individuals sampled as well as the cells sampled::

    <samples>
        <sample gen="100" lineage="Sp1">
            <individuals_per_cell>10</individuals_per_cell>
            <cells>
                <cell x="5" y="5" />
                <cell x="6" y="5" />
                :
            </cells>
        </sample>
        :
    </samples>

In some cases, all that might required is the occurrence or incidence data (i.e., the number of individuals of a particular species in each cell of the landscape at particular generation), and not the genealogies.
As the calculation of genealogies are usually relatively time-consuming and computationally-expensive, this procedure can be skipped when not neccessary by setting the ``trees`` attribute of the ``<sample>`` element to "``False``"::

    <samples>
        <sample gen="0" lineage="Sp1" trees="False" />
        <sample gen="50" lineage="Sp1" trees="False" />
        <sample gen="100" lineage="Sp1">
            <individuals_per_cell>10</individuals_per_cell>
            <cells>
                <cell x="5" y="5" />
                <cell x="6" y="5" />
                :
            </cells>
        </sample>
        <sample gen="200" lineage="Sp1" trees="False" />
        :
    </samples>

By default, all output files produced by a ``<sample>`` directive will have a filename prefix in the following form:

    ``<GINKGO-TITLE>_G<GENERATION#>_<LINEAGE-ID>``

So, for example, with if the ``title`` attribute of the ``<ginkgo>`` element is "run1", the following sampling directive::

    <sample gen="10000" lineage="Sp1" />

would produce the following files:

    ``run1_G00010000_Sp1_occurrences.asc``
        The ESRI ASCII grid showing the number of individuals of lineage "Sp1" on the landscape.
    ``run1_G00010000_Sp1.haploid.tre``
        The genealogy for the haploid locus of organisms of lineage "Sp1".
    ``run1_G00010000_Sp1.diploid1.tre``
        The genealogy for the 10 diploid locii of organisms of lineage "Sp1".
    ``run1_G00010000_Sp1.traits.nex``
        The character matrix summarizing the fitness trait values as well as fitness score for organisms of lineage "Sp1".

The ``<sample>`` element takes a ``label`` attribute that gets added to the filename prefix.
Thus, if the sampling directive given above were modified to include a ``label`` attribute set to "pre-climate-change"::

    <sample gen="10000" lineage="Sp1" label="pre-climate-change" />

then the following files would be produced instead:

    - ``run1_G00010000_Sp1_occurrences_pre-climate-change.asc``
    - ``run1_G00010000_Sp1_pre-climate-change.haploid.tre``
    - ``run1_G00010000_Sp1_pre-climate-change.diploid1.tre``
    - ``run1_G00010000_Sp1_pre-climate-change.traits.nex``

Running Ginkgo
==============

Once a configuration file is ready, executing Ginkgo is simply a matter of invoking the executable with the configuration file as an argument::

    $ ginkgo scenario.xml

