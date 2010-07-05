AM_CPPFLAGS = @CPPFLAGS@ -DGINKGO_NUM_LOCI=@NUM_NEUTRAL_LOCI@
AM_LDFLAGS  = @LDFLAGS@

bin_PROGRAMS = ginkgo@PROGRAMLABEL@

BUILT_SOURCES = ginkgo_info.h
EXTRA_DIST = echo_ginkgo_info.sh win_config.h

noinst_HEADERS = \
	ginkgo_defs.hpp \
	asciigrid.hpp \
	cell.hpp \
	cmdopt.hpp \
	confsys.hpp \
	convert.hpp \
	filesys.hpp \
	landscape.hpp \
	logging.hpp \
	memcheck.hpp \
	organism.hpp \
	population.hpp \
	randgen.hpp \
	species.hpp \
	textutil.hpp \
	tree.hpp \
	world.hpp \
	xmlParser.h

ginkgo@PROGRAMLABEL@_SOURCES = \
    ginkgo.cpp \
    randgen.cpp \
    organism.cpp \
    species.cpp \
    population.cpp \
    landscape.cpp \
    cell.cpp \
    world.cpp \
    logging.cpp \
    cmdopt.cpp \
    textutil.cpp \
    tree.cpp \
    confsys.cpp \
    asciigrid.cpp  \
    filesys.cpp \
    xmlParser.cpp \
    memcheck.cpp

ginkgo_info.h: $(HEADERS) $(SOURCES)
	$(top_srcdir)/ginkgocc/echo_ginkgo_info.sh > ginkgo_info.h
