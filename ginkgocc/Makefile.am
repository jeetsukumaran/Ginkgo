AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

bin_PROGRAMS = ginkgo@PROGRAMLABEL@

BUILT_SOURCES = ginkgo_build_id.h
EXTRA_DIST = echo_build_id.sh win_config.h

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

ginkgo_build_id.h: $(HEADERS) $(SOURCES)
	$(top_srcdir)/ginkgocc/echo_build_id.sh > ginkgo_build_id.h
