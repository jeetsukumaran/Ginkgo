AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

bindir = ${prefix}/@PACKAGE@/bin
BUILDTIMESTAMP="`date`"
BUILDDESC="`git symbolic-ref HEAD 2> /dev/null | cut -b 12-`-`git log --pretty=format:\"%h\" -1`"

bin_PROGRAMS = ginkgo@PROGRAMLABEL@

noinst_HEADERS = \
	ginkgo_defs.hpp

ginkgo@PROGRAMLABEL@_CPPFLAGS = -DBUILDDESC=\"${BUILDDESC}\" -DBUILDTIMESTAMP=\"${BUILDTIMESTAMP}\"
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
