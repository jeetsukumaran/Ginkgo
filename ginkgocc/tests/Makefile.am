AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

ginkgotestsdir = $(prefix)/@PACKAGE@/test

dist_ginkgotests_SCRIPTS = ginkgo_tests.py \
                   run_test_cmdopt.py \
                   run_test_tree.py \
                   run_test_textutil.py \
                   run_test_filesys.py \
                   run_test_asciigrid.py

ginkgotests_PROGRAMS = test_cmdopt \
               test_tree_from_parent_indexes \
               test_tree_from_genealogies \
               test_organism_states \
               test_textutil \
               test_asciigrid \
               test_filesys \
               test_convert \
               test_world_log

test_cmdopt_SOURCES = \
    test_cmdopt.cpp \
    ../cmdopt.cpp \
    ../textutil.cpp \
    ../filesys.cpp

test_tree_from_parent_indexes_SOURCES = \
    test_tree_from_parent_indexes.cpp \
    ../textutil.cpp \
    ../organism.cpp \
    ../species.cpp \
    ../population.cpp \
    ../tree.cpp \
    ../memcheck.cpp

test_tree_from_genealogies_SOURCES = \
    test_tree_from_genealogies.cpp \
    ../organism.cpp \
    ../species.cpp \
    ../population.cpp \
    ../tree.cpp \
    ../memcheck.cpp

test_organism_states_SOURCES = \
    test_organism_states.cpp \
    ../organism.cpp \
    ../species.cpp \
    ../population.cpp

test_textutil_SOURCES = \
    test_textutil.cpp \
    ../textutil.cpp

test_filesys_SOURCES = \
    test_filesys.cpp \
    ../filesys.cpp

test_asciigrid_SOURCES = \
    test_asciigrid.cpp \
    ../asciigrid.cpp \
    ../textutil.cpp

test_convert_SOURCES = \
    test_convert.cpp \
    ../convert.cpp \
    ../textutil.cpp

test_world_log_SOURCES = \
    test_world_log.cpp \
    ../world.cpp \
    ../randgen.cpp \
    ../organism.cpp \
    ../species.cpp \
    ../population.cpp \
    ../landscape.cpp \
    ../cell.cpp \
    ../cmdopt.cpp \
    ../textutil.cpp \
    ../tree.cpp \
    ../confsys.cpp \
    ../asciigrid.cpp  \
    ../filesys.cpp \
    ../xmlParser.cpp \
    ../memcheck.cpp
