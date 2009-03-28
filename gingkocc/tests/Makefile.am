AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

gingkotestsdir = $(prefix)/@PACKAGE@/test
gingkotestsdatadir = $(gingkotestsdir)/data

dist_gingkotestsdata_DATA = data/simple1.grd \
                           data/simple2.grd
 
dist_gingkotests_SCRIPTS = gingko_tests.py \
                   run_test_cmdopts.py \
                   run_test_tree.py \
                   run_test_conf_parser.py \
                   run_test_textutils.py

gingkotests_PROGRAMS = test_cmdopts \
               test_tree_from_parent_indexes \
               test_tree_from_genealogies \
               test_parse_conf_block \
               test_textutils \
               test_asciigrid
               
test_cmdopts_SOURCES = \
    test_cmdopts.cpp \
    ../cmdopts.cpp \
    ../textutils.cpp 
    
test_tree_from_parent_indexes_SOURCES = \
    test_tree_from_parent_indexes.cpp \
    ../textutils.cpp \
    ../tree.cpp
    
test_tree_from_genealogies_SOURCES = \
    test_tree_from_genealogies.cpp \
    ../biosys.cpp \
    ../tree.cpp
    
test_parse_conf_block_SOURCES = \
    test_parse_conf_block.cpp \
    ../textutils.cpp \
    ../conf.cpp
    
test_textutils_SOURCES = \
    test_textutils.cpp \
    ../textutils.cpp 

test_asciigrid_SOURCES = \
    test_asciigrid.cpp \
    ../asciigrid.cpp  \    
    ../textutils.cpp
        