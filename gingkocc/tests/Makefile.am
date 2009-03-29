AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

gingkotestsdir = $(prefix)/@PACKAGE@/test

nobase_dist_gingkotests_DATA =  data/dummy1/dummy.conf               \
                                   data/dummy1/grids/cc_0-200.grd       \
                                   data/dummy1/grids/cc_0-50.grd        \
                                   data/dummy1/grids/env1_U005.grd      \
                                   data/dummy1/grids/env1_default.grd   \
                                   data/dummy1/grids/env2_U-005.grd     \
                                   data/dummy1/grids/env2_default.grd   \
                                   data/dummy1/grids/env3_default.grd   \
                                   data/dummy1/grids/spA_barrier.grd    \
                                   data/dummy1/grids/spA_continuous.grd \
                                   data/dummy1/grids/spB_barrier.grd    \
                                   data/dummy1/grids/spB_continuous.grd
 
dist_gingkotests_SCRIPTS = gingko_tests.py \
                   run_test_cmdopts.py \
                   run_test_tree.py \
                   run_test_conf_parser.py \
                   run_test_textutils.py \
                   run_test_filesys.py \
                   run_test_asciigrid.py

gingkotests_PROGRAMS = test_cmdopts \
               test_tree_from_parent_indexes \
               test_tree_from_genealogies \
               test_parse_conf_block \
               test_textutils \
               test_asciigrid \
               test_parse_conf_file \
               test_filesys
               
test_cmdopts_SOURCES = \
    test_cmdopts.cpp \
    ../cmdopts.cpp \
    ../textutils.cpp \
    ../filesys.cpp
    
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
    ../worldconf.cpp
    
test_textutils_SOURCES = \
    test_textutils.cpp \
    ../textutils.cpp
    
test_filesys_SOURCES = \
    test_filesys.cpp \
    ../filesys.cpp 
    
test_asciigrid_SOURCES = \
    test_asciigrid.cpp \
    ../asciigrid.cpp \
    ../textutils.cpp
        
test_parse_conf_file_SOURCES = \
    test_parse_conf_file.cpp \
    ../textutils.cpp \
    ../worldconf.cpp
