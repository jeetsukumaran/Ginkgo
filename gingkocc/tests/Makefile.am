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
                   run_test_cmdopt.py \
                   run_test_tree.py \
                   run_test_conf_parser.py \
                   run_test_textutil.py \
                   run_test_filesys.py \
                   run_test_asciigrid.py

gingkotests_PROGRAMS = test_cmdopt \
               test_tree_from_parent_indexes \
               test_tree_from_genealogies \
               test_parse_conf_block \
               test_textutil \
               test_asciigrid \
               test_parse_conf_file \
               test_filesys \
               test_convert
               
test_cmdopt_SOURCES = \
    test_cmdopt.cpp \
    ../cmdopt.cpp \
    ../textutil.cpp \
    ../filesys.cpp
    
test_tree_from_parent_indexes_SOURCES = \
    test_tree_from_parent_indexes.cpp \
    ../textutil.cpp \
    ../tree.cpp
    
test_tree_from_genealogies_SOURCES = \
    test_tree_from_genealogies.cpp \
    ../biosys.cpp \
    ../tree.cpp
    
test_parse_conf_block_SOURCES = \
    test_parse_conf_block.cpp \
    ../textutil.cpp \
    ../confsys.cpp \
    ../convert.cpp
    
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
        
test_parse_conf_file_SOURCES = \
    test_parse_conf_file.cpp \
    ../textutil.cpp \
    ../confsys.cpp \
    ../filesys.cpp \
    ../convert.cpp
    
test_convert_SOURCES = \
    test_convert.cpp \
    ../convert.cpp
    