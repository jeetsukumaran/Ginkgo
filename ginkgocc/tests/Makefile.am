AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

ginkgotestsdir = $(prefix)/@PACKAGE@/test
 
dist_ginkgotests_SCRIPTS = ginkgo_tests.py \
                   run_test_cmdopt.py \
                   run_test_tree.py \
                   run_test_conf_parser.py \
                   run_test_textutil.py \
                   run_test_filesys.py \
                   run_test_asciigrid.py

ginkgotests_PROGRAMS = test_cmdopt \
               test_tree_from_parent_indexes \
               test_tree_from_genealogies \
               test_parse_conf_block \
               test_textutil \
               test_asciigrid \
               test_parse_conf_file \
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
    ../biosys.cpp \
    ../tree.cpp \
    ../memcheck.cpp
    
test_tree_from_genealogies_SOURCES = \
    test_tree_from_genealogies.cpp \
    ../biosys.cpp \
    ../tree.cpp \
    ../memcheck.cpp
    
test_parse_conf_block_SOURCES = \
    test_parse_conf_block.cpp \
    ../textutil.cpp \
    ../confsys.cpp \
    ../convert.cpp \
    ../randgen.cpp \
    ../world.cpp  \
    ../cell.cpp \
    ../landscape.cpp \
    ../biosys.cpp   \
    ../asciigrid.cpp \
    ../filesys.cpp\
    ../tree.cpp \
    ../xmlParser.cpp \
    ../memcheck.cpp    
        
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
    ../convert.cpp \
    ../randgen.cpp \
    ../world.cpp \
    ../cell.cpp \
    ../landscape.cpp \
    ../biosys.cpp \
    ../asciigrid.cpp\
    ../tree.cpp \
    ../xmlParser.cpp \
    ../memcheck.cpp
    
test_convert_SOURCES = \
    test_convert.cpp \
    ../convert.cpp \
    ../textutil.cpp 
    
test_world_log_SOURCES = \
    test_world_log.cpp \
    ../world.cpp \
    ../randgen.cpp \
    ../biosys.cpp \
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
