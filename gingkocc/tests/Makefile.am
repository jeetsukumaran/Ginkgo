AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

bin_PROGRAMS = test_cmdopts \
               test_tree_from_parent_indexes \
               test_tree_from_genealogies \
               test_conf_parse_block
               
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
    
test_conf_parse_block_SOURCES = \
    test_conf_parse_block.cpp \
    ../textutils.cpp \
    ../conf.cpp