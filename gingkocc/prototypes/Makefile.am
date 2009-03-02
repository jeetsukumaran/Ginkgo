AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

# Sometimes a Makefile variable name is derived from some text the
# maintainer supplies. For instance, a program name listed in `_PROGRAMS'
# is rewritten into the name of a `_SOURCES' variable. In cases like this,
# Automake canonicalizes the text, so that program names and the like do
# not have to follow Makefile variable naming rules. All characters in the
# name except for letters, numbers, the strudel (@), and the underscore
# are turned into underscores when making variable references.
# For example, if your program is named sniff-glue, the derived variable
# name would be sniff_glue_SOURCES, not sniff-glue_SOURCES. Similarly the
# sources for a library named libmumble++.a should be listed in the
# libmumble___a_SOURCES variable.

bin_PROGRAMS = version2

version2_SOURCES = \
    version2.cpp
