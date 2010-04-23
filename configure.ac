# This file is processed by autoconf to create a Makefile for the GARLI
# This file was written by hand (by mth using PAUP's configure.ac and http://autotoolset.sourceforge.net/tutorial.html as guides).

AC_PREREQ(2.59)
AC_INIT([ginkgo], [3.9.0], [jeet@ku.edu])

# AC_CONFIG_SRCDIR (UNIQUE-FILE-IN-SOURCE-DIR) UNIQUE-FILE-IN-SOURCE-DIR
# is some file that is in the package's source directory; `configure'
# checks for this file's existence to make sure that the directory that it
# is told contains the source code in fact does.
AC_CONFIG_SRCDIR([ginkgocc/ginkgo.cpp])
AC_CONFIG_HEADER([config.h])

#	Directory that contains install-sh and other auxiliary files
AC_CONFIG_AUX_DIR([config])

#   Set default install prefix
AC_PREFIX_DEFAULT([/opt])

################################################################################
# Checks for build-platform and target info
#	this defines the "target" variable that is used later in this file
################################################################################
AC_CANONICAL_TARGET

################################################################################
#	According to (http://www.mail-archive.com/autoconf@gnu.org/msg14232.html)
#		this macro should be after AC_INIT but before AM_INIT_AUTOMAKE
################################################################################
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.9 foreign dist-tarZ tar-ustar filename-length-max=299])

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP

# Initialize CXXFLAGS to prevent it from defaulting to "-g -O2"
USER_CXXFLAGS="$CXXFLAGS"
AC_PROG_CXX
CXXFLAGS="$USER_CXXFLAGS"

AC_PROG_INSTALL
# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h malloc.h stddef.h stdlib.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memmove memset pow sqrt strchr strdup strtol])


#--------------------------------------------------------------------------------------------------#
#  Set CXXFLAGS                                                                                    #
#--------------------------------------------------------------------------------------------------#

# Check for debugging mode.
AC_ARG_ENABLE(debugging,
	AC_HELP_STRING([--enable-debugging],[build for debugging]), , [enable_debugging=no])
if test "$enable_debugging" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: debugging is enabled; optimization is suppressed!
	   ])
fi

# Check for memcheck mode.
AC_ARG_ENABLE(memcheck,
	AC_HELP_STRING([--enable-memcheck],[build for memory-leak checking]), , [enable_memcheck=no])
if test "$enable_memcheck" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: memory checks will be run
	   ])
fi

# Check whether asserts should be allowed.
AC_ARG_ENABLE(asserts,
	AC_HELP_STRING([--enable-asserts],[build with asserts on (NDEBUG not defined)]), , [enable_asserts=no])
if test "$enable_asserts" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: compiling with assertions on (NDEBUG not defined)
	   ])
fi

# Check for mini mode.
AC_ARG_ENABLE(mini,
	AC_HELP_STRING([--enable-mini],[build mini version, landscape size limited to 255 cells]), , [enable_mini=no])
if test "$enable_mini" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: MINI mode: landscape limited to 255 cells
	   ])
fi

# program name decorator

PROGRAMLABEL=""
AC_ARG_WITH([label],
	AC_HELP_STRING([--with-label],[decorate program name with label]), , [with_label=""])
if test -n "$with_label"; then
	AC_MSG_NOTICE([
	   *** NOTE: adding "$with_label" to program name
	   ])
   PROGRAMLABEL="${with_label}"
fi
AC_SUBST(PROGRAMLABEL)

# Initialize optimization flag in case it doesn't get set below.
CXXFLAGS_OPTIM_SPEED="-O"

if test  "$enable_asserts" = no; then
	CXXFLAGS="$CXXFLAGS -DNDEBUG"
fi

if test  "$enable_memcheck" = yes; then
	CXXFLAGS="$CXXFLAGS -DMEMCHECK"
fi

if test  "$enable_mini" = yes; then
	CXXFLAGS="$CXXFLAGS -DMINI"
fi

if test "$CC" = "icc" -o "$CC" = "icc" ; then
#	Intel C compiler for Linux
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-O2 -ip -funroll-loops -fno-alias"
	fi
elif test "$CC" = "ccc"; then
#	Compaq C compiler for Linux
	if test "x$arch" = "x"; then
		arch="host"
	fi
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-fast -inline speed -arch $arch"
	fi
elif test "$CC" = "xlc"; then
#	IBM XL C compiler
	CCFLAGS="$CXXFLAGS -qsourcetype=c++ -qenablevmx -qchars=signed"
	if test "x$arch" = "x"; then
		arch="auto"
	fi
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-O3 -qarch=$arch -qtune=$arch -qalias=ansi -qunroll=yes"
	fi
elif test "x$GCC" = "xyes" ; then
	CXXFLAGS="$CXXFLAGS -Wno-uninitialized"
	if test "$enable_debugging" = yes; then
		CXXFLAGS_OPTIM_SPEED="-O0 -Wimplicit"
	else
		CXXFLAGS_OPTIM_SPEED="-O3 -ffast-math -funroll-loops -fstrict-aliasing"
	fi
	case "$build_os" in
	 	darwin*) CXXFLAGS="$CXXFLAGS" ;;
	 	*) CXXFLAGS="$CXXFLAGS -fsigned-char";;
	esac
fi

#single precision
# if test  "$enable_single_prec" = yes; then
# 	CXXFLAGS="$CXXFLAGS -DSINGLE_PRECISION_FLOATS"
# 	if test "$CC" = "gcc" ; then
# 		CXXFLAGS="$CXXFLAGS -fno-caller-saves"
# 	fi
# fi

#                                                        #
##########################################################

# -DMARKUP_STDCONV = flag for CxmlParser

AC_LANG(C++)
CXXFLAGS="$CXXFLAGS $CXXFLAGS_OPTIM_SPEED -DMARKUP_STDCONV"

#vcs_revision_desc="$(git symbolic-ref HEAD 2> /dev/null | cut -b 12-)-$(git log --pretty=format:'%h, %ad' -1)"
#if test -n "${vcs_revision_desc}" ; then
#    AC_DEFINE_UNQUOTED(VCS_REVISION_DESC,"${vcs_revision_desc}", [VCS revision info])
#else
#    AC_DEFINE(VCS_REVISION_DESC, "", [VCS revision info])
#fi

LIBS="$LIBS"
AC_SUBST([LDFLAGS])

# Variables to substitute in Makefile:
AC_SUBST(CPPFLAGS)

AC_CONFIG_FILES([
	Makefile 			\
	ginkgocc/Makefile   \
	])

AC_OUTPUT
