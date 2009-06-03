# This file is processed by autoconf to create a Makefile for the GARLI
# This file was written by hand (by mth using PAUP's configure.ac and http://autotoolset.sourceforge.net/tutorial.html as guides).

AC_PREREQ(2.59)
AC_INIT([gingko], [2.5], [jeet@ku.edu])

# AC_CONFIG_SRCDIR (UNIQUE-FILE-IN-SOURCE-DIR) UNIQUE-FILE-IN-SOURCE-DIR
# is some file that is in the package's source directory; `configure'
# checks for this file's existence to make sure that the directory that it
# is told contains the source code in fact does.
AC_CONFIG_SRCDIR([gingkocc/gingko.cpp])
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

# Check whether asserts should be allowed.
AC_ARG_ENABLE(asserts,
	AC_HELP_STRING([--enable-asserts],[build with asserts on (NDEBUG not defined)]), , [enable_asserts=no])
if test "$enable_asserts" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: compiling with assertions on (NDEBUG not defined)
	   ])
fi

# Compile openMP multithreaded version
#AC_ARG_ENABLE(openmp,
#	AC_HELP_STRING([--enable-openmp],[build OpenMP multithreaded version]), , [enable_openmp=no])
#if test "$enable_openmp" = yes; then
#	AC_MSG_NOTICE([
#	   *** NOTE: compiling OpenMP multithreaded version
#	   ])
#fi

# Mpi run forking version
#AC_ARG_ENABLE(mpi,
#	AC_HELP_STRING([--enable-mpi],[build MPI run distributing version]), , [enable_mpi=no])
#if test "$enable_mpi" = yes; then
#	AC_MSG_NOTICE([
#	   *** NOTE: compiling MPI run distributing version
#	   ])
#fi

# old Mpi wrapper version - run series of configs named run0.conf, run1.conf, etc
#AC_ARG_ENABLE(oldmpi,
#        AC_HELP_STRING([--enable-oldmpi],[build old MPI batch run version (YOU PROBABLY DON'T WANT THIS)]), , [enable_oldmpi=no])
#if test "$enable_oldmpi" = yes; then
#        AC_MSG_NOTICE([
#           *** NOTE: compiling old MPI batch run version
#           ])
#fi

# single precision version (still being tested)
# AC_ARG_ENABLE(single-prec,
#         AC_HELP_STRING([--enable-single-prec],[use single precision floating point variables (EXPERIMENTAL)]), , [enable_single_prec=no])
# if test "$enable_single_prec" = yes; then
#         AC_MSG_NOTICE([
#            *** NOTE: compiling single precision floating point version
#            ])
# fi

# Initialize optimization flag in case it doesn't get set below.
CXXFLAGS_OPTIM_SPEED="-O"

# "-g" may not work with some compilers, but end users shouldn't be
if test "$enable_debugging" = yes; then
	CXXFLAGS_OPTIM_SPEED="-O0"
	CXXFLAGS="$CXXFLAGS -Wall -pedantic -ansi -g -Wreturn-type -Wunused -Wredundant-decls -Wcast-align -Wcomment -Wextra"
#	CXXFLAGS="$CXXFLAGS -Werror"
fi

#DJZ allow turning off NDEBUG only
if test  "$enable_asserts" = no; then
	CXXFLAGS="$CXXFLAGS -DNDEBUG"
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

AC_LANG(C++)
CXXFLAGS="$CXXFLAGS $CXXFLAGS_OPTIM_SPEED"


################################################################################
#	By default, use the bundled ncl code from source, 
#		but let the --with-ncl argument override this.
################################################################################
#NCL_INC_DIR="/usr/local/include"
#NCL_LIB_DIR="/usr/local/lib"
#AC_ARG_WITH(
#	[ncl],
#	AC_HELP_STRING(
#		[--with-ncl=DIR],
#		[Specify the root directory for ncl library (parent of the include/ncl and lib directories). Omit this flag to build and use the bundled version of ncl.]
#		),
#	[
#	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
#		NCL_INC_DIR="$withval/include"
#		NCL_LIB_DIR="$withval/lib"
#	fi
#	])

#if ! test -d "$NCL_INC_DIR/ncl" ; then
#	AC_MSG_ERROR([NCL is a prerequisite for building Garli. Expecting a directory called $NCL_INC_DIR/ncl, verify that the argument to --with-ncl is correct. You may need to manually make the directory $NCL_INC_DIR/ncl copy all of the ncl headers into that directory if you are using NCL 2.0])
#else
#	if ! test -f "$NCL_INC_DIR/ncl/nxsdefs.h" ; then
#		AC_MSG_ERROR([NCL is a prerequisite for building Garli. Expecting a directory called $NCL_INC_DIR/ncl, verify that the argument to --with-ncl is correct. You may need to manually make the directory $NCL_INC_DIR/ncl copy all of the ncl headers into that directory if you are using NCL 2.0])
#	fi
#fi

#if ! test -d "$NCL_LIB_DIR" ; then
#	AC_MSG_ERROR([NCL is a prerequisite for building Garli. Expecting a directory called $NCL_LIB_DIR, verify that the argument to --with-ncl is correct.])
#fi

LIBS="$LIBS"
AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])


# Variables to substitute in Makefile:
AC_SUBST(CPPFLAGS)

AC_CONFIG_FILES([
	Makefile 			\
	gingkocc/Makefile   \
	gingkocc/tests/Makefile		\
	])

AC_OUTPUT
