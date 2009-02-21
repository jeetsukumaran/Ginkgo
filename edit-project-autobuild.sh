#! /bin/sh

# (1) Adding a File to the Build
# ------------------------------
# 1. In source directory's "Makefile.am"
#    (a) If file is a new executable, add name to "bin_PROGRAMS" list 
#        (space separated)
#    (b) Under "<canonical-progname>_SOURCES", add file. 'canonical-progname'
#        is the same as the <progname> enterd in the "bin_PROGRAMS" list, but
#        with all non ['A'-'Z', 'a'-'z', '0'-'9', '@'] characters substituted
#        by an underscore.
# 2. Run "bootstrap.sh".
#
# (2) Adding a Directory to the Build
# -----------------------------------
# 1. In "configure.ac" of project root, add "rel/path/to/new/dir/Makefile" 
#    under "AC_CONFIG_FILES" list.
# 2. In "Makefile.am" of project root, add "rel/path/to/new/dir" to "SUBDIRS"
#    list.
# 3. Add "Makefile.am" to new source directory.
# 4. Run "bootstrap.sh".

FILES=$(find . -name "configure.ac" -or -name "Makefile.am")

if [ -e "/usr/bin/bbedit" ] 
then
    ED="bbedit --new-window"
elif [ -n $EDITOR ]
then
    ED=$EDITOR
else
    ED=vi
fi

if [ -n "$FILES" ] 
then
    $ED $FILES
else
    echo "no autoconf/automake directive files found"
fi
