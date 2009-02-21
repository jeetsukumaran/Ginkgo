#! /bin/sh

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
