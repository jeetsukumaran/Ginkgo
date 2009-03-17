#! /usr/bin/env python

###############################################################################
##
## GINGKO Biogeographical Evolution Simulator.
##
## Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
## This program is free software; you can redistribute it and#or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License along
## with this program. If not, see <http:##www.gnu.org#licenses#>.
##
###############################################################################

import os
import logging
import unittest
import re

###############################################################################
## LOGGER SET UP

_LOGGING_LEVEL_ENVAR="GINGKO_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR="GINGKO_LOGGING_FORMAT"

def get_logging_level():
    if _LOGGING_LEVEL_ENVAR in os.environ:
        if os.environ[_LOGGING_LEVEL_ENVAR].upper() == "NOTSET":
            level = logging.NOTSET
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "DEBUG":
            level = logging.DEBUG 
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "INFO":
            level = logging.INFO
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "WARNING":
            level = logging.WARNING
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "ERROR":
            level = logging.ERROR
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
    else:
        level = logging.NOTSET  
    return level        

def get_logger(name="gingko"):
    """
    Returns a logger with name set as given, and configured 
    to the level given by the environment variable _LOGGING_LEVEL_ENVAR.
    """
    logger_set = False    
    logger = logging.getLogger(name)            
    if not logger_set:    
        level = get_logging_level()
        rich_formatter = logging.Formatter("[%(asctime)s] %(filename)s (%(lineno)d): %(levelname) 8s: %(message)s")
        simple_formatter = logging.Formatter("%(levelname) 8s: %(message)s")
        raw_formatter = logging.Formatter("%(message)s")
        default_formatter = simple_formatter
        logging_formatter = default_formatter
        if _LOGGING_FORMAT_ENVAR in os.environ:
            if os.environ[_LOGGING_FORMAT_ENVAR].upper() == "RICH":
                logging_formatter = rich_formatter
            elif os.environ[_LOGGING_FORMAT_ENVAR].upper() == "SIMPLE":
                logging_formatter = simple_formatter
            elif os.environ[_LOGGING_FORMAT_ENVAR].upper() == "NONE":
                logging_formatter = None
            else:
                logging_formatter = default_formatter
        else:
            logging_formatter = default_formatter
        if logging_formatter is not None:            
            logging_formatter.datefmt='%H:%M:%S'
        logger.setLevel(level)
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(logging_formatter)
        logger.addHandler(ch)
    return logger

###############################################################################
## TESTS RUN    

def get_test_suite():
    """
    Creates a unittest.TestSuite from all of the modules in
    `dendropy.tests`. Right now, assumes (a) no subdirectories (though
    this can easily be accomodated) and (b) every test to be run is
    sitting in a module with a file name of 'test*.py', and, conversely,
    every file with a name of 'test*.py' has test(s) to be run.
    """
    # get list of test file names'
    path = os.path.dirname(__file__)  
    files = os.listdir(path)                               
    test_file_pattern = re.compile("test.*\.py$", re.IGNORECASE)   
    test_files = []
    for f in files:
        if test_file_pattern.search(f):
            test_files.append(os.path.splitext(f)[0])
            
    # extract the tests            
    tests = unittest.defaultTestLoader.loadTestsFromNames(test_files)
    
    # return the suite
    return unittest.TestSuite(tests) 

def run():
    "Runs all of the unittests"
    runner = unittest.TextTestRunner()
    runner.run(get_test_suite())

if __name__ == "__main__":
    run()
    