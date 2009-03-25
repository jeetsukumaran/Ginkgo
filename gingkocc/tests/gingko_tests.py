#! /usr/bin/env python

###############################################################################
##
## GINGKO Biogeographical Evolution Simulator.
##
## Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
## This prog_name is free software; you can redistribute it and#or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This prog_name is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License along
## with this prog_name. If not, see <http:##www.gnu.org#licenses#>.
##
###############################################################################

import os
import logging
import unittest
import re
import subprocess
import sys

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
        default_formatter = logging.Formatter("%(name)s: (%(levelname)s) %(message)s")
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
## FILE PATHS

GINGKO_BIN_PATH_ENVAR = "GINGKO_BIN_PATH"

# from: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(prog_name):
    import os
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(prog_name)
    if fpath:
        if is_exe(prog_name):
            return prog_name
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, prog_name)
            if is_exe(exe_file):
                return exe_file
    return None

def get_gingko_bin_path():
    gingko_bin_path = None
    if GINGKO_BIN_PATH_ENVAR in os.environ:
        gingko_bin_path = os.environ[GINGKO_BIN_PATH_ENVAR]
        if not os.path.exists(gingko_bin_path):
            sys.stderr.write('Path "%s" specified by environmental variable "%s" does not exist.\n' % (gingko_bin_path, GINGKO_BIN_PATH_ENVAR))
            sys.exit(1)
    else:
        gingko_bin_path = os.path.dirname(__file__)
    return gingko_bin_path            
            
def get_gingko_program_path(prog_name):
    bin_path = get_gingko_bin_path()
    if bin_path is not None:
        prog_path = os.path.join(get_gingko_bin_path(), prog_name)
        if not os.path.exists(prog_path):
            sys.stderr.write('Program not found: "%s"\n' % prog_path)
            sys.exit(1)
        elif os.path.isdir(prog_path):
            sys.stderr.write('Expecting binary file but found directory: "%s"\n' % prog_path)
            sys.exit(1)        
    else:
        prog_path = which(prog_name)
        if prog_path is None:
            sys.stderr.write('Could not find "%s" on system path, and environmental variable "%s" specifying directory not set.\n' % (prog_name, GINGKO_BIN_PATH_ENVAR))
            sys.exit(1)
    return prog_path
    
def run_program(cmd, log=None):
    if log is not None:
        log.info('Invoking command: "%s"' % cmd)
    p1 = subprocess.Popen([cmd],
                           shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    stdout, stderr = p1.communicate()
    return stdout, stderr, p1.returncode

###############################################################################
## TESTS RUN    

def get_test_suite():
    """
    Creates a unittest.TestSuite from all of the modules in
    `dendropy.tests`. Right now, assumes (a) no subdirectories (though
    this can easily be accomodated) and (b) every test to be run is
    sitting in a module with a file name of 'test*.py', and, conversely,
    every file with a name of 'run_test_*.py' has test(s) to be run.
    """
    # get list of test file names'
    path = os.path.dirname(__file__)
    if not path:
        path = '.'       
    files = os.listdir(path)                               
    test_file_pattern = re.compile("run_test_.*\.py$", re.IGNORECASE)   
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
    