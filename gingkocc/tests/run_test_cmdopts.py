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

import unittest
import subprocess

from run_tests import get_logger
from run_tests import get_gingko_program_path

_LOG = get_logger("test_cmdopts")

class SumTreesTest(unittest.TestCase):
    def setUp(self):
        self.prog_path = get_gingko_program_path("test_cmdopts")
        
    def testDefaultArgs(self):
        cmd = self.prog_path    
        _LOG.info('Testing default arguments using command: "%s"' % cmd)
        p1 = subprocess.Popen([cmd],
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        stdout, stderr = p1.communicate()
        stdout = stdout.split("\n")
        assert p1.returncode == 0, "Program exited with error code %d" % p1.returncode
        expected_strings = ["1000", "1000", "1000", "1000", "0.1", "default 1", "0"]
        for idx, expected in enumerate(expected_strings):
            _LOG.info('Line %d: %s (correct: "%s")' % (idx+1, stdout[idx], expected))        
            assert stdout[idx] == expected, \
                'Expecting "%s", but found "%s" in line %d of stdout' \
                % (expected, stdout[idx], idx+1)            
            

if __name__ == "__main__":
    unittest.main()
    