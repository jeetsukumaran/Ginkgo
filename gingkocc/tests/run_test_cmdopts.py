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
from run_tests import run_program

_LOG = get_logger("test_cmdopts")

class CmdOptsTest(unittest.TestCase):
    def setUp(self):
        self.prog_path = get_gingko_program_path("test_cmdopts")
        
    def check_opts_parsing(self, cmd, expected_strings):  
        _LOG.info('Testing arguments parsing: "%s"' % cmd)
        stdout, stderr, returncode = run_program(cmd)
        assert returncode == 0, "Program exited with error code %d:\n%s" % (p1.returncode, stderr)
        stdout = stdout.split("\n")        
        for idx, expected in enumerate(expected_strings):
            _LOG.info('Line %d: %s (correct = "%s")' % (idx+1, stdout[idx], expected))        
            assert stdout[idx] == expected, \
                'Expecting "%s", but found "%s" in line %d of stdout' \
                % (expected, stdout[idx], idx+1)
                
    def check_wrong_opts_parsing(self, wrong_flag):  
        cmd = self.prog_path + " " + wrong_flag
        _LOG.info('Testing for correct rejection of wrong or invalid options: %s' % cmd)
        stdout, stderr, returncode = run_program(cmd)
        assert returncode != 0, "Program did not exit with error on wrong flag:\n%s\n%s" % (stdout, stderr)
        
    def testDefaultArgs(self):
        self.check_opts_parsing(self.prog_path, ["1000", "1000", "1000", "1000", "0.1", "default 1", "0"])
        
    def testShortFlagsArgs(self):
        cmd = self.prog_path + ' -a -10 -b 10 -c -100 -d 100 -e 2.718 -f "the quick brown fox" -g'
        self.check_opts_parsing(cmd, ["-10", "10", "-100", "100", "2.718", "the quick brown fox", "1"])        
        
    def testLongFlagsArgs(self):
        cmd = self.prog_path + ' --seta -10 --setb 10 --setc -100 --setd 100 --sete 2.718 --setf "the quick brown fox" --setg'
        self.check_opts_parsing(cmd, ["-10", "10", "-100", "100", "2.718", "the quick brown fox", "1"]) 
        
    def testWrongFlag(self):
        self.check_wrong_opts_parsing("-x")
        self.check_wrong_opts_parsing("--hello")
       
    def testInvalidValue(self):
        _LOG.info('Testing for correct rejection of invalid values')
        stdout, stderr, returncode = run_program(self.prog_path + " -a xx", _LOG)
        assert returncode != 0
        stdout, stderr, returncode = run_program(self.prog_path + " -a 1x", _LOG)
        assert returncode != 0
        stdout, stderr, returncode = run_program(self.prog_path + " -a x1", _LOG)
        assert returncode != 0
        stdout, stderr, returncode = run_program(self.prog_path + " -a 1.1", _LOG)
        assert returncode != 0
        stdout, stderr, returncode = run_program(self.prog_path + " --seta xx", _LOG)
        assert returncode != 0
        stdout, stderr, returncode = run_program(self.prog_path + " --seta 1x", _LOG)
        assert returncode != 0
        stdout, stderr, returncode = run_program(self.prog_path + " --seta x1", _LOG)
        assert returncode != 0
        stdout, stderr, returncode = run_program(self.prog_path + " --seta 1.1", _LOG)
        assert returncode != 0          

if __name__ == "__main__":
    unittest.main()
    