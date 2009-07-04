#! /usr/bin/env python

###############################################################################
##
## GINKGO Biogeographical Evolution Simulator.
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
from ginkgo_tests import get_logger
from ginkgo_tests import get_ginkgo_program_path
from ginkgo_tests import run_program
from ginkgo_tests import run_external_tests

_LOG = get_logger("test_asciigrid")

class TextUtilsTest(unittest.TestCase):

    def setUp(self):
        self.prog_path = get_ginkgo_program_path("test_asciigrid")

    def testTextUtils(self):
        run_external_tests(self.prog_path, _LOG, "ASCII Grid parsing tests")       
            
if __name__ == "__main__":
    unittest.main()
    