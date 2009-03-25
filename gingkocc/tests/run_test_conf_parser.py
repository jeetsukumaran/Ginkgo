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

from gingko_tests import get_logger
from gingko_tests import get_gingko_program_path
from gingko_tests import run_program

_LOG = get_logger("test_conf_parser")

class ParseConfBlockTest(unittest.TestCase):
    def setUp(self):
        self.prog_path = get_gingko_program_path("test_parse_conf_block")
            
if __name__ == "__main__":
    unittest.main()
    