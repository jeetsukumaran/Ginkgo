///////////////////////////////////////////////////////////////////////////////
//
// GINGKO Biogeographical Evolution Simulator.
//
// Copyright 2009 Jeet Sukumaran and Mark T. Holder.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License along
// with this program. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <cassert>

#include "../asciigrid.h"

using namespace gingko;

const char * TEST_GRID1 =
        "ncols         4\n"
        "nrows         6\n"
        "xllcorner     0.0\n"
        "yllcorner     0.0\n"
        "cellsize      50.0\n"
        "NODATA_value  -9999\n"
        "-9999 -9999 5 2\n"
        "-9999 20 100 36\n"
        "3 8 35 10\n"
        "32 42 50 6\n"
        "88 75 27 9\n"
        "13 5 1 -9999\n";

const char * TEST_GRID2 =
        "ncols         10\n"
        "nrows         6\n"
        "xllcorner     0.0\n"
        "yllcorner     0.0\n"
        "cellsize      50.0\n"
        "0 1 2 3 4 5 6 7 8 9\n"
        "10 11 12 13 14 15 16 17 18 19\n"
        "20 21 22 23 24 25 26 27 28 29\n"
        "30 31 32 33 34 35 36 37 38 39\n"
        "40 41 42 43 44 45 46 47 48 49\n"
        "50 51 52 53 54 55 56 57 58 59\n";
        
bool check_against_expected_grid1(asciigrid::AsciiGrid& ag) {
    assert(ag.get_ncols() == 4);
    assert(ag.get_nrows() == 6);
    assert(ag.has_size(4,6));
    std::cout << "Grid 1: reported size OK" << std::endl;
    std::vector<long> cell_vals = ag.get_cell_values();
    assert(cell_vals.size() == 24);
    std::cout << "Grid 1: actual size OK" << std::endl;    
    long expected[] = {-9999, -9999,   5,     2, 
                       -9999,    20, 100,    36, 
                           3,     8,  35,    10,
                          32,    42,  50,     6, 
                          88,    75,  27,     9,
                          13,     5,   1, -9999};
    for (int i=0; i<24; ++i) {
        assert(cell_vals[i] == expected[i]);
    }
    std::cout << "Grid 1: values OK" << std::endl;
    return true;
}
 
bool check_against_expected_grid2(asciigrid::AsciiGrid& ag) {
    assert(ag.get_ncols() == 10);
    assert(ag.get_nrows() == 6);
    assert(ag.has_size(10,6));
    std::cout << "Grid 2: reported size OK" << std::endl;    
    std::vector<long> cell_vals = ag.get_cell_values();
    assert(cell_vals.size() == 60);
    std::cout << "Grid 2: actual size OK" << std::endl;    
    long expected[] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 
                       10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
                       20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
                       30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
                       40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 
                       50, 51, 52, 53, 54, 55, 56, 57, 58, 59,};
    for (int i=0; i<24; ++i) {
        assert(cell_vals[i] == expected[i]);
    }
    std::cout << "Grid 2: values OK" << std::endl;    
    return true;
}
 
void run_internal_tests() {
    std::cout << "Testing Grid 1 (from string)" << std::endl;
    std::istringstream g1(TEST_GRID1);
    asciigrid::AsciiGrid ag1(g1);
    assert(check_against_expected_grid1(ag1));
    std::cout << "Grid 1 (from string): PASS" << std::endl;
    
    std::cout << "Testing Grid 2 (from string)" << std::endl;    
    std::istringstream g2(TEST_GRID2);
    asciigrid::AsciiGrid ag2(g2);
    assert(check_against_expected_grid2(ag2));
    std::cout << "Grid 2 (from string): PASS" << std::endl;
}

int main(int argc, char * argv[]) {
    if (argc == 1) {
        run_internal_tests();
    } else {       
        asciigrid::AsciiGrid ag(argv[1]);
        std::vector<long> v = ag.get_cell_values();
        for (std::vector<long>::const_iterator i = v.begin(); i != v.end(); ++i) {
            std::cout << *i << std::endl;
        }
    }   
}

