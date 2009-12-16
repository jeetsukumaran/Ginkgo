///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Biogeographical Evolution Simulator.
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

#include "../asciigrid.hpp"

using namespace ginkgo;

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

bool check_against_expected_grid1(asciigrid::AsciiGrid<long>& ag) {
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

bool check_against_expected_grid2(asciigrid::AsciiGrid<long>& ag) {
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

const char * TEST_GRID3 =
        "ncols         10\n"
        "nrows         6\n"
        "xllcorner     0.0\n"
        "yllcorner     0.0\n"
        "cellsize      50.0\n"
        "0.1 1.2 2.3 3.4 4.5 5.6 6.7 7.8 8.9 9.1\n"
        ".10 .11 .12 .13 .14 .15 .16 .17 .18 .19\n"
        "0.20 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29\n"
        "3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9\n"
        "40 41 42 43 44 45 46 47 48 49\n"
        "50 51 52 53 54 55 56 57 58 59\n";

bool check_against_expected_grid3(asciigrid::AsciiGrid<float>& ag) {
    assert(ag.get_ncols() == 10);
    assert(ag.get_nrows() == 6);
    assert(ag.has_size(10,6));
    std::cout << "Grid 3: reported size OK" << std::endl;
    std::vector<float> cell_vals = ag.get_cell_values();
    assert(cell_vals.size() == 60);
    std::cout << "Grid 3: actual size OK" << std::endl;
    float expected[] = { 0.1,  1.2,  2.3,  3.4,  4.5,  5.6,  6.7,  7.8,  8.9,  9.1,
                       .10, .11, .12, .13, .14, .15, .16, .17, .18, .19,
                       0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
                       3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
                       40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                       50, 51, 52, 53, 54, 55, 56, 57, 58, 59,};
    for (int i=0; i<24; ++i) {
        std::cout << cell_vals[i] << " <-> " << expected[i] << std::endl;
        assert(cell_vals[i] == expected[i]);
    }
    std::cout << "Grid 2: values OK" << std::endl;
    return true;
}

void run_internal_tests() {
    std::cout << "Testing Grid 1 (from string)" << std::endl;
    std::istringstream g1(TEST_GRID1);
    asciigrid::AsciiGrid<long> ag1(g1);
    assert(check_against_expected_grid1(ag1));
    std::cout << "Grid 1 (from string): PASS" << std::endl;

    std::cout << "Testing Grid 1 (round-trip)" << std::endl;
    std::vector<long> vals1 = ag1.get_cell_values();
    std::ostringstream sa1;
    asciigrid::write_grid(vals1, ag1.get_ncols(), ag1.get_nrows(), std::cout);
    asciigrid::write_grid(vals1, ag1.get_ncols(), ag1.get_nrows(), sa1);
    std::istringstream sa1i(sa1.str());
    asciigrid::AsciiGrid<long> ag1b(sa1i);
    assert(check_against_expected_grid1(ag1b));
    std::cout << "Grid 1 (round-trip): PASS" << std::endl;

    std::cout << "Testing Grid 2 (from string)" << std::endl;
    std::istringstream g2(TEST_GRID2);
    asciigrid::AsciiGrid<long> ag2(g2);
    assert(check_against_expected_grid2(ag2));
    std::cout << "Grid 2 (from string): PASS" << std::endl;

    std::cout << "Testing Grid 2 (round-trip)" << std::endl;
    std::vector<long> vals2 = ag2.get_cell_values();
    std::ostringstream sa2;
    asciigrid::write_grid(vals2, ag2.get_ncols(), ag2.get_nrows(), std::cout);
    asciigrid::write_grid(vals2, ag2.get_ncols(), ag2.get_nrows(), sa2);
    std::istringstream sa2i(sa2.str());
    asciigrid::AsciiGrid<long> ag2b(sa2i);
    assert(check_against_expected_grid2(ag2b));
    std::cout << "Grid 2 (round-trip): PASS" << std::endl;

    std::cout << "Testing Grid 3 (from string)" << std::endl;
    std::istringstream g3(TEST_GRID3);
    asciigrid::AsciiGrid<float> ag3(g3);
    assert(check_against_expected_grid3(ag3));
    std::cout << "Grid 3 (from string): PASS" << std::endl;

    std::cout << "Testing Grid 3 (round-trip)" << std::endl;
    std::vector<float> vals3 = ag3.get_cell_values();
    std::ostringstream sa3;
    asciigrid::write_grid(vals3, ag3.get_ncols(), ag3.get_nrows(), std::cout);
    asciigrid::write_grid(vals3, ag3.get_ncols(), ag3.get_nrows(), sa3);
    std::istringstream sa3i(sa3.str());
    asciigrid::AsciiGrid<float> ag3b(sa3i);
    assert(check_against_expected_grid3(ag3b));
    std::cout << "Grid 3 (round-trip): PASS" << std::endl;
}

int main(int, char *) {
    run_internal_tests();
//    if (argc == 1) {
//        run_internal_tests();
//    } else {
//        asciigrid::AsciiGrid ag(argv[1]);
//        std::vector<long> v = ag.get_cell_values();
//        for (std::vector<long>::const_iterator i = v.begin(); i != v.end(); ++i) {
//            std::cout << *i << std::endl;
//        }
//    }
}

