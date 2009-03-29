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

#include "../filesys.h"

using namespace gingko::filesys;

void test_get_filename() {
    char * input[] = {"/home/user/filename", "user/filename", "filename", ""};
    char ** check = input;
    std::cout << "Testing \"get_filename\" ..." << std::endl;
    while (strcmp(*check, "") != 0) {
        std::cout << "\"" << *check << "\" => ";
        std::string result = get_filename(*check);
        std::cout << "\"" << result << "\"" << std::endl;
        assert(strcmp(result.c_str(), "filename") == 0);
        ++check;
    }
    std::cout << "\"get_filename\" tests: PASS" << std::endl;
}

int main(int, char * []) {
    test_get_filename();
}
