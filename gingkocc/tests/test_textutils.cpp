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

#include "../textutils.h"

using namespace gingko;

void test_extract_filename_from_path() {
    char * input[] = {"/home/user/filename", "user/filename", "filename", ""};
    char ** check = input;
    std::cout << "Testing \"extract_filename_from_path\" ..." << std::endl;
    while (strcmp(*check, "") != 0) {
        std::cout << "\"" << *check << "\" => ";
        std::string result = extract_filename_from_path(*check);
        std::cout << "\"" << result << "\"" << std::endl;
        assert(strcmp(result.c_str(), "filename") == 0);
        ++check;
    }
    std::cout << "\"extract_filename_from_path\" tests: PASS" << std::endl;
}

bool check_split(std::vector<std::string> result,
                 char * expected[],
                 unsigned expected_count) {
    std::cout << "Testing string splitting ..." << std::endl;
    std::cout << "Expecting size: " << expected_count << std::endl;
    std::cout << "Resulting size: " << result.size() << std::endl;
    assert(result.size() == expected_count);
    std::cout << "Expecting values:";
    for (unsigned i = 0; i < expected_count; ++i) {
        std::cout << " " << expected[i];
    }
    std::cout << std::endl;
    std::cout << "Resulting values:";
    for (unsigned i = 0; i < expected_count; ++i) {
        std::cout << " " << result[i];
    }
    std::cout << std::endl;    
    for (unsigned i = 0; i < expected_count; ++i) {
        assert(result[i] == expected[i]);
    }
    
    return true;
}

int main(int, char * []) {
    test_extract_filename_from_path();
    
    char * e1[] = {"the", "quick", "brown", "fox", "jumps", "over", "the", "lazy", "dog"};
//     split("the quick brown fox jumps over the lazy dog", " ", 0, true);
    check_split( split("the quick brown fox jumps over the lazy dog", " ", 0, true),
                e1, 9);
}


