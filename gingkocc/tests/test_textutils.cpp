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
    if (result.size() != expected_count) {
        std::cerr << " *** INCORRECT RESULT SIZE: ";
        for (std::vector<std::string>::iterator r = result.begin(); r != result.end(); ++r) {
            std::cerr << " [" << *r << "] ";
        }
        std::cerr << std::endl;
    }    
    assert(result.size() == expected_count);
    std::cout << "Expecting values: ";
    for (unsigned i = 0; i < expected_count; ++i) {
        std::cout << " [" << expected[i] << "] ";
    }
    std::cout << std::endl;
    std::cout << "Resulting values: ";
    for (unsigned i = 0; i < expected_count; ++i) {
        std::cout << " [" << result[i] << "] ";
    }
    std::cout << std::endl;    
    for (unsigned i = 0; i < expected_count; ++i) {
        assert(result[i] == expected[i]);
    }
    
    return true;
}

void test_splits() {
    char * e1[] = {"the", "quick", "brown", "fox", "jumps", "over", "the", "lazy", "dog"};
    check_split( split("the quick brown fox jumps over the lazy dog", " ", 0, true),
                e1, 9);
       
    char * e2[] = {"the", "quick brown fox jumps over the lazy dog"};
    check_split( split("the quick brown fox jumps over the lazy dog", " ", 1, true),
                e2, 2);    
                
    char * e3[] = {"the", "quick", "brown fox jumps over the lazy dog"};
    check_split( split("the quick brown fox jumps over the lazy dog", " ", 2, true),
                e3, 3);
                
    char * e4 [] = {"", "", "a", "", "b", "", "", "c", "", ""};
    check_split( split(",,a,,b,,,c,,", ",", 0, true), e4, 10 );
    
    char * e5 [] = {"a", "b", "c"};
    check_split( split(",,a,,b,,,c,,", ",", 0, false), e5, 3 );
    
    char * e6 [] = {"", ",a,,b,,,c,,"};
    check_split( split(",,a,,b,,,c,,", ",", 1, true), e6, 2 );
    
    char * e7 [] = {"a", ",b,,,c,,"};
    check_split( split(",,a,,b,,,c,,", ",", 1, false), e7, 2 );
        
    char * e8 [] = {"", "", "a,,b,,,c,,"};
    check_split( split(",,a,,b,,,c,,", ",", 2, true), e8, 3 );
    
    char * e9 [] = {"a", "b", ",,c,,"};
    check_split( split(",,a,,b,,,c,,", ",", 2, false), e9, 3 );    
}

void test_split_on_any() {
    char * e1[] = {"a", "b", "c", "d", "e f g h", "i"};
    check_split( split_on_any("a+b*c&d%e f g h&i", "+*&%,!", 0, true), e1, 6);
    char * e2[] = {"a", "b", "c", "", "", "d", "e f g h", "i"};
    check_split( split_on_any("a+b*c&&&d%e f g h&i", "+*&%,!", 0, true), e2, 8);
    char * e3[] = {"a", "b*c&&&d%e f g h&i"};
    check_split( split_on_any("a+b*c&&&d%e f g h&i", "+*&%,!", 1, true), e3, 2);    
}

int main(int, char * []) {
    test_extract_filename_from_path();
    test_splits();
    test_split_on_any();
}


