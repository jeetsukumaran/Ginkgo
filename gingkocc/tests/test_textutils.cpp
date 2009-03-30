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

#include "../textutils.hpp"

using namespace gingko::textutils;

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
    char * e4[] = {"a", "b", "**c&&&d%e f g h&i"};
    check_split( split_on_any("a+b***c&&&d%e f g h&i", "+*&%,!", 2, true), e4, 3);
    char * e5[] = {"a", "b", "", "*c&&&d%e f g h&i"};
    check_split( split_on_any("a+b***c&&&d%e f g h&i", "+*&%,!", 3, true), e5, 4);    
    
    char * a1[] = {"a", "b", "c", "d", "e f g h", "i"};
    check_split( split_on_any("a+b*c&d%e f g h&i", "+*&%,!", 0, false), a1, 6);
    char * a2[] = {"a", "b", "c", "d", "e f g h", "i"};
    check_split( split_on_any("a+b*c&&&d%e f g h&i", "+*&%,!", 0, false), a2, 6);
    char * a3[] = {"a", "b*c&&&d%e f g h&i"};
    check_split( split_on_any("a+b*c&&&d%e f g h&i", "+*&%,!", 1, false), a3, 2);
    char * a4[] = {"a", "b", "**c&&&d%e f g h&i"};
    check_split( split_on_any("a+b***c&&&d%e f g h&i", "+*&%,!", 2, false), a4, 3);
    char * a5[] = {"a", "b", "c", "&&d%e f g h&i"};
    check_split( split_on_any("a+b***c&&&d%e f g h&i", "+*&%,!", 3, false), a5, 4);     
}

void assert_str_equal(const std::string& result, const char * expected) {
    std::cout << "Expecting: \"" << expected << "\", received \"" << result << "\"" << std::endl;
    assert(result == expected);
}

void test_strip() {
    std::cout << "Testing string stripping ... " << std::endl;
    assert_str_equal(strip("    ", " "), "");
    assert_str_equal(strip("  a  ", " "), "a");
    assert_str_equal(strip("  a  ", " a"), "");
    assert_str_equal(strip("   hello    "), "hello");
    assert_str_equal(strip("   hello"), "hello");
    assert_str_equal(strip("hello"), "hello"); 
    assert_str_equal(strip("he llo"), "he llo");
    assert_str_equal(strip(" he llo  "), "he llo");
    assert_str_equal(strip("\n\th\tello\n\n\n \t\n"), "h\tello");
    assert_str_equal(strip("### hello ###", "#"), " hello ");
    assert_str_equal(strip("!### hello ###", "#"), "!### hello ");
    assert_str_equal(strip("!### hello ###  ", "!"), "### hello ###  ");
    assert_str_equal(strip("!### hello ###  ", "!#"), " hello ###  ");
    assert_str_equal(strip("!### hello ###  ", "!# "), "hello");    
}

void test_casing() {
    assert_str_equal(lower(" 12abcdABCD21 "), " 12abcdabcd21 ");
    assert_str_equal(upper(" 12abcdABCD21 "), " 12ABCDABCD21 "); 
}

int main(int, char * []) {
    test_splits();
    test_split_on_any();
    test_strip();
}


