///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Phylogeographical Evolution Simulator.
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

#include "../filesys.hpp"

using namespace ginkgo::filesys;

void check_path_op( std::string (*path_op_func)(const char *),  const char * path, const char * expected) {
    std::string result = path_op_func(path);
    std::cout << "-- \"" << path << "\" => \"" << result << "\"";
    std::cout << " (expecting \"" << expected << "\")" << std::endl;
    assert(result == expected);
}

void test_get_path_leaf() {
    std::cout << "Testing \"get_path_leaf\" ..." << std::endl;
    check_path_op(get_path_leaf, "/home/user/filename", "filename");
    check_path_op(get_path_leaf, "user/filename", "filename");
    check_path_op(get_path_leaf, "filename", "filename");
    check_path_op(get_path_leaf, "/home/user/dir/", "dir");
    check_path_op(get_path_leaf, "/home/user/filename", "filename");
    check_path_op(get_path_leaf, "", "");
    std::cout << "\"get_path_leaf\" tests: PASS" << std::endl;
}

void test_get_path_parent() {
    std::cout << "Testing \"test_get_path_parent\" ..." << std::endl;
    check_path_op(get_path_parent, "/home/user/filename", "/home/user");
    check_path_op(get_path_parent, "/home/user/filename/", "/home/user");    
    check_path_op(get_path_parent, "user/filename", "user");
    check_path_op(get_path_parent, "user/dir/filename", "user/dir");    
    check_path_op(get_path_parent, "filename", "");    
    std::cout << "\"test_get_path_parent\" tests: PASS" << std::endl;
}

int main(int, char * []) {
    test_get_path_leaf();
    test_get_path_parent();
}
