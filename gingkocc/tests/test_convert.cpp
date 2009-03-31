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
#include <cassert>
#include "../convert.hpp"

template <typename T>
void check_val(T expected, T result) {
    std::cout << "Expecting: " << expected << ", Result: " << result << std::endl;
    assert(expected == result);
}

template <typename T>
void check_error(std::string src) {
    try {
        T x = gingko::convert::to_scalar<T>(src);
        std::cout << "failed to throw value conversion error (\"" << src << "\"";
        std::cout << " => " << x << ")" << std::endl;
        assert(false);
    } catch (gingko::convert::ValueError& e) {
        std::cout << "value error correctly thrown (\"" << src << "\")" << std::endl;
    }    
}


int main(int, char *[]) {

    check_val(gingko::convert::to_scalar<int>("1"), 1);
    check_val(gingko::convert::to_scalar<int>("100"), 100); 
    check_val(gingko::convert::to_scalar<int>("-1"), -1);
    check_val(gingko::convert::to_scalar<double>("0.5"), 0.5);
    check_val(gingko::convert::to_scalar<double>("-0.5"), -0.5);

    check_error<int>("abc");
    check_error<int>("abc1");
    check_error<int>("1abc");
    check_error<int>("a2b");
    check_error<int>("2 2");
    check_error<unsigned>("-1");
    
}
