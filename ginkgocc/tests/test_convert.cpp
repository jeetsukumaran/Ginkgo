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
#include <cassert>
#include "../convert.hpp"

template <typename T>
void check_scalar(T expected, T result) {
    std::cout << "Expecting: " << expected << ", Result: " << result << std::endl;
    assert(expected == result);
}

template <typename T>
void check_scalar_error(std::string src) {
    try {
        T x = ginkgo::convert::to_scalar<T>(src);
        std::cout << "Failed to throw value conversion error (\"" << src << "\"";
        std::cout << " => " << x << ")" << std::endl;
        assert(false);
    } catch (ginkgo::convert::ValueError& e) {
        std::cout << "Value error correctly thrown (\"" << src << "\")" << std::endl;
    }    
}

void check_vector_int(const std::vector<int>& result) {
    bool ok=true;
    std::cout << "Expecting: {1,2,3,4}, Result: {";
    for (std::vector<int>::const_iterator i = result.begin(); i != result.end(); ++i) {
        int k = *i;
        std::cout << k;
        if ( (i-result.begin()) < 3) {
            std::cout << ",";
        }
        if (k != ((i-result.begin())+1)) {
            ok = false;
        }
    }
    std::cout << "}" << std::endl;
    assert(ok);
}

void check_vector_float(const std::vector<float>& result) {
    bool ok=true;
    std::cout << "Expecting: {0.1,0.2,0.3,0.4}, Result: {";
    for (std::vector<float>::const_iterator i = result.begin(); i != result.end(); ++i) {
        float k = *i;
        std::cout << k;
        if ( (i-result.begin()) < 3) {
            std::cout << ",";
        }
        if (k != static_cast<float>(((i-result.begin())+1))/10) {
            ok = false;
        }
    }
    std::cout << "}" << std::endl;
    assert(ok);
}

void check_vector_string(const std::vector<std::string>& result, const char * expected[]) {
    bool ok=true;
    std::cout << "Result: {";
    for (std::vector<std::string>::const_iterator i = result.begin(); i != result.end(); ++i) {
        std::string k = *i;
        std::cout << k;
        if ( (i-result.begin()) < 3) {
            std::cout << " ";
        }
        const char * e = expected[(i-result.begin())];
        if (k != e) {
            ok = false;
        }
    }
    std::cout << "}" << std::endl;
    assert(ok);
}

int main(int, char *[]) {

    std::cout << "Testing scalar literal conversion ..." << std::endl;
    check_scalar(ginkgo::convert::to_scalar<int>("1"), 1);
    check_scalar(ginkgo::convert::to_scalar<int>("100"), 100); 
    check_scalar(ginkgo::convert::to_scalar<int>("-1"), -1);
    check_scalar(ginkgo::convert::to_scalar<double>("0.5"), 0.5);
    check_scalar(ginkgo::convert::to_scalar<double>("-0.5"), -0.5);
    std::cout << "Scalar literal conversion: PASS" << std::endl;

    std::cout << "Testing scalar literal conversion errors ... " << std::endl;
    check_scalar_error<int>("abc");
    check_scalar_error<int>("abc1");
    check_scalar_error<int>("1abc");
    check_scalar_error<int>("a2b");
    check_scalar_error<int>("2 2");
    check_scalar_error<unsigned>("-1");
    std::cout << "Scalar literal conversion errors: PASS" << std::endl;    
    
    std::cout << "Testing vector literal conversion ..." << std::endl;
    check_vector_int(ginkgo::convert::to_vector<int>("1 2 3   4"));
    check_vector_float(ginkgo::convert::to_vector<float>("0.1   0.2 0.3 0.4"));
    const char * e1[] = {"a", "b", "c", "d"};
    check_vector_string(ginkgo::convert::to_vector<std::string>("a   b c d"), e1);
    const char * e2[] = {"1,1", "2,2", "3,3", "4,4"};
    check_vector_string(ginkgo::convert::to_vector<std::string>("1,1  2,2 3,3 4,4"), e2);    
    std::cout << "Vector literal conversion: PASS" << std::endl;
}
