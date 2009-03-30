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

#include "../cmdopt.hpp"
#include <iostream>
#include <string>

using namespace gingko;

int main(int argc, char* argv[]) {
    long          a = 1000;
    unsigned long b = 1000;
    int           c = 1000;
    unsigned int  d = 1000;
    double        e = 0.10;
    std::string   f = "default 1";
    bool          g = false;
    

    gingko::OptionParser parser = gingko::OptionParser();
    parser.add_option<long>(&a, "-a", "--seta", 
                                     "set value of a", "#");
    parser.add_option<unsigned long>(&b, "-b", "--setb", 
                                     "set value of a", "#");
    parser.add_option<int>(&c, "-c", "--setc", 
                                     "set value of a", "#");
    parser.add_option<unsigned int>(&d, "-d", "--setd", 
                                     "set value of a", "#");
    parser.add_option<double>(&e, "-e", "--sete", 
                                     "set value of e", "#");
    parser.add_option<std::string>(&f, "-f", "--setf", 
                                     "set value of a", "#");
    parser.add_switch(&g, "-g", "--setg", "set value of e", "#");
    parser.parse(argc, argv);
    
    std::cerr << "--Flagged options:" << std::endl;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << c << std::endl;
    std::cout << d << std::endl;    
    std::cout << e << std::endl;
    std::cout << f << std::endl;    
    std::cout << g << std::endl;
    std::cerr << "--Position arguments:" << std::endl;
    OptionParser::PosArgs args = parser.get_args();
    for (OptionParser::PosArgs::iterator a = args.begin(); a != args.end(); ++a) {
        std::cout << *a << std::endl;
    }
    std::cerr << "(done)" << std::endl;
    
}
