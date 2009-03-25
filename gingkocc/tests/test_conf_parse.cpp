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
#include <fstream>
#include <istream>
#include <string>
#include <sstream>

#include "../conf.h"

using namespace gingko;

// -- tests -- 
bool test_empty_block() {
    std::cerr << "(testing empty block)" << std::endl;
    bool exception_thrown = false;
    std::istringstream in0("");
    try {
        ConfigurationBlockParser species_block(in0);
        exception_thrown = false;
    } catch (const ConfigurationParseError& e) {
        exception_thrown = true;    
    }
    assert(exception_thrown);
    return exception_thrown;
}

bool test_missing_block_open(std::string s) {
    std::cerr << "(testing missing block open)" << std::endl;
    bool exception_thrown = false;
    std::istringstream in0(s);
    try {
        ConfigurationBlockParser species_block(in0);
        exception_thrown = false;
    } catch (const ConfigurationParseError& e) {
        exception_thrown = true;    
    }
    assert(exception_thrown);
    return exception_thrown;
}

bool test_missing_block_close(std::string s) {
    std::cerr << "(testing missing block close)" << std::endl;
    bool exception_thrown = false;                                        
    std::istringstream in0(s);   
    try {
        ConfigurationBlockParser species_block(in0);
        exception_thrown = false;
    } catch (const ConfigurationParseError& e) {
        exception_thrown = true;    
    }
    assert(exception_thrown);
    return exception_thrown;
}

int main(int, char * []) {
//     const char * src = "setup.conf";
//     std::ifstream in(src);
//     std::string s;
//     while (in) {
//         std::getline(in, s);
//         std::cout << "'" << s << "'" << std::endl;
//     }
    
    const char * species_block_cstr = " @species Sp1 { \n"
                                        "   selection = (1,1,1,1); \n"
                                        "   genotype = (0,0,0,0);\n" 
                                        "   mutation-rate = 0.01;\n"
                                        "   max-mutation-size = 1;\n"
                                        "   fecundity = 8;\n"
                                        "   fecundity-evolution-size = 2;\n"
                                        "   movement-capacity = 10;\n"
                                        "}\n";
    std::string species_block_str(species_block_cstr);
    
    // confirm that exception will be thrown if block is empty
    test_empty_block();    

    // confirm that exception will be thrown if block does not start with 
    // correct leading character ("@")
    test_missing_block_open(species_block_str.substr(2, species_block_str.size()));
    
    // confirm that exception will be thrown if block does not start with 
    // correct block close ("}")
    test_missing_block_close(species_block_str.substr(0, species_block_str.size()-2));
    
    
    std::istringstream in1(species_block_str);
    ConfigurationBlockParser species_block(in1);   
}


