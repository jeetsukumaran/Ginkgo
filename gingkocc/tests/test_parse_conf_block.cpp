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

bool catch_block_parse_exception(const char * message, std::string s) {
    std::cerr << message << std::endl;
    bool exception_thrown = false;
    std::istringstream in0(s);
    try {
        ConfigurationBlockParser species_block(in0);
        exception_thrown = false;
    } catch (const ConfigurationParseError& e) {
        exception_thrown = true;
        std::cerr << "Exception correctly thrown: \"" << e.what() << "\"" << std::endl;
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
    

    catch_block_parse_exception("(testing empty block error)", "");
    
    catch_block_parse_exception("(testing missing block open error)", 
        species_block_str.substr(2, species_block_str.size()));
        
    catch_block_parse_exception("(testing missing block close error)", 
        species_block_str.substr(0, species_block_str.size()-2));

    std::string s1 = species_block_str;
    std::string::size_type s1_pos = s1.find("{");
    s1[s1_pos] = 'X';
    catch_block_parse_exception("(testing missing block open error)", s1);          
      
    s1.insert(s1_pos, "{{");    
    catch_block_parse_exception("(testing multiple block open error)", s1);
    
    s1 = species_block_str;
    s1_pos = s1.find(" Sp1");
    s1[s1_pos] = '_';
    catch_block_parse_exception("(testing single elemnt in head lack of label/name separation)", s1);
    
    s1 = species_block_str;
    s1_pos = s1.find(" Sp1");
    s1.insert(s1_pos, " extra");
    catch_block_parse_exception("(testing too many elements in head)", s1);    
    
    std::istringstream in1(species_block_str);
    ConfigurationBlockParser species_block(in1);   
}


