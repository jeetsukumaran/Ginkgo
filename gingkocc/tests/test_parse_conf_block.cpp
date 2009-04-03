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
#include <cassert>

#include "../confsys.hpp"

using namespace gingko::confsys;

const char * SPECIES_BLOCK_CSTR = "\n\n@species Sp1 { \n"
                                    "   selection = 1 1 1 1 ; \n"
                                    "   genotype =  0 0 0 0 ;  \n" 
                                    "   mutation-rate = 0.01 ;\n"
                                    "   max-mutation-size = 1 \n"
                                    "   fecundity =    8 ; \n"
                                    "   fecundity-evolution-size = 2;\n"
                                    "   movement-capacity = 10;\n"
                                    "   movement-costs = /Users/user/data/grid.asc \n"
                                    "}\n\n";

bool catch_block_parse_exception(const char * message, std::string s) {
    std::cout << message << std::endl;
    bool is_exception = false;
    std::istringstream in0(s);
    try {
        ConfigurationBlock species_block(in0);
        is_exception = false;
    } catch (const ConfigurationSyntaxError& e) {
        is_exception = true;
        std::cout << "Exception correctly thrown: \"" << e.what() << "\"" << std::endl;
    }
    assert(is_exception);
    return is_exception;
}

void test_parse_errors() {
    std::cout << "Testing configuration block parse error detection ..." << std::endl;
    std::string species_block_str(SPECIES_BLOCK_CSTR);    
//     catch_block_parse_exception("(testing empty block error)", "");    
    catch_block_parse_exception("(testing missing block open error)", 
        species_block_str.substr(3, species_block_str.size()));        
    catch_block_parse_exception("(testing missing block close error)", 
        species_block_str.substr(0, species_block_str.size()-3));
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
    s1 = species_block_str;
    s1_pos = s1.find("=");
    s1[s1_pos] = 'X';
    catch_block_parse_exception("(testing incomplete key-val entry)", s1);      
}    

void test_parse_correct() {
    std::cout << "Testing configuration block parse ..." << std::endl;
    std::string species_block_str(SPECIES_BLOCK_CSTR);
    std::istringstream in0(species_block_str);
    ConfigurationBlock species_block = ConfigurationBlock(in0);
    assert(species_block.get_entry<std::string>("selection") == "1 1 1 1");
    assert(species_block.get_entry<std::string>("genotype") == "0 0 0 0");
    assert(species_block.get_entry<float>("mutation-rate") == 0.01);
    assert(species_block.get_entry<int>("max-mutation-size") == 1);
    assert(species_block.get_entry<unsigned>("fecundity") == 8);
    assert(species_block.get_entry<long>("movement-capacity") == 10);    
    assert(species_block.get_entry<std::string>("movement-costs") == "/Users/user/data/grid.asc");   
    assert(species_block.get_keys().size() == 8);
}

void test_parse_empty() {
    std::cout << "Testing empty configuration block parse ..." << std::endl;
    std::string s = "";
    std::istringstream in0(s);
    ConfigurationBlock cb(in0);
    assert(cb.is_block_set() == false);
}

int main(int, char * []) {
    test_parse_correct();
    test_parse_empty();
    test_parse_errors();
    std::cout << "\nPASS\n" << std::endl;
}


