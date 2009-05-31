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

#include "../tree.hpp"
#include "../textutil.hpp"
#include <iostream>
#include <vector>
#include <cstdlib>

int main(int argc, char* argv[]) {

    if (argc == 1) {
        std::cerr << "No parent indexes specified" << std::endl;
        exit(1);
    }
    
    gingko::Tree tree(NULL, true);
    const char * label;
    for (int i = 1; i < argc; ++i) {
        std::vector< std::string > parts = gingko::textutil::split(argv[i], ":");
        unsigned long idx = atol(parts[0].c_str());            
        if (parts.size() > 1) {
            label = parts[1].c_str();
        } else {
            label = NULL;
        }
        tree.add_indexed_node(idx, label);
    }
    tree.dump(std::cerr);
    std::cerr << std::endl;
    tree.write_newick_tree(std::cout);
    std::cerr << std::endl;
}
