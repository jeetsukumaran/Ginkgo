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

#include "../tree.h"
#include "../cmdopts.h"
#include "../textutils.h"
#include <iostream>
#include <vector>
#include <cstdlib>

int main(int argc, char* argv[]) {

    gingko::OptionParser parser = gingko::OptionParser("Tree Testing",
            "Given a list of parent indexes (with optional labels), returns a NEWICK representation of the tree", 
            "%prog [options] <PARENT INDEX>[:<LABEL>] <PARENT INDEX>[:<LABEL>] [<PARENT INDEX>[:<LABEL>] ... ");

    parser.parse(argc, argv);
    std::vector< std::string > args = parser.get_args();
    if (args.size() == 0) {
        std::cerr << "No parent indexes specified" << std::endl;
        exit(1);
    }
    
    gingko::Tree tree(true);
    const char * label;
    for (std::vector< std::string >::iterator arg_iter = args.begin();
         arg_iter != args.end();
         ++arg_iter) {
        std::vector< std::string > parts = gingko::split(*arg_iter, ":");
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
