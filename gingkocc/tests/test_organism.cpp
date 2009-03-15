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

#include <ctime>
#include <iostream>
#include <vector>
#include "../biosys.h"
#include "../randgen.h"
#include "../tree.h"

using namespace gingko;

int main(int argc, char* argv[]) {
    long seed = 0;
    if (argc > 1) {
        seed = atol(argv[1]);
    } else {
        seed = time(0);
    }
    
    std::cout << "\nSeed = " << seed << "\n";
    RandomNumberGenerator rng(seed);
    Species sp(0, "gecko", 4, rng);
    
    Organism g0_1 = sp.new_organism();

    Organism g1_1 = sp.new_organism(g0_1, g0_1);
    Organism g1_2 = sp.new_organism(g0_1, g0_1);

    Organism g2_1 = sp.new_organism(g1_1, g1_1);
    Organism g2_2 = sp.new_organism(g1_2, g1_2);
    Organism g2_3 = sp.new_organism(g1_2, g1_2);
    
    Organism g3_1 = sp.new_organism(g2_1, g2_1);
    Organism g3_2 = sp.new_organism(g2_2, g2_2);
    Organism g3_3 = sp.new_organism(g2_2, g2_2);
    Organism g3_4 = sp.new_organism(g2_3, g2_3);
    Organism g3_5 = sp.new_organism(g2_3, g2_3);      
    
    Tree tree;
    tree.process_node(g3_1.haploid_marker().node());
    tree.process_node(g3_2.haploid_marker().node());
    tree.process_node(g3_3.haploid_marker().node());
    tree.process_node(g3_4.haploid_marker().node());
    tree.process_node(g3_5.haploid_marker().node());    
    
    tree.dump(std::cerr);
    std::cerr << "\n---\n";
    tree.write_newick_tree(std::cout);
    std::cerr << std::endl;    
     
}