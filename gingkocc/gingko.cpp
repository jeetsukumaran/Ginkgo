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

#include "gingko_defs.h"
#include "biosys.h"
#include "world.h"
#include "tree.h"
#include "cmdopts.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <set>

int main(int argc, char* argv[]) {
    gingko::OptionParser parser = gingko::OptionParser();
    parser.add_option<std::string>("-x", "--dim-x", "size of landscape in the x-dimension", "X");
    parser.add_option<std::string>("-y", "--dim-y", "size of landscape in the y-dimension", "Y");
    parser.add_option<std::string>("-c", "--carrying-capacity", "maximum carrying capacity of each cell", "K");
    parser.add_option<std::string>("-n", "--num-gens","number of generations to run", "#GENERATIONS");  
    parser.add_option<std::string>("-f", "--num-fitness", "number of fitness factors", "#FACTORS");
    parser.add_option<std::string>("-z", "--random-seed", "random number seed", "SEED");                      


    parser.parse(argc, argv);
       

//     long seed = 0;
// 	std::vector<unsigned long> args_as_longs;
// 	for (int i = 1; i < argc; ++i) {
// 		std::string a(argv[i]);
// 		if (a.length() > 2 && a[0] == '-') {
// 			if (a[1] == 's') {
// 				seed = atol(a.c_str() + 2);
// 			}
// 			else {
// 				std::cerr << "Unknown flag: " << a << '\n';
// 				exit(2);
// 			}
// 		}
// 		else {
// 			args_as_longs.push_back(atol(a.c_str()));
// 		}
// 	}
//     if (args_as_longs.size() < 5) {
//         std::cerr << "usage: " << argv[0] <<  " <DIM-X> <DIM-Y> <CELL-CARRYING-CAPACITY> <NUM-CELLS-TO-POPULATE> <NUM-GENS>\n";
//         exit(1);
//     }

//     long size_x = args_as_longs[0];
//     long size_y = args_as_longs[1];
//     long cc = args_as_longs[2];
//     long num_cells_init = args_as_longs[3];
//     long num_gens = args_as_longs[4];
//     int num_env_factors = 1;
//     
//     
//     if (seed < 1)
//     	seed = time(0);
//     std::cerr << "Using seed of " << seed << '\n';
//     gingko::World world(seed);
// 
//     std::cerr << "(generating landscape)\n";
//     
// 	world.generate_landscape(size_x, size_y, num_env_factors);
// 	
//     std::cerr << "(setting carrying capacity)\n";
// 
// 	world.set_cell_carrying_capacity(cc);
// 
//     std::cerr << "(adding species)\n";
// 	
// 	gingko::Species& sp1 = world.new_species("gecko");
// 	
// 	std::vector<int> costs;
// 	costs.assign(size_x * size_y, 1);
// 	sp1.set_movement_costs(costs);
// 	sp1.set_movement_capacity(1);
// 	
// 	std::vector< gingko::FitnessFactorType > genotype;
// 	genotype.reserve(num_env_factors);
// 	for (int i = 0; i < num_env_factors; ++i) {
// 	    genotype.push_back(static_cast< gingko::FitnessFactorType >(world.rng().uniform_int(-10, 10)));
// 	}
//     
//     unsigned long max_index = (size_x * size_y)-1;
//     gingko::CellIndexType cell_index = 0;
//     for (std::set< gingko::CellIndexType > seeded; num_cells_init > 0; --num_cells_init) {
//         do {
//             cell_index = world.rng().uniform_int(0, max_index);
//         } while ((seeded.find(cell_index) != seeded.end()) and seeded.size() < max_index+1);
//         if (seeded.size() >= max_index+1) {
//             break;
//         }
//         std::cerr << "(seeding cell " << cell_index << " with " << cc << " individuals" << ")" << std::endl;
//         seeded.insert(cell_index);
//         world.seed_population(cell_index, sp1.get_index(), cc);
//     }
//     
//     std::cerr << "(running cycles)\n";
//     world.run(num_gens);
// 
// 
//     std::cerr << "\n#### FINAL STATUS ####\n"; 
//     world.landscape().dump(std::cerr);
//     
//     std::cerr << "\n\n#### TREE(S) ####\n";
//     gingko::Tree tree;
//     for (unsigned long i = (size_x * size_y); i != 0; --i) {
//         gingko::OrganismVector& ov = world.landscape()[i-1].organisms();
//         for (gingko::OrganismVector::iterator oiter = ov.begin();
//                 oiter != ov.end();
//                 ++oiter) {
//             tree.process_node(oiter->haploid_marker().node());
//         }                
//     }    
//     
// 
//     tree.dump(std::cerr);
//     std::cerr << "\n---\n\n";
//     tree.write_newick_tree(std::cout);
//     std::cerr << std::endl;
    
}
