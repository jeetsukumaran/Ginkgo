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
#include <ctime>

int main(int argc, char* argv[]) {
    unsigned long size_x = 1000;
    unsigned long size_y = 1000;
    unsigned long cc = 100;
    unsigned long num_gens = 100000;
    unsigned int num_fitness = 10;
    unsigned long rand_seed = time(0);
    unsigned int num_cells_init = 1;

    gingko::OptionParser parser = gingko::OptionParser("Gingko Biogeographical Evolution Simulator", 
            "Version 0.01", 
            "%prog [options]");
    
    parser.add_option<unsigned long>(&size_x, "-x", "--dim-x", 
                                     "size of landscape in the x-dimension (default = %default)", "DIM-X");
    parser.add_option<unsigned long>(&size_y, "-y", "--dim-y", 
                                     "size of landscape in the y-dimension (default = %default)", "DIM-Y");
    parser.add_option<unsigned long>(&cc, "-c", "--carrying-capacity", 
                                     "maximum carrying capacity of each cell (default = %default)", "K");
    parser.add_option<unsigned long>(&num_cells_init, "-i", "--init-cells", 
                                     "number of cells to seed with initial population (default = %default)", "#");
    parser.add_option<unsigned long>(&num_gens, "-g", "--num-gens",                                     
                                     "number of generations to run (default = %default)", "#GENERATIONS");  
    parser.add_option<unsigned int>(&num_fitness, "-f", "--num-fitness", 
                                    "number of fitness factors (default = %default)", "#FACTORS");
    parser.add_option<unsigned long>(&rand_seed, "-z", "--random-seed", 
                                     "random number seed (default = %default)", "SEED");                                
                                     
    parser.parse(argc, argv);       
    
    std::cerr << "           Landscape size: (" << size_x << ", " << size_y << ")\n";
    std::cerr << "   Cell carrying capacity: " << cc << '\n';
    std::cerr << "    Number of generations: " << num_gens << '\n';
    std::cerr << "Number of fitness factors: " << num_fitness << '\n';
    std::cerr << "       Random number seed: " << rand_seed << std::endl;
    
    gingko::World world(rand_seed);

    std::cerr << "(generating landscape)\n";    
	world.generate_landscape(size_x, size_y, num_fitness);	
    std::cerr << "(setting carrying capacity)\n";
	world.set_cell_carrying_capacity(cc);
    std::cerr << "(adding species)\n";	
	gingko::Species& sp1 = world.new_species("Sp1");
	
	std::vector<int> costs;
	costs.assign(size_x * size_y, 1);
	sp1.set_movement_costs(costs);
	sp1.set_movement_capacity(1);
	
	std::vector< gingko::FitnessFactorType > genotype;
	genotype.reserve(num_fitness);
	for (unsigned i = 0; i < num_fitness; ++i) {
	    genotype.push_back(static_cast< gingko::FitnessFactorType >(world.rng().uniform_int(-10, 10)));
	}
    
    unsigned long max_index = (size_x * size_y)-1;
    gingko::CellIndexType cell_index = 0;
    for (std::set< gingko::CellIndexType > seeded; num_cells_init > 0; --num_cells_init) {
        do {
            cell_index = world.rng().uniform_int(0, max_index);
        } while ((seeded.find(cell_index) != seeded.end()) and seeded.size() < max_index+1);
        if (seeded.size() >= max_index+1) {
            break;
        }
        std::cerr << "(seeding cell " << cell_index << " with " << cc << " individuals" << ")" << std::endl;
        seeded.insert(cell_index);
        world.seed_population(cell_index, sp1.get_index(), cc);
    }
    
    std::cerr << "(running cycles)\n";
    world.run(num_gens);


    std::cerr << "\n#### FINAL STATUS ####\n"; 
    world.landscape().dump(std::cerr);
    
    std::cerr << "\n\n#### TREE(S) ####\n";
    gingko::Tree tree;
    for (unsigned long i = (size_x * size_y); i != 0; --i) {
        gingko::OrganismVector& ov = world.landscape()[i-1].organisms();
        for (gingko::OrganismVector::iterator oiter = ov.begin();
                oiter != ov.end();
                ++oiter) {
            std::string label = sp1.get_organism_label(*oiter);
            tree.process_node(oiter->haploid_marker().node(), &label);
        }                
    }    
    

//     tree.dump(std::cerr);
    std::cerr << "\n---\n\n";
    tree.write_newick_tree(std::cout);
    std::cout << std::endl;
    
}
