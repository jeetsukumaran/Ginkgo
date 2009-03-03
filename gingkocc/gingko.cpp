///////////////////////////////////////////////////////////////////////////////
//
// GINGKO Biogeographical and Evolution Simulator.
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

#include "random.h"
#include <iostream>

int main(int argc, char* argv[]) {
    long seed = 0;
	std::vector<int> args_as_ints;
	for (int i = 1; i < argc; ++i) {
		std::string a(argv[i]);
		if (a.length() > 2 && a[0] == '-') {
			if (a[1] == 's') {
				seed = atoi(a.c_str() + 2);
			}
			else {
				std::cerr << "Unknown flag: " << a << '\n';
				exit(2);
			}
		}
		else {
			args_as_ints.push_back(atoi(a.c_str()));
		}
	}
    if (args_as_ints.size() < 5) {
        std::cerr << "usage: " << argv[0] <<  " <DIM-X> <DIM-Y> <CELL-CARRYING-CAPACITY> <NUM-CELLS-TO-POPULATE> <NUM-GENS>\n";
        exit(1);
    }

    int size_x = args_as_ints[0];
    int size_y = args_as_ints[1];
    int cc = args_as_ints[2];
    int num_cells_init = args_as_ints[3];
    int num_gens = args_as_ints[4];
    int num_env_factors = 4;
    if (seed < 1)
    	seed = time(0);
    std::cerr << "Using seed of " << seed << '\n';
    World   world(seed);

//##DEBUG##
DEBUG_BLOCK( std::cerr << "(generating landscape)\n"; )
    
	world.generate_landscape(size_x, size_y, num_env_factors);
	
//##DEBUG##
DEBUG_BLOCK( std::cerr << "(setting carrying capacity)\n"; )

	world.set_cell_carrying_capacity(cc);

//##DEBUG##
DEBUG_BLOCK( std::cerr << "(adding species)\n"; )	
	
	Species& sp1 = world.new_species("gecko");
	
	std::vector<int> costs;
	costs.assign(size_x * size_y, 1);
	sp1.set_movement_costs(costs);
	sp1.set_movement_capacity(1);
	
	std::vector<FitnessFactorType> genotype;
	genotype.reserve(num_env_factors);
	for (int i = 0; i < num_env_factors; ++i) {
	    genotype.push_back(static_cast<FitnessFactorType>(world.rng().randint(-10, 10)));
	}
		
//##DEBUG##
DEBUG_BLOCK( std::cerr << "(seeding populations)\n"; )
    
    unsigned long max_index = (size_x * size_y)-1;
    CellIndexType cell_index = 0;
    for (std::set<CellIndexType> seeded; num_cells_init > 0; --num_cells_init) {
        do {
            cell_index = world.rng().randint(0, max_index);
        } while ((seeded.find(cell_index) != seeded.end()) and seeded.size() < max_index+1);
        if (seeded.size() >= max_index+1) {
            break;
        }
        seeded.insert(cell_index);
        world.seed_population(cell_index,
                              sp1.get_index(),
                              cc);
    }
    
//##DEBUG##
DEBUG_BLOCK( std::cerr << "(running cycles)\n"; )
    world.run(num_gens);


    std::cerr << "\n#### FINAL STATUS ####\n\n"; 
    world.landscape().dump(std::cerr);
}