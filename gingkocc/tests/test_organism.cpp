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
#include "../random.h"

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
    
    std::cout << "\nGenerating gen1 ...\n";
    std::vector<Organism> gen1;
    for (int i = 0; i < 4; ++i) {
        gen1.push_back(sp.new_organism());
    }
    for (int i = 0; i < 4; ++i) {
        gen1[i].dump(std::cout);
        std::cout << std::endl;
    }   
    
    std::cout << "\nGenerating gen2 ...\n";
    std::vector<Organism> gen2;
    gen2.push_back(sp.new_organism(gen1[0], gen1[1]));
    gen2.push_back(sp.new_organism(gen1[1], gen1[2]));
    gen2.push_back(sp.new_organism(gen1[2], gen1[3]));
    
    std::cout << "\nDestroying gen1 ...\n";    
    gen1.erase(gen1.begin(), gen1.end());

    for (int i = 0; i < 3; ++i) {
        gen2[i].dump(std::cout);
        std::cout << std::endl;
    }    
    
    std::cout << "\nGenerating gen3 ...\n";
    std::vector<Organism> gen3;
    gen3.push_back(sp.new_organism(gen2[1], gen2[2]));
    
    std::cout << "\nDestroying gen2 ...\n";       
    gen2.erase(gen2.begin(), gen2.end());

    for (int i = 0; i < 1; ++i) {
        gen3[i].dump(std::cout);
        std::cout << std::endl;
    }    
    
    
}