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

#include "gingko_defs.hpp"
#include "confsys.hpp"
#include "biosys.hpp"
#include "world.hpp"
#include "tree.hpp"
#include "cmdopt.hpp"
#include "filesys.hpp"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <set>
#include <ctime>
#include <string>
#include <sstream>


int main(int argc, char* argv[]) {

    std::string config_filepath;
    std::string output_dir = ".";
    std::string replicate_id;
    bool validate_config_only = false;
    unsigned long rand_seed = 0;

    gingko::OptionParser parser = gingko::OptionParser("Gingko 0.1.3",
            "Gingko Biogeographical Evolution Simulator",
            "%prog [options] <CONFIGURATION-FILEPATH>");
    
    parser.add_option<std::string>(&replicate_id, "-i", "--replicate-id", 
                               "optional identifier to add to names of output files produced during this run");
    parser.add_option<std::string>(&output_dir, "-o", "--output-dir", 
                                   "directory to which to save output files (default = current)");
    parser.add_switch(&validate_config_only, "-v", "--validate",
                      "load and process configuration file to check for errors, but do not actually execute run");
    parser.add_option<unsigned long>(&rand_seed, "-z", "--random-seed", 
                               "random seed for this run (overrides seed set in configuration file)");

    parser.parse(argc, argv);       
    
    std::vector<std::string> args = parser.get_args();
    
    if (args.size() == 0) {
        parser.write_usage(std::cout);
        exit(1);
    }
    
    gingko::World world;
    world.set_replicate_id(replicate_id);
    world.set_output_dir(output_dir);
       
    gingko::confsys::configure_world(world, args[0]);
    if (validate_config_only) {
        std::cout << "World configured using: \"" + args[0] + "\"" << std::endl; 
        std::cout << "Configuration file validates." << std:: endl;
    } else {
        world.open_logs();
        world.log_info("World configured using: \"" + args[0] + "\"");
        
        std::ostringstream seed_info;
        if (parser.is_set("--random-seed")) {
            seed_info << "Setting random seed from command line: " << rand_seed << ".";
            world.set_random_seed(rand_seed);
        } else {
            seed_info << "Using random seed: " << world.get_random_seed() << "."; 
        }
        world.log_info(seed_info.str());
        world.run();
        world.close_logs();
    }
    
}
