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


int main(int argc, char* argv[]) {

    std::string config_filepath;
    std::string output_dir = ".";
    bool dry_run;

    gingko::OptionParser parser = gingko::OptionParser("Gingko 0.01",
            "Gingko Biogeographical Evolution Simulator",
            "%prog [options] <CONFIGURATION-FILEPATH>");
    
    parser.add_option<std::string>(&output_dir, "-o", "--output-dir", 
                                   "directory to which to save output files (default = current)");
    parser.add_switch(&dry_run, "-n", "--dry-run",
                      "load configuration file, but do not actually produce output (can be used to validate configuration file)");                               
                                     
    parser.parse(argc, argv);       
    
    std::vector<std::string> args = parser.get_args();
    
    if (args.size() == 0) {
        parser.write_usage(std::cout);
        exit(1);
    }
    
    gingko::World world;
    world.set_output_dir(output_dir);
    if (not dry_run) {
        world.open_logs();
        world.log_info("Processing configuration file: \"" + args[0] + "\"."); 
    }        
    gingko::confsys::configure_world(world, args[0]);
    if (not dry_run) {
        world.log_info("World configured.");
        world.log_info("Starting simulation.");
        
        world.log_info("Ending simulation.");
    }
    world.log_info("Exiting simulation.");
}
