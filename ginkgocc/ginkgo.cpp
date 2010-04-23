///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Phylogeographical Evolution Simulator.
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
#include "ginkgo_defs.hpp"
#include "confsys.hpp"
#include "organism.hpp"
#include "population.hpp"
#include "species.hpp"
#include "world.hpp"
#include "tree.hpp"
#include "cmdopt.hpp"
#include "filesys.hpp"
#include "randgen.hpp"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <set>
#include <ctime>
#include <string>
#include <sstream>
#include <cstring>
#ifdef HAVE_CONFIG_H
#	include <config.h>
#   include "ginkgo_build_id.h"
#else
#	include "win_config.h"
#endif

std::string get_program_identification() {
    std::ostringstream s;
    s << PACKAGE_NAME << " v" << PACKAGE_VERSION;
#if defined(BUILDDESC) || defined(BUILDTIMESTAMP)
    s << " (";
    #ifdef BUILDDESC
        s << BUILDDESC;
        #ifdef BUILDTIMESTAMP
            s << ", ";
        #endif
    #endif
    #ifdef BUILDTIMESTAMP
        s << BUILDTIMESTAMP;
    #endif
    s << ")";
#endif
    return  s.str();
}

int main(int argc, char* argv[]) {

    std::string config_filepath;
    std::string output_dir = ".";
    std::string replicate_id;
    bool validate_config_only = false;
    unsigned long rand_seed = 0;
    unsigned long log_freq = 10;
    std::string program_ident = get_program_identification();

    ginkgo::OptionParser parser = ginkgo::OptionParser(
             program_ident.c_str(),
            "GINKGO Phylogeographical Evolution Simulator",
            "%prog [options] <CONFIGURATION-FILEPATH>");

    parser.add_option<std::string>(&replicate_id, "-i", "--replicate-id",
                               "optional identifier to add to names of output files produced during this run");
    parser.add_option<std::string>(&output_dir, "-o", "--output-dir",
                                   "directory to which to save output files (default = current)");
    parser.add_switch(&validate_config_only, "-v", "--validate",
                      "load and process configuration file to check for errors, but do not actually execute run");
    parser.add_option<unsigned long>(&rand_seed, "-z", "--random-seed",
                               "random seed for this run (overrides seed set in configuration file)");
    parser.add_option<unsigned long>(&log_freq, "-l", "--log-frequency",
                               "number of generations between each routine log message (default = 10)");

    parser.parse(argc, argv);

    std::vector<std::string> args = parser.get_args();

    if (args.size() == 0) {
        std::cerr << program_ident << std::endl << std::endl;
        parser.write_usage(std::cerr);
        exit(1);
    }

    ginkgo::World& world = ginkgo::World::get_instance();
    world.set_replicate_id(replicate_id);
    world.set_output_dir(output_dir);

    ginkgo::confsys::configure_world(world, args[0]);
    if (validate_config_only) {
        std::cout << program_ident << std::endl;
        std::cout << "World configured using: \"" + args[0] + "\"" << "." << std::endl;
        std::cout << "Configuration file validates." << std:: endl;
    } else {
        world.init_logger();
        world.logger().hide_elapsed_time();
        world.logger().info("Starting " + program_ident);
        world.logger().info("World configured using: \"" + args[0] + "\"");

        std::ostringstream seed_info;
        if (parser.is_set("--random-seed")) {
            seed_info << "Setting random seed from command line: " << rand_seed << ".";
            world.set_random_seed(rand_seed);
        } else {
            seed_info << "Using random seed: " << world.get_random_seed() << ".";
        }
        world.logger().info(seed_info.str());

        if (parser.is_set("--log-frequency")) {
            world.set_log_frequency(log_freq);
        }

        world.run();

        world.logger().hide_elapsed_time();
        std::ostringstream bbye;
        bbye << "Ending " << program_ident << ", run on '" << args[0] << "'";
        bbye << " with total run time of " << world.logger().get_hours_since_launched() << " hours.";
        world.logger().info(bbye.str());
    }
}
