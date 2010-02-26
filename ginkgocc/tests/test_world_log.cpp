///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Biogeographical Evolution Simulator.
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

#include "../randgen.hpp"
#include "../world.hpp"

using namespace ginkgo;

int main(int, char * []) {
    World& world = World::get_instance();
    world.set_output_dir("/tmp");
    world.init_logger();
    world.logger().debug("This is a DEBUG level message");
    world.logger().info("This is a INFO level message");
    world.logger().warning("This is a WARNING level message");
    world.logger().error("This is a ERROR level message");
    world.logger().critical("This is a CRITICAL level message");
}
