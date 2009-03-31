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

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <cassert>

#include "../world.hpp"
#include "../confsys.hpp"
#include "../filesys.hpp"

using namespace gingko;

void test_parse_dummy1(const char * prog_path) {
    World   world;
    std::string srcf = filesys::compose_path(filesys::get_path_parent(prog_path), 
                                             "data/dummy1/dummy.conf");
    ConfigurationFile cf(srcf);
    cf.configure(world);
}

int main(int argc, char * argv[]) {
    assert(argc);
    test_parse_dummy1(argv[0]);
}


