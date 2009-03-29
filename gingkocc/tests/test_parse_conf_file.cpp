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

#include "../worldconf.h"
#include "../filesys.h"

using namespace gingko;

void test_parse_dummy1(const char * prog_path) {
    std::string srcf = filesys::compose_path(prog_path, "data/dummy1/dummy1.conf");
    ConfigurationFile cf(srcf);
}

int main(int argc, char * argv[]) {
    assert(argc);
    test_parse_dummy1(argv[0]);
}


