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

#define assertm(cond, out, msg) if (!(cond)) { out << msg; abort(); }

#include "../tree.h"
#include "../cmdopts.h"
#include "../textutils.h"
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {

    gingko::OptionParser parser = gingko::OptionParser("Tree Testing",
            NULL, 
            "%prog [options] <NODE>[:<LABEL>] [<NODE>[:<LABEL>] [<NODE>[:<LABEL>] [<NODE>[:<LABEL>] ... ]]]");

    parser.parse(argc, argv);
    std::vector< std::string > args = parser.get_args();
    if (args.size() == 0) {
    
    }

}
