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

#include <string>

#if !defined(GINGKO_FILESYS_H)
#define GINGKO_FILESYS_H

namespace gingko {
namespace filesys {
     
/**
 * Returns filename (and extension) from supplied path.
 *
 * Technically, returns the final element of a "/"-separated path.
 * @param  path path to file
 * @return      filename and extension
 */
std::string extract_filename_from_path(const char * path);

/**
 * Returns the current working directory.
 * @return      current working directory
 */
std::string current_working_dir();

} // filesys
} // gingko

#endif
