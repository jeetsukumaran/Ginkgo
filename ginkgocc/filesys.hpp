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

#include <string>

#if !defined(GINKGO_FILESYS_H)
#define GINKGO_FILESYS_H

namespace ginkgo {
namespace filesys {

///////////////////////////////////////////////////////////////////////////////
// PATH TEXT/STRING OPERATIONS

/**
 * Returns filename (and extension) from supplied path.
 *
 * Technically, returns the final element of a "/"-separated path.
 * @param  path path to file
 * @return      filename and extension
 */
std::string get_path_leaf(const std::string& path);

/**
 * Returns filename (and extension) from supplied path.
 *
 * Technically, returns the final element of a "/"-separated path.
 * @param  path path to file
 * @return      filename and extension
 */
std::string get_path_parent(const std::string& path);


/**
 * Returns filename (and extension) from supplied path.
 *
 * Technically, returns the final element of a "/"-separated path.
 * @param  path path to file
 * @return      filename and extension
 */
std::string get_path_leaf(const char * path);

/**
 * Returns filename (and extension) from supplied path.
 *
 * Technically, returns the final element of a "/"-separated path.
 * @param  path path to file
 * @return      filename and extension
 */
std::string get_path_parent(const char * path);

/**
 * Joins elements of a path.
 *
 * @param parent  parent path element
 * @param child   child path element
 * @return        fully composed path
 */
std::string compose_path(const std::string& parent, const std::string& child);

/**
 * Given a relative path assumed to be rooted in the current working directory,
 * return an absolute path.
 *
 * @param rel_path  path relative to current working directory
 * @return          absolute path
 */
std::string abs_path_from_cwd(const std::string& rel_path);

/**
 * Returns <code>true</code> if path is absolute.
 *
 * @param path      file path
 * @return          <code>true</code> if path is absolute
 */
bool is_abs_path(const std::string& path);

} // filesys

} // ginkgo

#endif
