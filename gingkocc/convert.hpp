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

#if !defined(GINGKO_CONVERT_H)
#define GINGKO_CONVERT_H

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace gingko { 
namespace convert {

/**
 * Exception thrown when conversion fails.
 */
class ValueError : public std::runtime_error {
    
    public:
        ValueError() : std::runtime_error("invalid literal for value type") {}
        ValueError(std::string value) : std::runtime_error("invalid literal for value type: \"" + value +"\"") {}
};

/**
 * Converts from one simple streamable type to another.
 * @param from  representation of value (e.g. "3.14")
 * @return      value represented in type T
 */
template <typename T, typename U>
T to_scalar(U from) {
    std::ostringstream o;
    T target;
    o << from;
    std::istringstream i(o.str());
    i >> target;
    if (i.fail() or not i.eof()) {
        throw ValueError(o.str());
    }
    return target;
}

/**
 * Converts from one simple streamable type to vector of streamable types.
 * @param from          representation of vector (e.g. "3 32 1 3 4", 
 *                      "1.2,1.3,3.1", "dda;adf;da" etc.)
 * @param separator     element delimiter (e.g, ",", " ", etc.)
 * @return              value represented in type T
 */
template <typename T, typename U>
T to_vector(U from, std::string separator) {
    std::ostringstream o;
    T target;
    o << from;
    std::istringstream i(o.str());
    i >> target;
    if (i.fail() or not i.eof()) {
        throw ValueError(o.str());
    }
    return target;
}

} // convert
} // gingko

#endif
