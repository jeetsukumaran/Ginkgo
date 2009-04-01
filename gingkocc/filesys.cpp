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

#include <cstdlib>
#include <string>
#include "filesys.hpp"

namespace gingko {
namespace filesys {

const char * PATH_SEPARATOR = "/";

///////////////////////////////////////////////////////////////////////////////
// PATH TEXT/STRING OPERATIONS

// extracts filenames from path
std::string get_path_leaf(const std::string& path) {
    std::string::size_type last_path_char = path.find_last_of(PATH_SEPARATOR);
    if (last_path_char == std::string::npos) {
        return path;
    } else if (last_path_char == path.size()-1) {
        return get_path_leaf(path.substr(0, path.size()-1).c_str());
    } else {
        if (last_path_char >= path.size()) {
            return std::string("");
        } else {
            return path.substr(last_path_char+1);
        }
    }    
} 

// extracts directory from path
std::string get_path_parent(const std::string& path) {
    std::string::size_type last_path_char = path.find_last_of(PATH_SEPARATOR);
    if (last_path_char == std::string::npos) {
        return std::string("");
    } else if (last_path_char == path.size()-1) {
        return get_path_parent(path.substr(0, path.size()-1).c_str());
    } else {
        if (last_path_char >= path.size()) {
            return std::string("");
        } else {
            return path.substr(0, last_path_char);
        }
    }     
}

std::string get_path_leaf(const char * path) {
    return get_path_leaf(std::string(path));
}

std::string get_path_parent(const char * path) {
    return get_path_parent(std::string(path));
}


// put together file path
std::string compose_path(const std::string& parent, const std::string& child) {
    return parent + PATH_SEPARATOR + child;
}

// join rel_path to current working directory
std::string abs_path_from_cwd(const std::string& rel_path) {
    return current_path() + PATH_SEPARATOR + rel_path;
}

///////////////////////////////////////////////////////////////////////////////
// OPERATING/FILE SYSTEM INTERACTIONS

// returns current working directory
std::string current_path() {
    char cwd[1024];
    ::getcwd(cwd, 1024);
    return std::string(cwd);
}

} // filesys
} // gingko
