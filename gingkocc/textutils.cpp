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
#include <cassert>
#include <string>
#include <sstream>
#include <vector>
#include <cctype>
#include "textutils.h"

namespace gingko {
namespace textutils {


// Wraps text (preferably at word boundaries).
std::string textwrap(const std::string& source, 
        unsigned line_width,
        unsigned first_line_indent, 
        unsigned subsequent_line_indent)  {
    std::string wrapped;
    unsigned col_count = 1;
    unsigned line_count = 1;
    std::string subsequent_line_indent_str(subsequent_line_indent, ' ');
    for (std::string::const_iterator s = source.begin();
            s != source.end();
            ++s, ++col_count) {

        if (*s == '\n') {
            wrapped += "\n";
            col_count = 0;
            line_count += 1;
            continue;
        }

        if (col_count > line_width) {
            std::string::size_type last_break = wrapped.find_last_of("\n");
            std::string::size_type wrap_pos = wrapped.find_last_of(" ");
            if (wrap_pos == std::string::npos or ((last_break != std::string::npos) and (last_break > wrap_pos))) {
                wrapped += "\n";
                col_count = 1;
            } else {                
                wrapped.replace(wrap_pos, 1, "\n" + subsequent_line_indent_str);             
                col_count = wrapped.size() - wrap_pos;                    
            }
        }
            
        if (col_count == 1 and line_count == 1 and first_line_indent > 0) {
            for (unsigned i = 0; i < first_line_indent; ++i) {
                wrapped += ' ';
            }
            col_count += first_line_indent;
        } else if (col_count == 1 and line_count > 1) {   
            wrapped += subsequent_line_indent_str;
            col_count += subsequent_line_indent;                    
        }
        wrapped += *s;

    }

    return wrapped;
} 


// Extracts filenames from path
std::string extract_filename_from_path(const char * path) {
    
    // copy of string
    std::string full_path = path;
    
    // replace dos/windows slashes with nix ones
    // if there are insane filepaths being passed here (specifically, a path
    // with backslash characters), this will give wrong results
    std::string::size_type p = full_path.find('\\');
    while (p != std::string::npos) {
        full_path.replace(p, 1, "/");
    }
    std::string::size_type last_path_char = full_path.find_last_of('/');
    if (last_path_char == std::string::npos) {
        return full_path;
    } else {
        if (last_path_char >= full_path.size()) {
            return std::string("");
        } else {
            return full_path.substr(last_path_char+1);
        }
    }    
} 


// Split a string by given separator delimiter
std::vector<std::string> split(const char * ssrc, 
                               const char * sep, 
                               unsigned max_splits, 
                               bool include_empty_tokens) {
    return split(std::string(ssrc), sep, max_splits, include_empty_tokens);
}


// Split a string by given separator delimiter
std::vector<std::string> split(const std::string& src, 
                               const char * sep, 
                               unsigned max_splits, 
                               bool include_empty_tokens) {
    std::vector< std::string > v;
    std::string::size_type start_pos = 0;
    std::string::size_type end_pos = src.find(sep, start_pos);
    unsigned num_splits = 0;
    while (end_pos != std::string::npos and (max_splits == 0 or num_splits < max_splits)) {        
        std::string result = src.substr(start_pos, end_pos-start_pos);
        if (result.size() != 0 or include_empty_tokens) {
            num_splits += 1;
            v.push_back(result);
        }            
        start_pos = end_pos+1;
        end_pos = src.find(sep, start_pos);
    }
    std::string result = src.substr(start_pos, std::string::npos);
    if (result.size() != 0 or include_empty_tokens) {
        v.push_back(result);
    } 
    return v;
}


// Split a string by any character in given list of separators/delimiter
std::vector<std::string> split_on_any(const char * ssrc, 
                                      const char * sep, 
                                      unsigned max_splits, 
                                      bool include_empty_tokens) {
    return split_on_any(std::string(ssrc), sep, max_splits, include_empty_tokens);
}


// Split a string by any character in given list of separators/delimiter
std::vector<std::string> split_on_any(const std::string& src, 
                                      const char * sep, 
                                      unsigned max_splits, 
                                      bool include_empty_tokens) {
    std::vector< std::string > v;
    std::string::size_type start_pos = 0;
    std::string::size_type end_pos = src.find_first_of(sep, start_pos);
    unsigned num_splits = 0;
    while (end_pos != std::string::npos and (max_splits == 0 or num_splits < max_splits)) {        
        std::string result = src.substr(start_pos, end_pos-start_pos);
        if (result.size() != 0 or include_empty_tokens) {
            num_splits += 1;
            v.push_back(result);
        }            
        start_pos = end_pos+1;
        end_pos = src.find_first_of(sep, start_pos);
    }
    std::string result = src.substr(start_pos, std::string::npos);
    if (result.size() != 0 or include_empty_tokens) {
        v.push_back(result);
    } 
    return v;
}


// Strip characters from a string
std::string strip(const std::string& s, const char * to_strip) {
    if (!s.empty()) {
        std::size_t start = s.find_first_not_of(to_strip);
        if (start == std::string::npos) {
            return "";
        }
        std::size_t end = s.find_last_not_of(to_strip);
        return s.substr(start, end-start+1);
    } else {
        return s;
    }
}

std::string lower(const std::string& s) {
    std::string result = s;
    std::transform(s.begin(), s.end(), result.begin(), (int(*)(int))std::tolower);
    return result;
}

std::string upper(const std::string& s) {
    std::string result = s;
    std::transform(s.begin(), s.end(), result.begin(), (int(*)(int))std::toupper);
    return result;
}

} // textutils namespace
} // gingko namespace
