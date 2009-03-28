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
#include <sstream>
#include <vector>
#include <algorithm>

#if !defined(GINGKO_TEXTUTILS_H)
#define GINGKO_TEXTUTILS_H

namespace gingko {
namespace textutils {

/**
 * Wraps a line of text to specified width.
 *
 * @param  source                   text to be wrapped
 * @param  line_width               width of wrapping
 * @param  first_line_indent        amount to indent first line
 * @param  subsequent_line_indent   amount to indent remaining lines
 * @return                          wrapped text (line breaks="\n")
 */
std::string textwrap(const std::string& source, 
        unsigned line_width=78,
        unsigned first_line_indent=0, 
        unsigned subsequent_line_indent=0);                         
     
/**
 * Returns filename (and extension) from supplied path.
 *
 * Technically, returns the final element of a "/"-separated path.
 * @param  source                   path to file
 * @return                          filename and extension
 */
std::string extract_filename_from_path(const char * path);

/**
 * Splits a const char source string into tokens as delimited by 
 * <code>sep</code>.
 *
 * @param src           source string to be split
 * @param sep           separator token (single or multiple characters)
 * @max_splits          maximum number of splits (0 = no limit)
 * @param include_empty_tokens treat consecutive delimiters as single?
 * @return      vector of std::string tokens
 */
std::vector<std::string> split(const char * ssrc, const char * sep = " ", unsigned max_splits=0, bool include_empty_tokens=true);

/**
 * Splits a std::string source string into tokens as delimited by 
 * <code>sep</code>.
 *
 * @param src   source string to be split
 * @param sep   separator token (single or multiple characters)
 * @max_splits          maximum number of splits (0 = no limit) 
 * @param include_empty_tokens treat consecutive delimiters as single?
 * @return      vector of std::string tokens
 */
std::vector<std::string> split(const std::string& src, const char * sep = " ", unsigned max_splits=0, bool include_empty_tokens=true);

/**
 * Splits a const char source string into tokens as delimited by any character
 * given in <code>sep</code>.
 *
 * @param src           source string to be split
 * @param sep           separator token (single or multiple characters)
 * @max_splits          maximum number of splits (0 = no limit)
 * @param include_empty_tokens treat consecutive delimiters as single?
 * @return      vector of std::string tokens
 */
std::vector<std::string> split_on_any(const char * ssrc, const char * sep = " ", unsigned max_splits=0, bool include_empty_tokens=true);

/**
 * Splits a const char source string into tokens as delimited by any character
 * given in <code>sep</code>.
 *
 * @param src   source string to be split
 * @param sep   separator token (single or multiple characters)
 * @max_splits          maximum number of splits (0 = no limit) 
 * @param include_empty_tokens treat consecutive delimiters as single?
 * @return      vector of std::string tokens
 */
std::vector<std::string> split_on_any(const std::string& src, const char * sep = " ", unsigned max_splits=0, bool include_empty_tokens=true);


/**
 * Returns copy of given string with specified characters stripped from
 * beginning and end.
 *
 * @param s         source string
 * @param to_strip  characters to remove
 * @return          copy of source with specified characters removed from
 *                  beginning and end
 */
std::string strip(const std::string& s, const char * to_strip=" \t\n\r");

/** 
 * Converts a string to lower case.
 *
 * @param s         source string
 * @return          copy of string in all lower case characters.
 */
std::string lower(const std::string& s);

/** 
 * Converts a string to upper case.
 *
 * @param s         source string
 * @return          copy of string in all upper case characters.
 */
std::string upper(const std::string& s);

} // namespace textutils
} // namespace gingko

#endif
