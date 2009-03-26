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

#if !defined(GINGKO_CONF_H)
#define GINGKO_CONF_H

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <map>
#include <utility>
#include <stdexcept>

namespace gingko {

/**
 * General configuration format error.
 */
class ConfigurationParseError : public std::runtime_error {
    public:
        ConfigurationParseError(const char * msg) : std::runtime_error(msg) {}
        ConfigurationParseError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

/**
 * Generic configuration block with a file.
 */
class ConfigurationBlockParser {

    public:
    
        /**
         * Default constructor.
         */
        ConfigurationBlockParser();
        
        /**
         * Constructor that populates a ConfigurationBlock object by calling
         * parse().
         *
         * @param in    input stream with data
         */
        ConfigurationBlockParser(std::istream& in);  
        
        /**
         * Default destructor.
         */
        ~ConfigurationBlockParser();       
        
        /**
         * Composes exception message.
         *
         * @param pos   position in stream
         * @param desc  description of error 
         */
        std::string compose_error_message(unsigned long pos, const char * desc);
         
        /**
         * Composes exception message.
         *
         * @param pos   position in stream
         * @param desc  description of error 
         */
        std::string compose_error_message(unsigned long pos, const std::string& desc);                      

        /**
         * Populates this block object by reading from the given input stream
         * up to the next END_BLOCK_BODY character.
         *
         * @param in    input stream with data
         */
        void parse(std::istream& in);

    private:
        /** The type of block (e.g. "species", "world", "generation") */
        std::string                             type_;
        /** The name of the block (e.g. "Sp1"). */
        std::string                             name_;
        /** Key-value pairs making up the block body. */
        std::map< std::string, std::string >    entries_;
}; // ConfigurationBlockParser


} // namespace gingko

#endif
