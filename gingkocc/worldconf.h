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

#if !defined(GINGKO_WORLDCONF_H)
#define GINGKO_WORLDCONF_H

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <map>
#include <utility>
#include <stdexcept>
#include <vector>

namespace gingko {

/**
 * General i/o error.
 */
class ConfigurationIOError : public std::runtime_error {
    public:
        ConfigurationIOError(const char * msg) : std::runtime_error(msg) {}
        ConfigurationIOError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

/**
 * General configuration format error.
 */
class ConfigurationSyntaxError : public std::runtime_error {
    public:
        ConfigurationSyntaxError(const char * msg) : std::runtime_error(msg) {}
        ConfigurationSyntaxError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

/**
 * Encapsulates parsing of a configuration file, and instantiation of 
 * corresponding World object.
 */
class ConfigurationFileParser {

    public:
        
        /**
         * Initializes metadata and binds to source stream.
         *
         * @param src   data source
         */
        ConfigurationFileParser(std::istream& src);
        
        /**
         * Initializes metadata and binds to source file.
         *
         * @param fpath filepath of data source
         */
        ConfigurationFileParser(const char * fpath);
        
        /**
         * Initializes metadata and binds to source file.
         *
         * @param fpath filepath of data source
         */
        ConfigurationFileParser(const std::string& fpath);        
        
        /**
         * Default no-op destructor.
         */
        ~ConfigurationFileParser();
        
    private: 
    
        /** Input (file) stream. */
        std::ifstream       fsrc_;
        
        /** Input stream. */
        std::istream&       src_;     
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
         * Returns type label of the block.
         *
         * @return type label of block
         */
        std::string get_type() const;

        /**
         * Returns type label of the block.
         *
         * @return type label of block
         */
        std::string get_name() const;
        
        /**
         * Returns value for given key.
         *
         * @param   key key for entry
         * @return      value for entry
         */
        std::string get_entry(const std::string& key) const;   
        
        /**
         * Returns vector of keys in entries.
         * @return      vector of keys in entries
         */
        std::vector<std::string> get_keys() const;
        
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
