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
#include <map>
#include <utility>
#include <stdexcept>
#include <vector>

#include "world.hpp"

#if !defined(GINGKO_CONFSYS_H)
#define GINGKO_CONFSYS_H

namespace gingko {

///////////////////////////////////////////////////////////////////////////////
// Client code should call one of the following to configure World objects.

/**
 * Build/populate a World object according to a given configuration source.
 * @param   world       World object to configure
 * @param   conf_src    input stream source of configuration file
 * @return              World object
 */
World& configure_world(World& world, std::istream& conf_src);

/**
 * Build/populate a World object according to a given configuration source file.
 * @param   world       World object to configure
 * @param   conf_fpath  path of configuration file
 * @return              World object
 */
World& configure_world(World& world, const char * conf_fpath);

/**
 * Build/populate a World object according to a given configuration source file.
 * @param   world       World object to configure 
 * @param   conf_fpath  path of configuration file
 * @return              World object
 */
World& configure_world(World& world, const std::string& conf_fpath);

///////////////////////////////////////////////////////////////////////////////
// Infrastructure supporting above functions.

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
 * Generic configuration block with a file.
 */
class ConfigurationBlock {

    public:
    
        /**
         * Default constructor.
         */
        ConfigurationBlock();
        
        /**
         * Constructor that populates a ConfigurationBlock object by calling
         * parse().
         *
         * @param in    input stream with data
         */
        ConfigurationBlock(std::istream& in);  
        
        /**
         * Default destructor.
         */
        ~ConfigurationBlock();
        
        /**
         * Clears existing data.
         */
        void clear();         
        
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
         * Returns <code>true</code> if the block was parsed and set.
         *
         * @return <code>true</code> if the block was parsed and set
         */
        bool is_block_set() const;        
        
        /**
         * Returns value for given key.
         *
         * @param       key key for entry
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
        
        /**
         * Allow for easy input semantics.
         * @param in        input stream         
         * @param cblock    block to populate
         * @return          input stream
         */
        friend std::istream& operator>> (std::istream& in, ConfigurationBlock& cblock);

    private:
        /** The type of block (e.g. "species", "world", "generation") */
        std::string                             type_;
        /** The name of the block (e.g. "Sp1"). */
        std::string                             name_;
        /** Key-value pairs making up the block body. */
        std::map< std::string, std::string >    entries_;
        /** Tracks whether or not the block was actually set. */
        bool                                    is_block_set_;
}; // ConfigurationBlock

/**
 * Encapsulates parsing of a configuration file, and populating of WorldConf,
 * SpeciesConf, GenerationConf, etc. objects.
 */
class ConfigurationFile {

    public:
        
        /**
         * Initializes metadata and binds to source stream.
         *
         * @param src   data source
         */
        ConfigurationFile(std::istream& src);
        
        /**
         * Initializes metadata and binds to source file.
         *
         * @param fpath filepath of data source
         */
        ConfigurationFile(const char * fpath);
        
        /**
         * Initializes metadata and binds to source file.
         *
         * @param fpath filepath of data source
         */
        ConfigurationFile(const std::string& fpath);        
        
        /**
         * Default no-op destructor.
         */
        ~ConfigurationFile();
        
        /**
         * Clears blocks.
         */
        void clear();        
        
        /**
         * Parses the configuration file, loading data into blocks.
         */
        void parse();
        
    private: 
    
        /** Input (file) stream. */
        std::ifstream       fsrc_;
        
        /** Input stream. */
        std::istream&       src_;
                
        /** Collection of World blocks. */
        ConfigurationBlock                world_;
        
        /** Collection of Species blocks. */
        std::vector<ConfigurationBlock>   species_;
        
        /** Collection of Generation blocks. */
        std::vector<ConfigurationBlock>   generations_;
        
};  

} // namespace gingko

#endif
