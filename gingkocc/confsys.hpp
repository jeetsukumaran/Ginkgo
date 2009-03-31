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
#include "convert.hpp"

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
 * General configuration error.
 */
class ConfigurationError : public std::runtime_error {
    public:
        ConfigurationError(const char * msg) : std::runtime_error(msg) {}
        ConfigurationError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

/**
 * General i/o error.
 */
class ConfigurationIOError : public ConfigurationError {
    public:
        ConfigurationIOError(const char * msg) : ConfigurationError(msg) {}
        ConfigurationIOError(const std::string& msg) : ConfigurationError(msg) {}
};

/**
 * General configuration format error.
 */
class ConfigurationSyntaxError : public ConfigurationError {
    public:
        ConfigurationSyntaxError(const char * msg) : ConfigurationError(msg) {}
        ConfigurationSyntaxError(const std::string& msg) : ConfigurationError(msg) {}
};

/**
 * General configuration content error.
 */
class ConfigurationIncompleteError : public ConfigurationError {
    public:
        ConfigurationIncompleteError(const char * msg) : ConfigurationError(msg) {}
        ConfigurationIncompleteError(const std::string& msg) : ConfigurationError(msg) {}
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
         * @return type label of block
         */
        std::string get_type() const;

        /**
         * Returns type label of the block.
         * @return type label of block
         */
        std::string get_name() const;
        
        /**
         * Returns <code>true</code> if the block was parsed and set.
         * @return <code>true</code> if the block was parsed and set
         */
        bool is_block_set() const;        
        
        /**
         * Returns value for given key.
         * @param   key key for entry
         * @return      value for entry
         */
        template <typename T>
        T get_entry(const std::string& key) const {
            std::map< std::string, std::string >::const_iterator val = this->entries_.find(key);
            if (val == this->entries_.end()) {
                throw ConfigurationIncompleteError("entry \"" + key + "\" not found");
            }
            return convert::to_type<T>(val->second);
        }
        
        /**
         * Returns value for given key (with specification of default value if 
         * not found.
         * @param   key             key for entry
         * @param   default_value   value to return if key not found
         * @return                  value for entry or default
         */
        template <typename T>
        T get_entry(const std::string& key, T default_value) const {
            std::map< std::string, std::string >::const_iterator val = this->entries_.find(key);
            if (val == this->entries_.end()) {
                return default_value;
            }
            return convert::to_type<T>(val->second);
        }
        
        /**
         * Returns <code>true</code> if has entry with specified key.
         * @param   key     key for entry
         * @return          <code>true</code> if has entry with specified key
         */
        bool has_key(const std::string key) const {
            std::map< std::string, std::string >::const_iterator val = this->entries_.find(key);
            return (val != this->entries_.end());
        }
        
        /**
         * Returns vector of keys in entries.
         * @return      vector of keys in entries
         */
        std::vector<std::string> get_keys() const;
        
        /**
         * Composes exception message.
         * @param pos   position in stream
         * @param desc  description of error 
         */
        std::string compose_error_message(unsigned long pos, const char * desc);
         
        /**
         * Composes exception message.
         * @param pos   position in stream
         * @param desc  description of error 
         */
        std::string compose_error_message(unsigned long pos, const std::string& desc);                      

        /**
         * Populates this block object by reading from the given input stream
         * up to the next END_BLOCK_BODY character.
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
 * Base class for configurators.
 */
class Configurator {
    
    public:
    
        /** 
         * Stores variables for error reporting. 
         * @param cb                configuration data parsed into 
         *                          ConfigurationBlock structure
         * @param block_start_pos   start position of this block in the stream 
         *                          that is the source of the configuration 
         *                          data (for error reporting)
         * @param block_end_pos     start position of this block in the stream 
         *                          that is the source of the configuration 
         *                          data (for error reporting)         
         */
        Configurator(const ConfigurationBlock& cb,
                     unsigned long block_start_pos, 
                     unsigned long block_end_pos);
                     
        virtual ~Configurator();                     
                     
        /** 
         * Takes the string fields of ConfigurationBlock and interprets values
         * as needed for a World object.   
         */        
        virtual void parse() = 0;
        
        /**
         * Retrieves value for specified key.
         * @param key   entry in the ConfigurationBlock dictionary
         * @return      value for key
         */
        template <typename T>
        T get_configuration_value(const std::string& key) {
            try {
                return this->configuration_block_.get_entry<T>(key);
            } catch (ConfigurationIncompleteError e) {
                throw this->build_exception(e.what());
            } catch (convert::ValueError e) {
                throw this->build_exception(std::string(e.what()) + " (specified for \"" + key + "\")");
            }            
        }
        
        /**
         * Retrieves value for specified key, with default value returned if key
         * is not found.
         * @param key           key for entry in the ConfigurationBlock dictionary
         * @param default_value value if key not found
         * @return              value for key in entries or default value
         */
        template <typename T>
        T get_configuration_value(const std::string& key, T default_value) {
            try {
                return this->configuration_block_.get_entry<T>(key, default_value);
            } catch (convert::ValueError e) {
                throw this->build_exception(std::string(e.what()) + " (specified \"" + key + "\")");
            }            
        }        
        
        /**
         * Builds and returns an exception, including file position in 
         * exception message.
         *
         * @param message   information regarding error
         * @return          exception object
         */
        ConfigurationError build_exception(const std::string& message) const;       

    private:
        /** Underlying configuration block. */        
        ConfigurationBlock            configuration_block_;
        /** For error messages. */
        unsigned long                 block_start_pos_;
        /** For error messages. */
        unsigned long                 block_end_pos_;        
        
}; // Configurator

/**
 * Takes a ConfigurationBlock assumed to be wrapped around World information, 
 * and parses/translates values appropriately.
 */
class WorldConfigurator : public Configurator {

    public:
    
        /** 
         * Constructs objects, and then passes ConfigurationBlock onto parse()
         * for processing. 
         * @param cb                a populated ConfigurationBlock object
         * @param block_start_pos   start position of this block in the stream 
         *                          that is the source of the configuration 
         *                          data (for error reporting)
         * @param block_end_pos     start position of this block in the stream 
         *                          that is the source of the configuration 
         *                          data (for error reporting)         
         */
        WorldConfigurator(const ConfigurationBlock& cb, 
                  unsigned long block_start_pos, 
                  unsigned long block_end_pos);

        /** 
         * Takes the string fields of ConfigurationBlock and interprets values
         * as needed for a World object.
         */        
        void parse();
        
        /**
         * Configures a World object according to settings.
         */
        void configure(World& world);         
           
    private:
        /** Size in x-dimension (number of columns). */
        CellIndexType   size_x_;
        /** Size in y-dimension (number of rows). */
        CellIndexType   size_y_;
        /** Number of fitness factors. */
        unsigned int    num_fitness_factors_;
        /** Random number seed. */        
        unsigned long   rand_seed_;

}; // WorldConfigurator


/**
 * Encapsulates parsing of a configuration file, and populating of WorldConfigurator,
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
                
        /** Collection of WorldConfigurator blocks. */
        std::vector<WorldConfigurator>      worlds_;
        
        /** Collection of Species blocks. */
        std::vector<ConfigurationBlock>   species_;
        
        /** Collection of Generation blocks. */
        std::vector<ConfigurationBlock>   generations_;
        
};  

} // namespace gingko

#endif
