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
            return convert::to_scalar<T>(val->second);
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
            return convert::to_scalar<T>(val->second);
        }
        
        /**
         * Returns <code>true</code> if has entry with specified key.
         * @param   key     key for entry
         * @return          <code>true</code> if has entry with specified key
         */
        bool has_key(const std::string key) const {
            return (this->entries_.find(key) != this->entries_.end());
        }
        
        /**
         * Returns vector of keys in entries.
         * @return      vector of keys in entries
         */
        std::vector<std::string> get_keys() const;
        
        /**
         * Returns vector of keys in entries that start with a particular 
         * string.
         * @param   key_start   string matching beginning of key
         * @return              vector of keys in entries
         */
        std::vector<std::string> get_keys(const std::string& key_start) const;
        
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
        
        /**
         * Returns the file character position of the beginning of this block 
         * in the source stream.
         * @return      file position of the start of the block
         */
        unsigned long get_block_start_pos() const;
        
        /**
         * Returns the file character position of the end of this block 
         * in the source stream.
         * @return      file position of the end of the block
         */
        unsigned long get_block_end_pos() const;        

    private:
        /** The type of block (e.g. "species", "world", "generation") */
        std::string                             type_;
        /** The name of the block (e.g. "Sp1"). */
        std::string                             name_;
        /** Key-value pairs making up the block body. */
        std::map< std::string, std::string >    entries_;
        /** Tracks whether or not the block was actually set. */
        bool                                    is_block_set_;
        /** Start of block in the source stream. */
        unsigned long                           block_start_pos_;
        /** End of block in the source stream. */
        unsigned long                           block_end_pos_;
        
}; // ConfigurationBlock

/**
 * Track distribution of organisms, either for seed populations or sampling 
 * regimes.
 */
struct OrganismDistribution {
    
    public:
        OrganismDistribution()
            : num_organisms(0) { }

    public:
        std::string     species_label;
        unsigned long   num_organisms;
        std::vector<CellIndexType>  x;
        std::vector<CellIndexType>  y;
        

}; // OrganismDistribution

/**
 * Base class for configurators.
 */
class Configurator {
    
    public:
    
        /** 
         * Stores variables for error reporting. 
         * @param cb                configuration data parsed into 
         *                          ConfigurationBlock structure      
         */
        Configurator(const ConfigurationBlock& cb);
                     
        virtual ~Configurator();                     
                     
        /** 
         * Takes the string fields of ConfigurationBlock and interprets values
         * as needed for a World object.   
         */        
        virtual void parse() = 0;
        
        /**
         * Returns name of underlying ConfigurationBlock.
         * @return  name of block
         */
        std::string get_name() const {
            return this->configuration_block_.get_name();
        }
        
        /**
         * Retrieves value for specified key.
         * @param key   entry in the ConfigurationBlock dictionary
         * @return      value for key
         */
        template <typename T>
        T get_configuration_scalar(const std::string& key) {
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
        T get_configuration_scalar(const std::string& key, T default_value) {
            try {
                return this->configuration_block_.get_entry<T>(key, default_value);
            } catch (convert::ValueError e) {
                throw this->build_exception(std::string(e.what()) + " (specified \"" + key + "\")");
            }            
        }
        
        /**
         * Retrieves vector of values for specified key.
         * @param key   entry in the ConfigurationBlock dictionary
         * @return      value for key
         */
        template <typename T>
        std::vector<T> get_configuration_vector(const std::string& key) {
            try {
                std::string s = this->configuration_block_.get_entry<std::string>(key);
                return convert::to_vector<T>(s, " ");
            } catch (ConfigurationIncompleteError e) {
                throw this->build_exception(e.what());
            } catch (convert::ValueError e) {
                throw this->build_exception(std::string(e.what()) + " (specified for \"" + key + "\")");
            }            
        }
        
        /**
         * Retrieves vector of cell positions from a configuration entry,
         * parses them, and adds them to the OrganismDistribution object.
         */
        void get_configuration_positions(const std::string& key, OrganismDistribution& od);         
        
        /**
         * Returns vector of keys in entries that start with a particular 
         * string.
         * @param   key_start   string matching beginning of key
         * @return              vector of keys in entries
         */
        std::vector<std::string> get_matching_configuration_keys(const std::string& key_start) const;        
        
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
         */
        WorldConfigurator(const ConfigurationBlock& cb);

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
        /** Number of generations to run. */        
        unsigned long   generations_to_run_;        
        /** Number of fitness factors. */
        unsigned int    num_fitness_factors_;
        /** Random number seed. */        
        unsigned long   rand_seed_;

}; // WorldConfigurator

/**
 * Takes a ConfigurationBlock assumed to be wrapped around Species information, 
 * and parses/translates values appropriately.
 */
class SpeciesConfigurator : public Configurator {

    public:
    
        /** 
         * Constructs objects, and then passes ConfigurationBlock onto parse()
         * for processing. 
         * @param cb                a populated ConfigurationBlock object     
         */
        SpeciesConfigurator(const ConfigurationBlock& cb);

        /** 
         * Takes the string fields of ConfigurationBlock and interprets values
         * as needed for a Species object.
         */        
        void parse();
        
        /**
         * Configures a Species object according to settings.
         */
        void configure(World& world);         
        
    private:
        
        /**
         * Special-case parsing of seed population statements.
         */
        void process_seed_populations();
           
    private:
        /** coefficients for the fitness functions */
        std::vector<float>                  selection_strengths_;
        /** rate of mutation for the genotypic fitness factors */
        float                               mutation_rate_;
        /** window for perturbations of fitness factor values */
        FitnessFactorType                   max_mutation_size_;        
        /** mean number of offspring per female */
        unsigned                            mean_reproductive_rate_;
        /** allowing for evolution in fecundity */
        unsigned                            reproductive_rate_mutation_size_;
        /** landscape migration potential for this species */
        std::vector<int>                    movement_costs_;
        /** movement potential of each organism at the start of each round */
        int                                 movement_capacity_;
        /** genotype for organisms created de novo */
        std::vector<FitnessFactorType>      default_genotypic_fitness_factors_;
        /** seed populations */
        std::vector<OrganismDistribution>   seed_populations_;

}; // SpeciesConfigurator

/**
 * Takes a ConfigurationBlock assumed to be wrapped around Generation
 * information, and parses/translates values appropriately.
 */
class GenerationConfigurator : public Configurator {

    public:
    
        /** 
         * Constructs objects, and then passes ConfigurationBlock onto parse()
         * for processing. 
         * @param cb                a populated ConfigurationBlock object     
         */
        GenerationConfigurator(const ConfigurationBlock& cb);

        /** 
         * Takes the string fields of ConfigurationBlock and interprets values
         * as needed for a Generation object.
         */        
        void parse();
        
        /**
         * Configures a Generation object according to settings.
         */
        void configure(World& world);         
           
    private:
        /** 
         * Environmental regimes that need to be changed/set (expressed as factor
         * indexes mapped to ESRI ASCII Grid file paths). 
         */
        std::map<unsigned, std::string>     environments_;
        
        /** 
         * Movement costs that need to be changed/set. (expressed as species labels
         * mapped to ESRI ASCII Grid file paths). 
         */
        std::map<std::string, std::string>  movement_costs_;
        
        /** 
         * Sampling regime.
         */
        std::map<std::string, SamplingRegime>   samples_;          

}; // GenerationConfigurator

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
         * Parses the configuration file, loading data into blocks.
         */
        void configure(World& world);
        
    private: 
    
        /** Input (file) stream. */
        std::ifstream       fsrc_;
        
        /** Input stream. */
        std::istream&       src_;
        
};  

///////////////////////////////////////////////////////////////////////////////
// helper functions
namespace confsys_detail {
    
    /**
     * Composes and returns and appropriate exception.
     * @param cb                ConfigurationBlock that has the error     
     * @param message           error message
     * @return                  ConfiguratonError exception to be thrown
     */
    ConfigurationError build_configuration_block_exception(const ConfigurationBlock& cb,
            const std::string& message);
    
} // confsys_detail

} // namespace gingko

#endif
