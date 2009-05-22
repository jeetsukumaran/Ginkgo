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
namespace confsys {

///////////////////////////////////////////////////////////////////////////////
// Client code should call one of the following to configure World objects.

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
// ERRORS AND EXCEPTIONS

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

///////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION DETAILS

namespace confsys_detail {

/**
 * Track distribution of organisms, either for seed populations or sampling 
 * regimes.
 */
struct OrganismDistribution {
    
    public:
        OrganismDistribution()
            : num_organisms_per_cell(0) { }

    public:
        std::string     species_label;
        unsigned long   num_organisms_per_cell;
        std::vector<CellIndexType>  x;
        std::vector<CellIndexType>  y;
        

}; // OrganismDistribution

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
         * Constructor that populates a ConfigurationBlock object by calling
         * parse(), also storing filepath for path manipulation.
         * @param in                input stream with data
         * @param config_filepath   path of configuration file         
         */
        ConfigurationBlock(std::istream& in, const std::string& config_filepath);            
        
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
         * Returns exactly as with multimap.equal_range.
         */
        typedef std::pair<std::multimap<std::string, std::string>::iterator, std::multimap<std::string, std::string>::iterator> MultiEntryIterator;
        MultiEntryIterator equal_range(const std::string& key);
        
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
         * Returns <code>true</code> if name is set.
         */
        bool has_name() const {
            return (this->name_.size() > 0);
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
        
        /**
         * Returns path of the source configuration file, or an empty string if
         * this is not defined.
         * @return      path of source file (if defined) or empty string
         */
        std::string get_config_filepath() const;

    private:
        /** The type of block (e.g. "species", "world", "generation") */
        std::string                             type_;
        /** The name of the block (e.g. "Sp1"). */
        std::string                             name_;
        /** Key-value pairs making up the block body. */
        std::multimap< std::string, std::string >    entries_;
        /** Tracks whether or not the block was actually set. */
        bool                                    is_block_set_;
        /** Start of block in the source stream. */
        unsigned long                           block_start_pos_;
        /** End of block in the source stream. */
        unsigned long                           block_end_pos_;
        /** Filepath source of configuration file. */
        std::string                             config_filepath_;
        
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
         */
        Configurator(const ConfigurationBlock& cb);
    
        /** No-op destructor. */                     
        virtual ~Configurator();                     
                     
        /** 
         * Takes the string fields of ConfigurationBlock and interprets values
         * as needed for a World object.   
         */        
        virtual void parse() = 0;
        
        /**
         * Returns type of underlying ConfigurationBlock.
         * @return  type of block
         */
        std::string get_type() {
            return this->configuration_block_.get_type();
        }        
        
        /**
         * Returns name of underlying ConfigurationBlock.
         * @return  name of block
         */
        std::string get_name() {
            std::string name = this->configuration_block_.get_name();
            if (name.size() == 0) {
                throw this->build_exception("name required for block but not specified");
            }
            return name;
        }
        
        /**
         * Returns name of underlying ConfigurationBlock.
         * @return  name of block
         */
        std::string get_name(const std::string& default_name) {
            std::string name = this->configuration_block_.get_name();
            if (name.size() == 0) {
                return default_name;
            } else {
                return name;
            }
        }
        
        /**
         * Returns <code>true</code> if name is set.
         */
        bool has_name() const {
            return this->configuration_block_.has_name();
        }            
        
        /**
         * Retrieves value for specified key.
         * @param key   entry in the ConfigurationBlock dictionary
         * @return      value for key
         */
        template <typename T>
        T get_configuration_scalar(const std::string& key) const {
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
        T get_configuration_scalar(const std::string& key, T default_value) const {
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
        std::vector<T> get_configuration_vector(const std::string& key) const {
            try {
                std::string s = this->configuration_block_.get_entry<std::string>(key);
                std::vector<T> x = convert::to_vector<T>(s, " ", true);
                return x;
            } catch (ConfigurationIncompleteError e) {
                throw this->build_exception(e.what());
            } catch (convert::ValueError e) {
                throw this->build_exception(std::string(e.what()) + " (specified for \"" + key + "\")");
            }            
        }
        
        /**
         * Returns <code>true</code> if entry exists.
         */
        bool has_configuration_entry(const std::string& key) const {
            return this->configuration_block_.has_key(key);
        }
        
        /** Converts x,y to cell index, with validation. */
        CellIndexType xy_to_index(World& world, CellIndexType x, CellIndexType y);
        
        /**
         * Parses string in the form of "x,y" to a pair of coordinates.
         */
        bool parse_position_coordinates(const std::string& pos, CellIndexType& x, CellIndexType &y);
        
        /**
         * Retrieves vector of cell positions from a configuration entry,
         * parses them, and adds them to the OrganismDistribution object.         
         * @param key               entry key
         * @param od                OrganismDistribution object to populate
         * @param allow_wildcards   allow '*' to stand in for all positions
         */
        void get_configuration_positions(const std::string& key, OrganismDistribution& od, bool allow_wildcard=false);         
        
        /**
         * Returns vector of keys in entries that start with a particular 
         * string.
         * @param   key_start   string matching beginning of key
         * @return              vector of keys in entries
         */
        std::vector<std::string> get_matching_configuration_keys(const std::string& key_start) const;        
        
        /**
         * Returns path of the source configuration file, or an empty string if
         * this is not defined.
         * @return      path of source file (if defined) or empty string
         */
        std::string get_config_filepath() const {
            return this->configuration_block_.get_config_filepath();
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
        /** Fitness factor scaling */
        unsigned int    fitness_factor_grain_;
        /** Random number seed. */        
        unsigned long   rand_seed_;
        /** Produce final output? */
        bool            produce_final_output_;        

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
        std::vector<float>                  selection_weights_;
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
         * as needed for a WorldSettings object.
         */        
        void parse();
        
        /**
         * Configures a WorldSettings according to specs.
         */
        void configure(World& world);
        
    private:
    
        /** 
         * Verifies a grid file, ensuring that it exists, can be loaded
         * and is the same dimensions as a given World, and, if successful
         * returns the validated filepath. Throws an exception if any of these 
         * conditions fail.
         * @param   grid_path   path to grid file
         * @param   world       world against which to validate
         * @return              validated file path
         */
        std::string get_validated_grid_path(const std::string& grid_path, const World& world);        
        
        /** Processes carrying capacity entries. */
        void process_carrying_capacity();
        
        /** Processes environment entries. */
        void process_environments();
        
        /** Processes movement costs. */
        void process_movement_costs();
        
        /** Processes dispersal specs */
        void process_dispersals();
           
    private:
    
        struct OrganismDispersal {
        
            public:
                OrganismDispersal() 
                    : num_organisms(0),
                      src_x(0),
                      src_y(0),
                      dest_x(0),
                      dest_y(0),
                      probability(0) { }
        
            public:
                /** Pointer to species. */ 
                std::string                 species_label;       
                /**  Number of organisms from cell to be sampled (0 = all). */
                unsigned long               num_organisms;
                /** Origin cell index. */
                CellIndexType               src_x;
                CellIndexType               src_y;
                /** Destination cell index. */
                CellIndexType               dest_x;
                CellIndexType               dest_y;                
                /** Probability of dispersal. */
                float                       probability;
                
        };    
    
    
        /** Generation #. */
        unsigned long                       generation_;
    
        /** Path to grid defining the carrying capacity. */
        std::string                         carrying_capacity_;
        
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
        std::vector<OrganismDispersal>      dispersals_;     

}; // GenerationConfigurator

/**
 * Takes a ConfigurationBlock assumed to be wrapped around a Tree sampling
 * directive block information, and parses/translates values appropriately.
 */
class TreeSamplingConfigurator : public Configurator {

    public:
    
        /** 
         * Constructs objects, and then passes ConfigurationBlock onto parse()
         * for processing. 
         * @param cb                a populated ConfigurationBlock object     
         */
        TreeSamplingConfigurator(const ConfigurationBlock& cb);

        /** 
         * Takes the string fields of ConfigurationBlock and interprets values
         * as needed for a Tree/Occurrence samples.
         */        
        void parse();
        
        /**
         * Configures a the sampling regime.
         */
        void configure(World& world);
        
    private:
        /** Generation to sample. */
        unsigned long           generation_;
        /** Sampling regime if custom positions are set. */
        OrganismDistribution    organism_sampling_;
        /** Random sampling regime. */
        unsigned long           random_sample_size_;
};

/**
 * Takes a ConfigurationBlock assumed to be wrapped around a Tree sampling
 * directive block information, and parses/translates values appropriately.
 */
class OccurrenceSamplingConfigurator : public Configurator {

    public:
    
        /** 
         * Constructs objects, and then passes ConfigurationBlock onto parse()
         * for processing. 
         * @param cb                a populated ConfigurationBlock object     
         */
        OccurrenceSamplingConfigurator(const ConfigurationBlock& cb);

        /** 
         * Takes the string fields of ConfigurationBlock and interprets values
         * as needed for a Tree/Occurrence samples.
         */        
        void parse();
        
        /**
         * Configures a the sampling regime.
         */
        void configure(World& world);
        
    private:
        /** Generation to sample. */
        unsigned long           generation_;
        /** Species label. */
        std::string             species_label_;
};

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
        
        /** Path to configuration file. */
        std::string         config_filepath_;    
    
        /** Input (file) stream. */
        std::ifstream       fsrc_;
        
        /** Input stream. */
        std::istream&       src_;
        
};  

///////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS
    
/**
 * Composes and returns and appropriate exception.
 * @param cb                ConfigurationBlock that has the error     
 * @param message           error message
 * @return                  ConfiguratonError exception to be thrown
 */
ConfigurationError build_configuration_block_exception(const ConfigurationBlock& cb,
        const std::string& message);

/**
 * Returns the next chunk of characters from the current file position to the 
 * first occurrence of the block terminator, skipping over characters commented
 * out.
 * @param   in                  input stream
 * @param   block_terminator    token signifying end of block
 * @return                      block
 */
std::string read_block_from_file(std::istream& is, char block_terminator);
    
} // confsys_detail
} // namespace confsys
} // namespace gingko

#endif
