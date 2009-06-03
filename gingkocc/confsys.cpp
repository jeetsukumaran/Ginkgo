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

#include <set>
#include <string>
#include <vector>
#include <ctime>

#include "gingko_defs.hpp"
#include "confsys.hpp"
#include "textutil.hpp"
#include "convert.hpp"
#include "world.hpp"
#include "biosys.hpp"
#include "asciigrid.hpp"
#include "filesys.hpp"

namespace gingko {
namespace confsys {

/* PARSE LOGIC:
 *
 * text => structured text [ConfigurationBlock] => structured values [XXXConf]
 *
 * (1) Data file parsed into [ConfigurationBlock] objects.
 * (2) Each [ConfigurationBlock] object = structured string representation of 
 *     data (e.g., name, type, and dictionary mapping string to strings).
 * (3) Each [ConfigurationBlock] object then is mapped into a corresponding
 *     WorldConfigurator, SpeciesConf, etc. objects, which take the strings and convert
 *     them into values of the appropriate type.
 * (4) The configure_world() functions then take the structured values and 
 *     populate/configure the World object correspondingly.
 */

///////////////////////////////////////////////////////////////////////////////
// Client code should call one of the following to configure World objects.

World& configure_world(World& world, std::istream& conf_src) {
    confsys_detail::ConfigurationFile cf(conf_src);
    cf.configure(world);
    return world;
}

World& configure_world(World& world, const char * conf_fpath) {
    std::string s(conf_fpath);
    return configure_world(world, s);
}

World& configure_world(World& world, const std::string& conf_fpath) {
    std::ifstream f(conf_fpath.c_str());
    confsys_detail::ConfigurationFile cf(conf_fpath);   
    cf.configure(world); 
    return world;
}

///////////////////////////////////////////////////////////////////////////////
// Supporting Classes and Constructs

namespace confsys_detail {

const char * BLOCK_START = "@";
const char * BLOCK_BODY_START = "{";
const char * BLOCK_BODY_END = "}";
const char * BLOCK_BODY_LINE_TERM = ";";
const char * BLOCK_BODY_KEY_VAL_SEP = "=";
const char * WHITESPACE = " \t\n\r";
const char * COMMENT_START = "/*";
const char * COMMENT_END = "*/";

///////////////////////////////////////////////////////////////////////////////
// ConfigurationBlock

// default constructor
ConfigurationBlock::ConfigurationBlock() 
    : is_block_set_(false) { }

// construct and parse
ConfigurationBlock::ConfigurationBlock(std::istream& in)
    : is_block_set_(false) {   
    this->parse(in);
}        

// construct and parse
ConfigurationBlock::ConfigurationBlock(std::istream& in, const std::string& config_filepath)
    : is_block_set_(false),
      config_filepath_(config_filepath) {   
    this->parse(in);
}

// default do-nothing destructor
ConfigurationBlock::~ConfigurationBlock() {}   

// clears/inits
void ConfigurationBlock::clear() {
    this->type_.clear();
    this->name_.clear();
    this->entries_.clear();
    this->is_block_set_ = false;
}

// return type
std::string ConfigurationBlock::get_type() const {
    return this->type_;
}

// return name
std::string ConfigurationBlock::get_name() const {
    return this->name_;
}

// return name
bool ConfigurationBlock::is_block_set() const {
    return this->is_block_set_;
}

// get keys
std::vector<std::string> ConfigurationBlock::get_keys() const {
    std::vector<std::string> keys;
    keys.reserve(this->entries_.size());
    for (std::map< std::string, std::string >::const_iterator e = this->entries_.begin();
            e != this->entries_.end();
            ++e) {
        keys.push_back(e->first);           
    }
    return keys;
}

// get keys matching pattern
std::vector<std::string> ConfigurationBlock::get_keys(const std::string& key_start) const {
    std::vector<std::string> keys;
    keys.reserve(this->entries_.size());
    for (std::map< std::string, std::string >::const_iterator e = this->entries_.begin();
            e != this->entries_.end();
            ++e) {
        if (textutil::startswith(e->first, key_start)) {            
            keys.push_back(e->first);
        }            
    }
    return keys;
}

// wrap up some of the tedium
std::string ConfigurationBlock::compose_error_message(unsigned long pos, const char * desc) {
    std::ostringstream msg;
    msg << "Block starting at character position " <<  pos+1 << ": ";
    msg << desc;
    return msg.str();
 }
 
 // using std string object
std::string ConfigurationBlock::compose_error_message(unsigned long pos, const std::string& desc) {
    return this->compose_error_message(pos, desc.c_str());
}

unsigned long ConfigurationBlock::get_block_start_pos() const {
    return this->block_start_pos_;
}

unsigned long ConfigurationBlock::get_block_end_pos() const {
    return this->block_end_pos_;
}    

std::string ConfigurationBlock::get_config_filepath() const {
    return this->config_filepath_;
}

// workhorse parser
void ConfigurationBlock::parse(std::istream& in) {
    
    this->clear();

    unsigned long start_pos(in.tellg());
    std::string raw = read_block_from_file(in, BLOCK_BODY_END[0]);
    
    if (raw.size() == 0) {
        return;
    }
    
    if (raw[0] != BLOCK_START[0]) {
        std::ostringstream msg;
        msg << "does not begin with '" << BLOCK_START << "'";
        throw ConfigurationSyntaxError(this->compose_error_message(start_pos, msg.str()));
    }
    
    if (in.eof()) {
        std::ostringstream msg;
        msg << "EOF before block body terminator ('" << BLOCK_BODY_END << "')";
        throw ConfigurationSyntaxError(this->compose_error_message(start_pos, msg.str()));                        
    }
    
    std::vector<std::string> parts = textutil::split(raw, BLOCK_BODY_START, false);
    
    if (parts.size() < 2) {
        std::ostringstream msg;
        msg << "missing block body initiator ('" << BLOCK_BODY_START << "')";    
        throw ConfigurationSyntaxError(this->compose_error_message(start_pos, msg.str()));
    }
    
    if (parts.size() > 2) {
        std::ostringstream msg;
        msg << "multiple block body initiators ('" << BLOCK_BODY_START << "')";    
        throw ConfigurationSyntaxError(this->compose_error_message(start_pos, msg.str()));
    }     
    
    std::vector<std::string> head_parts = textutil::split_on_any(textutil::strip(parts[0]), WHITESPACE, false);
    
    if (head_parts.size() > 2) {
        throw ConfigurationSyntaxError(this->compose_error_message(start_pos, "found multiple elements in block header, but expecting only two (type and name)"));
    }

    this->type_ = textutil::lower(textutil::strip(head_parts[0].substr(1)));
    if (head_parts.size() == 2) {
        this->name_ = textutil::strip(head_parts[1]);
    }    
                
    unsigned entry_count = 0;
    std::vector<std::string> body_parts = textutil::split_on_any(textutil::strip(parts[1]), BLOCK_BODY_LINE_TERM, false);
    for (std::vector<std::string>::const_iterator s = body_parts.begin(); s != body_parts.end(); ++s) {
        std::string entry = textutil::strip(*s, WHITESPACE);
        entry_count += 1;
        if (entry.size() > 0) {
            std::vector<std::string> entry_parts = textutil::split(entry, BLOCK_BODY_KEY_VAL_SEP, 1, false);
            if (entry_parts.size() < 2) {
                std::ostringstream msg;
                msg << "incomplete key-value specification in entry no. " << entry_count << " (missing \"=\"): ";
                msg << "\"" << entry << "\"";
                throw ConfigurationSyntaxError(this->compose_error_message(start_pos, msg.str()));
            }
            this->entries_.insert(std::make_pair(textutil::strip(entry_parts[0], WHITESPACE),
                                                 textutil::strip(entry_parts[1], WHITESPACE) ));
        }
    }
    this->is_block_set_ = true;
    this->block_start_pos_ = start_pos;
    this->block_end_pos_ = in.tellg();
}

template <>
std::string ConfigurationBlock::get_entry<std::string>(const std::string& key) const {
    std::map< std::string, std::string >::const_iterator val = this->entries_.find(key);
    if (val == this->entries_.end()) {
        throw ConfigurationIncompleteError("entry \"" + key + "\" not found");
    }
    return val->second;
}

ConfigurationBlock::MultiEntryIterator ConfigurationBlock::equal_range(const std::string& key) {
    return this->entries_.equal_range(key);
}

///////////////////////////////////////////////////////////////////////////////
// ConfigurationBlock inserter

std::istream& operator>> (std::istream& in, ConfigurationBlock& cblock) {
    cblock.parse(in);
    return in;
}

///////////////////////////////////////////////////////////////////////////////
// Configurator

Configurator::Configurator(const ConfigurationBlock& cb) 
        : configuration_block_(cb) { }
          
Configurator::~Configurator() { }

std::vector<std::string> Configurator::get_matching_configuration_keys(const std::string& key_start) const {
    return this->configuration_block_.get_keys(key_start);
}

bool Configurator::parse_position_coordinates(const std::string& pos, CellIndexType& x, CellIndexType &y) {
    std::string pos_clean = textutil::strip(pos, WHITESPACE);
    if (pos_clean.size() == 0) {
        return false;
    }
    std::vector<std::string> pos_parts = textutil::split(pos_clean, ",", 1, false);
    if (pos_parts.size() != 2) {
        throw this->build_exception("invalid or incomplete position specification: \"" + pos_clean + "\"");
    }
    try {
        x = convert::to_scalar<CellIndexType>(pos_parts[0]);
    } catch (const convert::ValueError& e) {
        throw this->build_exception("invalid x-coordinate value for position: \"" + pos_parts[0] + "\"");
    }
    try {
        y = convert::to_scalar<CellIndexType>(pos_parts[1]);
    } catch (const convert::ValueError& e) {
        throw this->build_exception("invalid y-coordinate value for position: \"" + pos_parts[1] + "\"");
    }
    return true;
}

void Configurator::get_configuration_positions(const std::string& key, OrganismDistribution& od, bool allow_wildcard) {
    std::vector<std::string> positions = this->get_configuration_vector<std::string>(key);
    
    // "*" = sample all cells
    // which is indicated by returning an empty list
    if (allow_wildcard and positions.size() == 1 and positions[0] == "*") {
        return;
    }
    
    od.x.reserve(od.x.size() + positions.size());
    od.y.reserve(od.y.size() + positions.size());
    CellIndexType x = 0;
    CellIndexType y = 0;
    for (std::vector<std::string>::const_iterator pos = positions.begin(); pos < positions.end(); ++pos) {
        if (this->parse_position_coordinates(*pos, x, y)) {
            od.x.push_back(x);
            od.y.push_back(y);
        }
    }
}

ConfigurationError Configurator::build_exception(const std::string& message) const {
    return build_configuration_block_exception(this->configuration_block_, message);
}

CellIndexType Configurator::xy_to_index(World& world, CellIndexType x, CellIndexType y) {
    if (x > world.landscape().size_x()-1) {
        std::ostringstream msg;
        msg << "maximum x-coordinate on landscape is " << world.landscape().size_x();
        msg << " but position specifies x-coordinate of " << x;
        throw this->build_exception(msg.str());
    }
    if (y > world.landscape().size_y()-1) {
        std::ostringstream msg;
        msg << "maximum x-coordinate on landscape is " << world.landscape().size_y();
        msg << " but position specifies x-coordinate of " << y;
        throw this->build_exception(msg.str());
    } 
    return world.landscape().xy_to_index(x, y);
}

///////////////////////////////////////////////////////////////////////////////
// WorldConfigurator

WorldConfigurator::WorldConfigurator(const ConfigurationBlock& cb) 
        : Configurator(cb),
          size_x_(0),
          size_y_(0),
          num_fitness_factors_(0),
          rand_seed_(0),
          produce_final_output_(true) {
    this->parse();
}

void WorldConfigurator::parse()  {
    this->size_x_ = this->get_configuration_scalar<CellIndexType>("ncols"); 
    this->size_y_ = this->get_configuration_scalar<CellIndexType>("nrows");
    this->generations_to_run_ = this->get_configuration_scalar<unsigned long>("ngens");
    this->num_fitness_factors_ = this->get_configuration_scalar<unsigned>("nfitness", MAX_FITNESS_FACTORS);
    if (this->num_fitness_factors_ > MAX_FITNESS_FACTORS) {
        std::ostringstream s;
        s << "maximum number of fitness factors allowed is " << MAX_FITNESS_FACTORS;        
        s << ", but requested " << this->num_fitness_factors_;
        throw this->build_exception(s.str());
    }
    this->fitness_factor_grain_ = this->get_configuration_scalar<unsigned>("fitness-grain", 1);
    if (this->fitness_factor_grain_ == 0) {
        throw this->build_exception("fitness factor grain must be > 0");  
    }
    try {
        this->rand_seed_ = this->get_configuration_scalar<unsigned>("rseed");
    } catch (const ConfigurationError& e) {
        // rely on rng constructor to take care of this
        // this->rand_seed_ = ctime();
        this->rand_seed_ = time(0);
    }
    this->produce_final_output_ = not this->get_configuration_scalar<bool>("suppress-final-output", 0);
}

void WorldConfigurator::configure(World& world)  {
    world.set_label(this->get_name("GingkoWorld"));        
    world.set_random_seed(this->rand_seed_);
    world.set_generations_to_run(this->generations_to_run_);
    world.set_num_fitness_factors(this->num_fitness_factors_);
    world.set_fitness_factor_grain(this->fitness_factor_grain_);
    world.generate_landscape(this->size_x_, this->size_y_);
    world.set_produce_final_output(this->produce_final_output_);
}

///////////////////////////////////////////////////////////////////////////////
// SpeciesConfigurator

SpeciesConfigurator::SpeciesConfigurator(const ConfigurationBlock& cb) 
        : Configurator(cb),
          mutation_rate_(0),
          max_mutation_size_(0),
          mean_reproductive_rate_(0),
          reproductive_rate_mutation_size_(0),
          movement_capacity_(0) {
    this->parse();
}

void SpeciesConfigurator::parse()  {
    try {
        this->selection_weights_ = this->get_configuration_vector<float>("selection-weights");
    } catch (ConfigurationIncompleteError& e) {
        this->selection_weights_.assign(MAX_FITNESS_FACTORS, 1);
    }
    try {
        this->default_genotypic_fitness_factors_ = this->get_configuration_vector<FitnessFactorType>("genotypic-fitness");
    } catch (ConfigurationIncompleteError& e) {
        this->default_genotypic_fitness_factors_.assign(MAX_FITNESS_FACTORS, 0);
    }
    
    this->mutation_rate_ = this->get_configuration_scalar<float>("mutation-rate", 0.01);
    this->max_mutation_size_ = this->get_configuration_scalar<unsigned>("max-mutation-size", 1);
    this->mean_reproductive_rate_ = this->get_configuration_scalar<unsigned>("fecundity", 8);    
    this->reproductive_rate_mutation_size_ = this->get_configuration_scalar<unsigned>("fecundity-evolution-size", 1);
    this->movement_capacity_ = this->get_configuration_scalar<unsigned>("movement-capacity", 1);
    this->process_seed_populations();
}

void SpeciesConfigurator::configure(World& world)  {
    std::string label = this->get_name();
    if (world.has_species(label)) {
        throw this->build_exception("species \"" + label + "\" has already been defined");
    }
    Species& sp = world.new_species(label);
    if (sp.get_num_fitness_factors() != this->selection_weights_.size()) {
        std::ostringstream msg;
        msg << "expecting " << sp.get_num_fitness_factors() << " default selection weights factors, but found " << this->selection_weights_.size() ;
        throw this->build_exception(msg.str());
    }    
    if (sp.get_num_fitness_factors() != this->default_genotypic_fitness_factors_.size()) {
        std::ostringstream msg;
        msg << "expecting " << sp.get_num_fitness_factors() << " default genotypic fitness factors, but found " << this->default_genotypic_fitness_factors_.size() ;
        throw this->build_exception(msg.str());
    }
    sp.set_default_genotypic_fitness_factors(this->default_genotypic_fitness_factors_);
    sp.set_mutation_rate(this->mutation_rate_);
    sp.set_max_mutation_size(this->max_mutation_size_);
    sp.set_mean_reproductive_rate(this->mean_reproductive_rate_);
//     sp.set_reproductive_rate_mutation_size(this->reproductive_rate_mutation_size_);
    sp.set_movement_capacity(this->movement_capacity_);
    for (std::vector<OrganismDistribution>::iterator odi = this->seed_populations_.begin(); odi != this->seed_populations_.end(); ++odi) {
        OrganismDistribution& od = *odi;
        assert(od.x.size() == od.y.size());
        for (unsigned i = 0; i < od.x.size(); ++i) {
            if (od.x[i] > world.landscape().size_x()-1) {
                std::ostringstream msg;
                msg << "maximum x-coordinate on landscape is " << world.landscape().size_x();
                msg << " but position specifies x-coordinate of " << od.x[i];
                throw this->build_exception(msg.str());
            }
            if (od.y[i] > world.landscape().size_y()-1) {
                std::ostringstream msg;
                msg << "maximum x-coordinate on landscape is " << world.landscape().size_y();
                msg << " but position specifies x-coordinate of " << od.y[i];
                throw this->build_exception(msg.str());
            }              
            CellIndexType cell_index = world.landscape().xy_to_index(od.x[i], od.y[i]);
            world.add_seed_population(cell_index, label, od.num_organisms_per_cell, od.ancestral_population_size, od.ancestral_generations);
        }
    }
}

void SpeciesConfigurator::process_seed_populations() {
    std::vector<std::string> keys = this->get_matching_configuration_keys("init");
    for (std::vector<std::string>::const_iterator k = keys.begin(); k != keys.end(); ++k) {
        this->seed_populations_.push_back(OrganismDistribution());
        OrganismDistribution& od = this->seed_populations_.back();        
        std::string key = *k;
        std::vector<std::string> key_parts = textutil::split(key, "(", 1, false);
        if (key_parts.size() < 2) {
            throw this->build_exception("full specification of seed population in the form of \"init(n:N*G)\" required, where: n = seed population size, N = ancestral seed population size, and G = ancestral bootstrap number of generations.");
        }
        std::vector<std::string> subparts = textutil::split_on_any(key_parts[1], ":*", 0, false);
        if (subparts.size() != 3) {
            throw this->build_exception("full specification of seed population in the form of \"init(n:N*G)\" required, where: n = seed population size, N = ancestral seed population size, and G = ancestral bootstrap number of generations.");
        }
        try {
            od.num_organisms_per_cell = convert::to_scalar<unsigned long>(subparts[0]);
        } catch (const convert::ValueError& e) {
            throw this->build_exception("invalid value for seed population size: \"" + subparts[0] + "\"");
        }
        try {
            od.ancestral_population_size = convert::to_scalar<unsigned long>(subparts[1]);
        } catch (const convert::ValueError& e) {
            throw this->build_exception("invalid value for seed population ancestral size: \"" + subparts[1] + "\"");
        }
        try {
            std::string g = textutil::strip(subparts[1], ")");
            od.ancestral_generations = convert::to_scalar<unsigned long>(g);
        } catch (const convert::ValueError& e) {
            throw this->build_exception("invalid value for seed population ancestral generations : \"" + subparts[2] + "\"");
        }           
        this->get_configuration_positions(key, od);
    }
}

///////////////////////////////////////////////////////////////////////////////
// GenerationConfigurator

GenerationConfigurator::GenerationConfigurator(const ConfigurationBlock& cb) 
        : Configurator(cb),
          generation_(0) {
    this->parse();
}

std::string GenerationConfigurator::get_validated_grid_path(const std::string& grid_path, const World& world) {   
    std::string root_filepath = filesys::get_path_parent(this->get_config_filepath());
    std::string full_grid_path;    
    if (filesys::is_abs_path(grid_path) or root_filepath.size() == 0) {
        full_grid_path = grid_path;      
    } else {        
        full_grid_path = filesys::compose_path(root_filepath, grid_path);      
    } 
    try {
        asciigrid::AsciiGrid grid(full_grid_path);
        std::vector<long> values = grid.get_cell_values();
        if (values.size() != world.size()) {
            std::ostringstream msg;
            msg << "landscape has " << world.size() << " cells, ";
            msg << "but grid \"" << full_grid_path << "\" describes " << values.size() << " cells";
            throw this->build_exception(msg.str());        
        }
        return full_grid_path;
    } catch (asciigrid::AsciiGridIOError e) {
        throw this->build_exception("I/O error reading grid \"" + full_grid_path + "\": " + e.what());
    } catch (asciigrid::AsciiGridFormatError e) {
        throw this->build_exception("format error reading grid \"" + full_grid_path + "\": " + e.what());
    }
}

void GenerationConfigurator::process_carrying_capacity() {
    this->carrying_capacity_ = this->get_configuration_scalar<std::string>("carrying-capacity", "");
}

void GenerationConfigurator::process_environments() {
    std::vector<std::string> keys = this->get_matching_configuration_keys("environment");
    for (std::vector<std::string>::iterator k = keys.begin(); k != keys.end(); ++k) {
        std::string& key = *k;
        std::vector<std::string> key_parts = textutil::split(key, ":", 1, false);
        if (key_parts.size() < 2) {
            throw this->build_exception("need to specify environment factor index");
        }
        unsigned environment_factor = 0;
        try {
            environment_factor = convert::to_scalar<unsigned long>(key_parts[1]);
        } catch (const convert::ValueError& e) {
            throw this->build_exception("invalid value for environment factor index: \"" + key_parts[1] + "\"");
        }
        if (environment_factor == 0) {
            throw this->build_exception("invalid value for environment factor index: 0 (environment factor indexes start at 1)"); 
        }
        this->environments_.insert(std::make_pair(environment_factor, this->get_configuration_scalar<std::string>(key)));
    }
}

void GenerationConfigurator::process_movement_costs() {
    std::vector<std::string> keys = this->get_matching_configuration_keys("movement");
    for (std::vector<std::string>::iterator k = keys.begin(); k != keys.end(); ++k) {
        std::string& key = *k;
        std::vector<std::string> key_parts = textutil::split(key, ":", 1, false);
        if (key_parts.size() < 2) {
            throw this->build_exception("need to specify species label for movement costs");
        }
        this->movement_costs_.insert(std::make_pair(key_parts[1], this->get_configuration_scalar<std::string>(key)));
    }
}

void GenerationConfigurator::process_dispersals() {
    std::vector<std::string> keys = this->get_matching_configuration_keys("disperse");
    for (std::vector<std::string>::iterator k = keys.begin(); k != keys.end(); ++k) {
        OrganismDispersal od;
        std::string& key = *k;
        std::vector<std::string> key_parts = textutil::split(key, ":", 1, false);
        if (key_parts.size() < 2) {
            throw this->build_exception("must specify species using \":\" token for dispersal event");
        }
        std::vector<std::string> species_num_parts = textutil::split(key_parts[1], "#", 1, false);
        od.species_label = species_num_parts[0];
        if (species_num_parts.size() == 2) {
            try {
                od.num_organisms = convert::to_scalar<unsigned long>(species_num_parts[1]);
            } catch (const convert::ValueError& e) {
                throw this->build_exception("invalid value for number of organisms: \"" + species_num_parts[1] + "\"");
            }        
        } else {
            od.num_organisms = 0;        
        }
        std::string disp_pos = this->get_configuration_scalar<std::string>(key);
        std::vector<std::string> disp_pos_parts = textutil::split(disp_pos, "|", 1, false);
        if (disp_pos_parts.size() == 1) {
            od.probability = 1.0;
        } else {
            try {
                od.probability = convert::to_scalar<float>(disp_pos_parts[1]);
            } catch (const convert::ValueError& e) {
                throw this->build_exception("invalid value for probability of dispersal: \"" + disp_pos_parts[1] + "\"");
            }             
        }
        std::vector<std::string> positions = textutil::split_on_any(disp_pos_parts[0], WHITESPACE, 1, false);
        if (positions.size() < 2) {
            throw this->build_exception("must specify source and destination positions for dispersal");
        }
        this->parse_position_coordinates(positions[0], od.src_x, od.src_y);
        this->parse_position_coordinates(positions[1], od.dest_x, od.dest_y);
        this->dispersals_.push_back(od);
    }
}

void GenerationConfigurator::parse()  {
    try {
        this->generation_ = convert::to_scalar<unsigned long>(this->get_name());
    } catch (convert::ValueError e) {
        throw this->build_exception("\"" + this->get_name() + "\" is an invalid value for a generation number");
    }  
    this->process_carrying_capacity();
    this->process_environments();
    this->process_movement_costs();
    this->process_dispersals();
}

void GenerationConfigurator::configure(World& world)  {
    WorldSettings world_settings;

    if (this->carrying_capacity_.size() > 0) {
        world_settings.carrying_capacity = this->get_validated_grid_path(this->carrying_capacity_, world);
    }
    for (std::map<unsigned, std::string>::iterator envi = this->environments_.begin(); envi != this->environments_.end(); ++envi) {
        if (envi->first > world.get_num_fitness_factors()) {
            std::ostringstream msg;
            msg << "invalid environment factor index: " << envi->first;
            msg << " (maximum valid index is " << world.get_num_fitness_factors();
            msg << ", given " << world.get_num_fitness_factors() << " defined factors)";
            throw this->build_exception(msg.str());           
        }
        world_settings.environments.insert(std::make_pair(envi->first-1, this->get_validated_grid_path(envi->second, world)));
    }
    for (std::map<std::string, std::string>::iterator mci = this->movement_costs_.begin(); mci != this->movement_costs_.end(); ++mci) {
        if (not world.has_species(mci->first)) {
            throw this->build_exception("movement costs: species \"" + mci->first + "\" not defined");
        }
        world_settings.movement_costs.insert(std::make_pair(mci->first, this->get_validated_grid_path(mci->second, world)));
    }
    for (std::vector<OrganismDispersal>::iterator odi = this->dispersals_.begin(); odi != this->dispersals_.end(); ++odi) {
        OrganismDispersal& od = *odi;
        DispersalEvent de;        
        if (world.has_species(od.species_label)) {
            de.species_ptr = world.get_species_ptr(od.species_label);
        } else {
            throw this->build_exception("Species \"" + od.species_label + "\" not defined");
        }
        de.num_organisms = od.num_organisms;
        de.source = this->xy_to_index(world, od.src_x, od.src_y);
        de.destination = this->xy_to_index(world, od.dest_x, od.dest_y);
        de.probability = od.probability;
        world_settings.dispersal_events.push_back(de);
    }
    world.add_world_settings(this->generation_, world_settings);
}

///////////////////////////////////////////////////////////////////////////////
// TreeSamplingConfigurator

TreeSamplingConfigurator::TreeSamplingConfigurator(const ConfigurationBlock& cb) 
        : Configurator(cb),
          random_sample_size_(0) {
    this->parse();
}

void TreeSamplingConfigurator::parse() {
    this->generation_ = this->get_configuration_scalar<unsigned long>("gen");
    this->organism_sampling_.species_label = this->get_configuration_scalar<std::string>("species");
    this->organism_sampling_.num_organisms_per_cell = this->get_configuration_scalar<unsigned long>("limit-per-cell", 0);    
    if (this->has_configuration_entry("cells")) {
        this->get_configuration_positions("cells", this->organism_sampling_);
    } else if (this->has_configuration_entry("random")) {
        // not yet implemented
    }
}

void TreeSamplingConfigurator::configure(World& world) {
    SamplingRegime world_sampling_regime;
    
    if (world.has_species(this->organism_sampling_.species_label)) {
        world_sampling_regime.species_ptr = world.get_species_ptr(this->organism_sampling_.species_label);
    } else {
        throw this->build_exception("Species \"" + this->organism_sampling_.species_label + "\" not defined");
    } 
    
    if (this->has_name()) {
        world_sampling_regime.label = this->get_name();
    }    
    
    world_sampling_regime.num_organisms_per_cell = this->organism_sampling_.num_organisms_per_cell;
    
    assert(this->organism_sampling_.x.size() == this->organism_sampling_.y.size());
    if (this->organism_sampling_.x.size() > 0) {
        for (unsigned i = 0; i < this->organism_sampling_.x.size(); ++i) {             
            world_sampling_regime.cell_indexes.insert(this->xy_to_index(world, this->organism_sampling_.x[i], this->organism_sampling_.y[i]));
        }
            
    } else if (this->random_sample_size_ > 0) {
        // not implemented
    } else {
        for (CellIndexType i = 0; i < world.landscape().size(); ++i) {
            world_sampling_regime.cell_indexes.insert(i);
        }
    }
    world.add_tree_sampling(this->generation_, world_sampling_regime);
}

///////////////////////////////////////////////////////////////////////////////
// OccurrenceSamplingConfigurator

OccurrenceSamplingConfigurator::OccurrenceSamplingConfigurator(const ConfigurationBlock& cb) 
        : Configurator(cb) {
    this->parse();
}

void OccurrenceSamplingConfigurator::parse() {
    this->generation_ = this->get_configuration_scalar<unsigned long>("gen");
    this->species_label_ = this->get_configuration_scalar<std::string>("species");
}

void OccurrenceSamplingConfigurator::configure(World& world) {
    Species * species_ptr;
    if (world.has_species(this->species_label_)) {
        species_ptr = world.get_species_ptr(this->species_label_);
    } else {
        throw this->build_exception("Species \"" + this->species_label_ + "\" not defined");
    }     
    world.add_occurrence_sampling(this->generation_, species_ptr);
}

///////////////////////////////////////////////////////////////////////////////
// ConfigurationFile

ConfigurationFile::ConfigurationFile(std::istream& src)
        : src_(src) {
    if (not this->src_) {
        throw ConfigurationIOError("invalid source stream");
    }
}

ConfigurationFile::ConfigurationFile(const char * fpath)
        : config_filepath_(fpath),
          fsrc_(fpath),
          src_(fsrc_) {
    if (not this->src_) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw ConfigurationIOError(msg.str());
    }          
}

ConfigurationFile::ConfigurationFile(const std::string& fpath) 
        : config_filepath_(fpath),
          fsrc_(fpath.c_str()),
          src_(fsrc_) {
    if (not this->src_) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw ConfigurationIOError(msg.str());
    }
}

ConfigurationFile::~ConfigurationFile() { }

void ConfigurationFile::configure(World& world) {
    assert(this->src_);
    unsigned num_worlds = 0;
    std::set<std::string> species_labels;
    std::set<unsigned long> generations;
    while (not this->src_.eof()) {        
        ConfigurationBlock cb(this->src_, this->config_filepath_);
//         this->src_ >> cb;
        if (cb.is_block_set()) {
            if (cb.get_type() == "world") {
                if (num_worlds != 0) {
                    throw build_configuration_block_exception(cb, "world must be defined before species are added");
                }                
                WorldConfigurator wcf(cb);
                wcf.configure(world);
                num_worlds += 1;
            }
            if (cb.get_type() == "species") {
                if (num_worlds == 0) {
                    throw build_configuration_block_exception(cb, "world must be defined before species are added");
                }
                std::string label = cb.get_name();
                if (species_labels.find(label) != species_labels.end()) {
                    throw build_configuration_block_exception(cb, "species \"" + label + "\" has already been defined");
                }
                SpeciesConfigurator spcf(cb);
                spcf.configure(world);
                species_labels.insert(label);
            }
            if (cb.get_type() == "generation") {
                if (num_worlds == 0) {
                    throw build_configuration_block_exception(cb, "world must be defined before world setting changes defined");
                }
                unsigned long generation = 0;
                try {
                    generation = convert::to_scalar<unsigned long>(cb.get_name());
                } catch (convert::ValueError e) {
                    throw build_configuration_block_exception(cb, "\"" + cb.get_name() + "\" is an invalid value for a generation number");
                }                
                if (generations.find(generation) != generations.end()) {
                    throw build_configuration_block_exception(cb, "generation \"" + cb.get_name() + "\" has already been defined");
                }
                GenerationConfigurator gcf(cb);
                gcf.configure(world);
            }
            if (cb.get_type() == "tree") {
                if (num_worlds == 0) {
                    throw build_configuration_block_exception(cb, "world must be defined before tree building directives defined");
                }
                TreeSamplingConfigurator scf(cb);
                scf.configure(world);
            }
            if (cb.get_type() == "occurs") {
                if (num_worlds == 0) {
                    throw build_configuration_block_exception(cb, "world must be defined before occurrence sampling directives defined");
                }
                OccurrenceSamplingConfigurator ocf(cb);
                ocf.configure(world);
            }              
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Implementation details
    
/**
 * Composes and returns and appropriate exception.
 * @param message           error message
 * @param cb                ConfigurationBlock that has the error
 * @return                  ConfiguratonError exception to be thrown
 */
ConfigurationError build_configuration_block_exception(const ConfigurationBlock& cb,
        const std::string& message) {
    std::ostringstream msg;
    msg << cb.get_type() << " block ";
    if (cb.has_name()) {        
        msg << "\"" << cb.get_name() << "\" ";
    }            
    msg << "(file position " << cb.get_block_start_pos() + 1;
    msg << " to position " << cb.get_block_end_pos() + 1 << ")";
    msg << ": " << message;
    return ConfigurationError(msg.str());    
}

std::string read_block_from_file(std::istream& is, char block_terminator) {
    std::ostringstream block;
    unsigned int comment_level = 0;
    char c = is.get();
    while (is.good() and (comment_level > 0 or c != block_terminator)) {
        if (comment_level == 0) {
            if ( c == COMMENT_START[0] ) {
                char c2 = is.get();
                if (is.good() and c2 == COMMENT_START[1]) {
                    comment_level = 1;
                } else {
                    block << c << c2;
                }
            } else {        
                block << c;
            }                
        } else {
            if ( c == COMMENT_END[0] ) {
                c = is.get();
                if (is.good() and c == COMMENT_END[1]) {
                    comment_level = 0;
                }
            }
        }
        c = is.get();
    }
    return textutil::strip(block.str());
}
    
} // confsys_detail

} // namespace confsys

} // namespace gingko
