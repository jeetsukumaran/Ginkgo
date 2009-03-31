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

#include "gingko_defs.hpp"
#include "confsys.hpp"
#include "textutil.hpp"
#include "convert.hpp"

namespace gingko {

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
    ConfigurationFile cf(conf_src);
    return world;
}

World& configure_world(World& world, const char * conf_fpath) {
    std::ifstream f(conf_fpath);
    return configure_world(world, f);
}

World& configure_world(World& world, const std::string& conf_fpath) {
    std::ifstream f(conf_fpath.c_str());
    return configure_world(world, f);
}

///////////////////////////////////////////////////////////////////////////////
// Supporting Classes and Constructs

const char * BLOCK_START = "@";
const char * BLOCK_BODY_START = "{";
const char * BLOCK_BODY_END = "}";
const char * BLOCK_BODY_LINE_TERM = ";\n";
const char * BLOCK_BODY_KEY_VAL_SEP = "=";
const char * WHITESPACE = " \t\n";

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

// get entry values by keys
// std::string ConfigurationBlock::get_entry(const std::string& key) const {
//     std::map< std::string, std::string >::const_iterator val = this->entries_.find(key);
//     assert(val != this->entries_.end());
//     return val->second;
// }

// get keys
std::vector<std::string> ConfigurationBlock::get_keys() const {
    std::vector<std::string> keys;
    keys.reserve(this->entries_.size());
    for (std::map< std::string, std::string >::const_iterator e = this->entries_.begin();
            e != this->entries_.end();
            ++e) {
        keys.push_back(e->second);           
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

// workhorse parser
void ConfigurationBlock::parse(std::istream& in) {
    
    this->clear();

    std::string raw;
    unsigned long start_pos(in.tellg());
    std::getline(in, raw, BLOCK_BODY_END[0]);
    raw = textutil::strip(raw);
    
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
    
    if (head_parts.size() < 2) {
        throw ConfigurationSyntaxError(this->compose_error_message(start_pos, "found only one element in block header, but expecting two (type and name)"));
    }
    
    if (head_parts.size() > 2) {
        throw ConfigurationSyntaxError(this->compose_error_message(start_pos, "found multiple elements in block header, but expecting only two (type and name)"));
    }
    
    this->type_ = textutil::lower(textutil::strip(head_parts[0].substr(1)));
    this->name_ = textutil::strip(head_parts[1]);
    
    unsigned entry_count = 0;
    std::vector<std::string> body_parts = textutil::split_on_any(textutil::strip(parts[1]), BLOCK_BODY_LINE_TERM, false);
    for (std::vector<std::string>::const_iterator s = body_parts.begin(); s != body_parts.end(); ++s) {
        std::string entry = textutil::strip(*s, WHITESPACE);
        entry_count += 1;
        if (entry.size() > 0) {
            std::vector<std::string> entry_parts = textutil::split(entry, BLOCK_BODY_KEY_VAL_SEP, 1, false);
            if (entry_parts.size() < 2) {
                std::ostringstream msg;
                msg << "incomplete key-value specification in entry #" << entry_count << " (missing \"=\")";            
                throw ConfigurationSyntaxError(this->compose_error_message(start_pos, msg.str()));
            }
            this->entries_[textutil::strip(entry_parts[0], WHITESPACE)] = textutil::strip(entry_parts[1], WHITESPACE) ;
        }
    }
    this->is_block_set_ = true;
}

///////////////////////////////////////////////////////////////////////////////
// ConfigurationBlock inserter

std::istream& operator>> (std::istream& in, ConfigurationBlock& cblock) {
    cblock.parse(in);
    return in;
}

///////////////////////////////////////////////////////////////////////////////
// Configurator

Configurator::Configurator(const ConfigurationBlock& cb, 
                           unsigned long block_start_pos, 
                           unsigned long block_end_pos) 
        : configuration_block_(cb),
          block_start_pos_(block_start_pos),
          block_end_pos_(block_end_pos) { }
          
Configurator::~Configurator() { }

ConfigurationError Configurator::build_exception(const std::string& message) const {
    std::ostringstream msg;
    msg << "block \"" << this->configuration_block_.get_name();
    msg << "\" (file position " << this->block_start_pos_ + 1;
    msg << " to position " << this->block_end_pos_ + 1 << ")";
    msg << ": " << message;
    return ConfigurationError(msg.str());
}

///////////////////////////////////////////////////////////////////////////////
// WorldConfigurator

WorldConfigurator::WorldConfigurator(const ConfigurationBlock& cb, 
                     unsigned long block_start_pos, 
                     unsigned long block_end_pos) 
        : Configurator(cb, block_start_pos, block_end_pos),
          size_x_(0),
          size_y_(0),
          num_fitness_factors_(0),
          rand_seed_(0) {
    this->parse();
}

void WorldConfigurator::parse()  {
    this->size_x_ = this->get_configuration_value<unsigned long>("nrows"); 
    this->size_y_ = this->get_configuration_value<unsigned long>("ncols");
    this->num_fitness_factors_ = this->get_configuration_value<unsigned>("nfitness", MAX_FITNESS_FACTORS);
    this->rand_seed_ = this->get_configuration_value<unsigned>("seed", 0);
}

void WorldConfigurator::configure(World& world)  {
    world.set_random_seed(this->rand_seed_);
    world.set_num_fitness_factors(this->num_fitness_factors_);
    world.generate_landscape(this->size_x_, this->size_y_);    
}

///////////////////////////////////////////////////////////////////////////////
// SpeciesConfigurator

SpeciesConfigurator::SpeciesConfigurator(const ConfigurationBlock& cb, 
                     unsigned long block_start_pos, 
                     unsigned long block_end_pos) 
        : Configurator(cb, block_start_pos, block_end_pos),
          mutation_rate_(0),
          max_mutation_size_(0),
          mean_reproductive_rate_(0),
          reproductive_rate_mutation_size_(0),
          movement_capacity_(0) {
    this->parse();
}

void SpeciesConfigurator::parse()  {

}

void SpeciesConfigurator::configure(Species& world)  {
   
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
        : fsrc_(fpath),
          src_(fsrc_) {
    if (not this->src_) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw ConfigurationIOError(msg.str());
    }          
}

ConfigurationFile::ConfigurationFile(const std::string& fpath) 
        : fsrc_(fpath.c_str()),
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
    ConfigurationBlock cb;
    unsigned long block_start_pos = 0;
    unsigned long block_end_pos = 0;
    unsigned num_worlds = 0;
    while (not this->src_.eof()) {
        block_start_pos = this->src_.tellg();
        this->src_ >> cb;
        block_end_pos = this->src_.tellg();
        if (cb.is_block_set()) {
            if (cb.get_type() == "world") {
                if (num_worlds != 0) {
                    throw ConfigurationIOError("multiple definitions of world found");
                }                
                WorldConfigurator wcf(cb, block_start_pos, block_end_pos);
                wcf.configure(world);                
            }
        }                    
    }
}

} // namespace gingko
