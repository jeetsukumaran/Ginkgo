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

#include "gingko_defs.hpp"
#include "confsys.hpp"
#include "textutil.hpp"
#include "convert.hpp"
#include "world.hpp"
#include "biosys.hpp"

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
const char * BLOCK_BODY_LINE_TERM = ";";
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
                msg << "incomplete key-value specification in entry no. " << entry_count << " (missing \"=\"): ";
                msg << "\"" << entry << "\"";
                throw ConfigurationSyntaxError(this->compose_error_message(start_pos, msg.str()));
            }
            this->entries_[textutil::strip(entry_parts[0], WHITESPACE)] = textutil::strip(entry_parts[1], WHITESPACE) ;
        }
    }
    this->is_block_set_ = true;
    this->block_start_pos_ = start_pos;
    this->block_end_pos_ = in.tellg();
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

void Configurator::get_configuration_positions(const std::string& key, OrganismDistribution& od) {
    std::vector<std::string> positions = this->get_configuration_vector<std::string>(key);
    od.x.reserve(od.x.size() + positions.size());
    od.y.reserve(od.y.size() + positions.size());
    for (std::vector<std::string>::const_iterator pos = positions.begin(); pos < positions.end(); ++pos) {
        std::string pos_clean = textutil::strip(*pos, WHITESPACE);
        if (pos_clean.size() == 0) {
            continue;
        }
        std::vector<std::string> pos_parts = textutil::split(pos_clean, ",", 1, false);
        if (pos_parts.size() != 2) {
            throw this->build_exception("invalid or incomplete position specification: \"" + pos_clean + "\"");
        }
        try {
            od.x.push_back(convert::to_scalar<unsigned long>(pos_parts[0]));
        } catch (const convert::ValueError& e) {
            throw this->build_exception("invalid x-coordinate value for position: \"" + pos_parts[0] + "\"");
        }
        try {
            od.y.push_back(convert::to_scalar<unsigned long>(pos_parts[1]));
        } catch (const convert::ValueError& e) {
            throw this->build_exception("invalid y-coordinate value for position: \"" + pos_parts[1] + "\"");
        }         
    }
}

ConfigurationError Configurator::build_exception(const std::string& message) const {
    return confsys_detail::build_configuration_block_exception(this->configuration_block_, message);
}

///////////////////////////////////////////////////////////////////////////////
// WorldConfigurator

WorldConfigurator::WorldConfigurator(const ConfigurationBlock& cb) 
        : Configurator(cb),
          size_x_(0),
          size_y_(0),
          num_fitness_factors_(0),
          rand_seed_(0) {
    this->parse();
}

void WorldConfigurator::parse()  {
    this->size_x_ = this->get_configuration_scalar<unsigned long>("nrows"); 
    this->size_y_ = this->get_configuration_scalar<unsigned long>("ncols");
    this->generations_to_run_ = this->get_configuration_scalar<unsigned long>("ngens");    
    this->num_fitness_factors_ = this->get_configuration_scalar<unsigned>("nfitness", MAX_FITNESS_FACTORS);
    this->rand_seed_ = this->get_configuration_scalar<unsigned>("rseed", 0);
}

void WorldConfigurator::configure(World& world)  {
    world.set_label(this->get_name());
    world.set_random_seed(this->rand_seed_);
    world.set_generations_to_run(this->generations_to_run_);
    world.set_num_fitness_factors(this->num_fitness_factors_);
    world.generate_landscape(this->size_x_, this->size_y_);    
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
    if (sp.get_num_fitness_factors() > this->default_genotypic_fitness_factors_.size()) {
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
            world.seed_population(cell_index, label, od.num_organisms);
        }
    }
}

void SpeciesConfigurator::process_seed_populations() {
    std::vector<std::string> keys = this->get_matching_configuration_keys("init");
    for (std::vector<std::string>::const_iterator k = keys.begin(); k != keys.end(); ++k) {
        this->seed_populations_.push_back(OrganismDistribution());
        OrganismDistribution& od = this->seed_populations_.back();
        std::string key = *k;
        std::vector<std::string> key_parts = textutil::split(key, "#", 1, false);
        if (key_parts.size() < 2) {
            throw this->build_exception("need to specify number of organisms in starting population using '#' token: \"" + key + "\"");
        }
        try {
            od.num_organisms = convert::to_scalar<unsigned long>(key_parts[1]);
        } catch (const convert::ValueError& e) {
            throw this->build_exception("invalid value for population size: \"" + key_parts[1] + "\"");
        }
        this->get_configuration_positions(key, od);
    }
}

///////////////////////////////////////////////////////////////////////////////
// GenerationConfigurator

GenerationConfigurator::GenerationConfigurator(const ConfigurationBlock& cb) 
        : Configurator(cb)  {
    this->parse();
}

void GenerationConfigurator::parse()  {

}

void GenerationConfigurator::configure(World& world)  {

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
    unsigned num_worlds = 0;
    std::set<std::string> species_labels;
    while (not this->src_.eof()) {
        this->src_ >> cb;
        if (cb.is_block_set()) {
            if (cb.get_type() == "world") {
                if (num_worlds != 0) {
                    throw confsys_detail::build_configuration_block_exception(cb, "world must be described before species are added");
                }                
                WorldConfigurator wcf(cb);
                wcf.configure(world);
                num_worlds += 1;
            }
            if (cb.get_type() == "species") {
                if (num_worlds == 0) {
                    throw confsys_detail::build_configuration_block_exception(cb, "world must be described before species are added");
                }
                std::string label = cb.get_name();
                if (species_labels.find(label) != species_labels.end()) {
                    throw confsys_detail::build_configuration_block_exception(cb, "species \"" + label + "\" has already been defined");
                }
                SpeciesConfigurator spcf(cb);
                spcf.configure(world);
                species_labels.insert(label);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Implementation details

namespace confsys_detail {
    
    /**
     * Composes and returns and appropriate exception.
     * @param message           error message
     * @param cb                ConfigurationBlock that has the error
     * @return                  ConfiguratonError exception to be thrown
     */
    ConfigurationError build_configuration_block_exception(const ConfigurationBlock& cb,
            const std::string& message) {
        std::ostringstream msg;
        msg << cb.get_type() << " block \"" << cb.get_name();
        msg << "\" (file position " << cb.get_block_start_pos() + 1;
        msg << " to position " << cb.get_block_end_pos() + 1 << ")";
        msg << ": " << message;
        return ConfigurationError(msg.str());    
    }
    
} // confsys_detail


} // namespace gingko
