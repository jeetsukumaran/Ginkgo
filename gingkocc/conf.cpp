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

#include "conf.h"
#include "textutils.h"

namespace gingko {

const char * BLOCK_START = "@";
const char * BLOCK_BODY_START = "{";
const char * BLOCK_BODY_END = "}";
const char * BLOCK_BODY_LINE_TERM = ";\n";
const char * BLOCK_BODY_KEY_VAL_SEP = "=";
const char * WHITESPACE = " \t\n";

///////////////////////////////////////////////////////////////////////////////
// ConfigurationBlockParser

ConfigurationFileParser::ConfigurationFileParser(std::istream& src)
        : src_(src) {
    if (not this->src_) {
        throw ConfigurationIOError("invalid source stream");
    }
}

ConfigurationFileParser::ConfigurationFileParser(const char * fpath)
        : fsrc_(fpath),
          src_(fsrc_) {
    if (not this->src_) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw ConfigurationIOError(msg.str());
    }          
}

ConfigurationFileParser::ConfigurationFileParser(const std::string& fpath) 
        : fsrc_(fpath.c_str()),
          src_(fsrc_) {
    if (not this->src_) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw ConfigurationIOError(msg.str());
    }
}

ConfigurationFileParser::~ConfigurationFileParser() { }

///////////////////////////////////////////////////////////////////////////////
// ConfigurationBlockParser

// default constructor
ConfigurationBlockParser::ConfigurationBlockParser() { }

// construct and parse
ConfigurationBlockParser::ConfigurationBlockParser(std::istream& in) {   
    this->parse(in);
}        

// default do-nothing destructor
ConfigurationBlockParser::~ConfigurationBlockParser() {}   

// return type
std::string ConfigurationBlockParser::get_type() const {
    return this->type_;
}

// return name
std::string ConfigurationBlockParser::get_name() const {
    return this->name_;
}

// get entry values by keys
std::string ConfigurationBlockParser::get_entry(const std::string& key) const {
    std::map< std::string, std::string >::const_iterator val = this->entries_.find(key);
    assert(val != this->entries_.end());
    return val->second;
}

// get keys
std::vector<std::string> ConfigurationBlockParser::get_keys() const {
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
std::string ConfigurationBlockParser::compose_error_message(unsigned long pos, const char * desc) {
    std::ostringstream msg;
    msg << "Block starting at character position " <<  pos+1 << ": ";
    msg << desc;
    return msg.str();
 }
 
 // using std string object
std::string ConfigurationBlockParser::compose_error_message(unsigned long pos, const std::string& desc) {
    return this->compose_error_message(pos, desc.c_str());
}

// workhorse parser
void ConfigurationBlockParser::parse(std::istream& in) {                
    std::string raw;
    unsigned long start_pos(in.tellg());
    std::getline(in, raw, BLOCK_BODY_END[0]);
    raw = textutils::strip(raw);
    
    if (raw.size() == 0) {
        throw ConfigurationSyntaxError(this->compose_error_message(start_pos, "empty block"));                
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
    
    std::vector<std::string> parts = textutils::split(raw, BLOCK_BODY_START, false);
    
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
    
    std::vector<std::string> head_parts = textutils::split_on_any(textutils::strip(parts[0]), WHITESPACE, false);
    
    if (head_parts.size() < 2) {
        throw ConfigurationSyntaxError(this->compose_error_message(start_pos, "found only one element in block header, but expecting two (type and name)"));
    }
    
    if (head_parts.size() > 2) {
        throw ConfigurationSyntaxError(this->compose_error_message(start_pos, "found multiple elements in block header, but expecting only two (type and name)"));
    }
    
    this->type_ = textutils::strip(head_parts[0]);
    this->name_ = textutils::strip(head_parts[1]);
    
    unsigned entry_count = 0;
    std::vector<std::string> body_parts = textutils::split_on_any(textutils::strip(parts[1]), BLOCK_BODY_LINE_TERM, false);
    for (std::vector<std::string>::const_iterator s = body_parts.begin(); s != body_parts.end(); ++s) {
        std::string entry = textutils::strip(*s, WHITESPACE);
        entry_count += 1;
        if (entry.size() > 0) {
            std::vector<std::string> entry_parts = textutils::split(entry, BLOCK_BODY_KEY_VAL_SEP, 1, false);
            if (entry_parts.size() < 2) {
                std::ostringstream msg;
                msg << "incomplete key-value specification in entry #" << entry_count << " (missing \"=\")";            
                throw ConfigurationSyntaxError(this->compose_error_message(start_pos, msg.str()));
            }
            this->entries_[textutils::strip(entry_parts[0], WHITESPACE)] = textutils::strip(entry_parts[1], WHITESPACE) ;
        }
    }
}

} // namespace gingko
