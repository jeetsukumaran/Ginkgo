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

// default constructor
ConfigurationBlockParser::ConfigurationBlockParser() { }

// construct and parse
ConfigurationBlockParser::ConfigurationBlockParser(std::istream& in) {   
    this->parse(in);
}        

// default do-nothing destructor
ConfigurationBlockParser::~ConfigurationBlockParser() {}   

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

// workhorse
void ConfigurationBlockParser::parse(std::istream& in) {                
    std::string raw;
    unsigned long start_pos(in.tellg());
    std::getline(in, raw, END_BLOCK_BODY);
    raw = strip(raw);
    
    if (raw.size() == 0) {
        throw ConfigurationParseError(this->compose_error_message(start_pos, "empty block"));                
    }
    
    if (raw[0] != BLOCK_START) {
        std::ostringstream msg;
        msg << "does not begin with '" << BLOCK_START << "'";
        throw ConfigurationParseError(this->compose_error_message(start_pos, msg.str()));
    }
    
    if (in.eof()) {
        std::ostringstream msg;
        msg << "EOF before block body terminator ('" << END_BLOCK_BODY << "')";
        throw ConfigurationParseError(this->compose_error_message(start_pos, msg.str()));                        
    }
    
    std::vector<std::string> parts = split(raw, "{", false);
    
    if (parts.size() < 2) {
        std::ostringstream msg;
        msg << "missing block body initiator ('" << END_BLOCK_BODY << "')";    
        throw ConfigurationParseError(this->compose_error_message(start_pos, msg.str()));
    }
    
    if (parts.size() > 2) {
        std::ostringstream msg;
        msg << "multiple block body initiators ('" << END_BLOCK_BODY << "')";    
        throw ConfigurationParseError(this->compose_error_message(start_pos, msg.str()));
    }     
    
    std::vector<std::string> head_parts = split(parts[0], " ", true);
}

} // namespace gingko
