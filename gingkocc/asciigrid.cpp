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

#include "textutils.h"
#include "asciigrid.h"
#include <sstream>
#include <cstdlib>

namespace gingko {
namespace asciigrid {

AsciiGrid::AsciiGrid() {
    this->ncols_ = 0;
    this->nrows_ = 0;
    this->xllcorner_ = 0;
    this->yllcorner_ = 0;
    this->xllcenter_ = 0;
    this->yllcenter_ = 0;    
    this->cell_size_ = 0;
    this->nodata_value_ = -99999;
}

AsciiGrid::~AsciiGrid() { }

// std::string AsciiGrid::get_token(std::istream& src, const char * eof_msg) {
//     std::string token;
//     src >> token;
//     if (src.eof() && src.fail()) {        
//         throw AsciiGridFormatEofError(eof_msg);
//     }
//     return token;
// }
// 
// std::string AsciiGrid::read_metadata_name(std::istream& src, const char * expected_token) {
//     std::string token = this->get_token(src, "unexpected EOF while try to read metadata row");
//     if (expected_token != NULL) {
//         if (token != expected_token) {
//             std::ostringstream msg;
//             msg << "expecting \"" << expected_token << "\" but found \"" << token << "\"";
//             throw AsciiGridFormatTokenError(msg.str());
//         }
//     }
//     return token;
// }
// 
// void AsciiGrid::read_metadata_value(std::istream& src, const char * expected_token, unsigned long& val) {
//     this->read_metadata_name(expected_token);
//     std::string val_str = 
// }
// 
// void AsciiGrid::read_metadata_value(std::istream& src, const char * expected_token, long& val) {
// 
// }
// 

void AsciiGrid::read_metadata(std::istream& src, std::string& metadata_name, long& metadata_value) {
    src >> metadata_name;
    if (src.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    src >> metadata_value;
    if (src.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    if (src.fail()) {    
        throw AsciiGridFormatValueError("Value error while reading metadata");
    }
}

void AsciiGrid::read_metadata(std::istream& src, std::string& metadata_name, unsigned long& metadata_value) {
    src >> metadata_name;
    if (src.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    src >> metadata_value;
    if (src.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    if (src.fail()) {    
        throw AsciiGridFormatValueError("Value error while reading metadata");
    }
}

void AsciiGrid::read_metadata(std::istream& src, std::string& metadata_name, float& metadata_value) {
    src >> metadata_name;
    if (src.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    src >> metadata_value;
    if (src.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    if (src.fail()) {    
        throw AsciiGridFormatValueError("value error while reading metadata");
    }
}

void AsciiGrid::parse_metadata(std::istream& src) {
    std::string name;
    
    this->read_metadata(src, name, this->ncols_);
    if (textutils::lower(name) != "ncols") {
        std::ostringstream msg;
        msg << "expecting \"ncols\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }    
    
    this->read_metadata(src, name, this->nrows_);
    if (textutils::lower(name) != "nrows") {
        std::ostringstream msg;
        msg << "expecting \"nrows\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }
    
    this->read_metadata(src, name, this->xllcorner_);
    if (textutils::lower(name) == "xllcorner") {
        // pass
    } else if (textutils::lower(name) == "xllcenter") {
        this->xllcenter_ = this->xllcorner_;
        this->xllcorner_ = 0;
    } else {
        std::ostringstream msg;
        msg << "expecting \"xllcorner\" or \"xllcenter\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }
    
    this->read_metadata(src, name, this->yllcorner_);
    if (textutils::lower(name) == "yllcorner") {
        // pass
    } else if (textutils::lower(name) == "yllcenter") {
        this->yllcenter_ = this->yllcorner_;
        this->yllcorner_ = 0;
    } else {
        std::ostringstream msg;
        msg << "expecting \"yllcorner\" or \"yllcenter\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }
    
    this->read_metadata(src, name, this->cell_size_);
    if (textutils::lower(name) != "cellsize") {
        std::ostringstream msg;
        msg << "expecting \"cellsize\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }       
}

void AsciiGrid::parse(std::istream& src) {     

    // process mandatory metadata
    this->parse_metadata(src);
    
    // reserve memory
    unsigned long x = 0;
    unsigned long y = 0;
    this->values_.clear();
    this->values_.reserve(this->nrows_ * this->ncols_);
    
    // process option nodata value
    std::string token;
    src >> token;
    if (src.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading data");
    }    
    if (textutils::lower(token) == "nodata_value") {
        src >> this->nodata_value_;
        if (src.eof()) {
            throw AsciiGridFormatEofError("unexpected EOF while reading data");
        }
        if (src.fail()) {
            throw AsciiGridFormatValueError("invalid value specified for NODATA_value");
        }           
    } else {
        std::istringstream s(token);
        long v = 0;
        s >> v;
        if (s.fail()) {
            std::ostringstream msg;
            msg << "invalid value specified for cell (0,0) in file character position " << src.tellg();
            throw AsciiGridFormatValueError(msg.str());
        }
        this->values_.push_back(v);
    }
    
    unsigned long cell_value;
    src >> cell_value;
    while (not src.eof()) {
        if (src.fail()) {
            std::ostringstream msg;
            msg << "invalid value specified for cell (" << x << ", " << y << ") in file character position " << src.tellg();
            throw AsciiGridFormatValueError(msg.str());
        }
        this->values_.push_back(cell_value);
        x += 1;
        if (x >= this->ncols_) {
            x = 0;
            y += 1;
        }
        src >> cell_value;
    }
    
}

} // asciigrid
} // gingko
