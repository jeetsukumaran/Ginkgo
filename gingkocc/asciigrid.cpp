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
#include <sstream>
#include <cstdlib>
#include <cassert>
#include "textutils.hpp"
#include "asciigrid.hpp"

namespace gingko {
namespace asciigrid {

AsciiGrid::AsciiGrid(std::istream& src)
        : src_(src) {
    if (not this->src_) {
        throw AsciiGridIOError("invalid source stream");
    }
    this->init_();
}

AsciiGrid::AsciiGrid(const char * fpath)
        : fsrc_(fpath),
          src_(fsrc_) {
    if (not this->src_) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw AsciiGridIOError(msg.str());
    }          
    this->init_();
}

AsciiGrid::AsciiGrid(const std::string& fpath) 
        : fsrc_(fpath.c_str()),
          src_(fsrc_) {
    if (not this->src_) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw AsciiGridIOError(msg.str());
    }              
    this->init_();
}

AsciiGrid::~AsciiGrid() { }

void AsciiGrid::init_() {
    this->is_metadata_loaded_ = false;
    this->is_cell_values_loaded_ = false;
    this->ncols_ = 0;
    this->nrows_ = 0;
    this->xllcorner_ = 0;
    this->yllcorner_ = 0;
    this->xllcenter_ = 0;
    this->yllcenter_ = 0;    
    this->cell_size_ = 0;
    this->nodata_value_ = -99999;
}

void AsciiGrid::read_metadata_(std::string& metadata_name, long& metadata_value) {
    assert(this->src_);
    this->src_ >> metadata_name;
    if (this->src_.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    this->src_ >> metadata_value;
    if (this->src_.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    if (this->src_.fail()) {    
        throw AsciiGridFormatValueError("value error while reading metadata (expecting long)");
    }
}

void AsciiGrid::read_metadata_(std::string& metadata_name, unsigned long& metadata_value) {
    assert(this->src_);
    this->src_ >> metadata_name;
    if (this->src_.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    this->src_ >> metadata_value;
    if (this->src_.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    if (this->src_.fail()) {    
        throw AsciiGridFormatValueError("value error while reading metadata (expecting unsigned long)");
    }
}

void AsciiGrid::read_metadata_(std::string& metadata_name, float& metadata_value) {
    assert(this->src_);
    this->src_ >> metadata_name;
    if (this->src_.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    this->src_ >> metadata_value;
    if (this->src_.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading metadata");
    }
    if (this->src_.fail()) {    
        throw AsciiGridFormatValueError("value error while reading metadata (expecting float)");
    }
}

void AsciiGrid::parse_metadata_() {
    std::string name;
    
    this->read_metadata_(name, this->ncols_);
    if (textutils::lower(name) != "ncols") {
        std::ostringstream msg;
        msg << "expecting \"ncols\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }    
    
    this->read_metadata_(name, this->nrows_);
    if (textutils::lower(name) != "nrows") {
        std::ostringstream msg;
        msg << "expecting \"nrows\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }
    
    this->read_metadata_(name, this->xllcorner_);
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
    
    this->read_metadata_(name, this->yllcorner_);
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
    
    this->read_metadata_(name, this->cell_size_);
    if (textutils::lower(name) != "cellsize") {
        std::ostringstream msg;
        msg << "expecting \"cellsize\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }
    
    this->is_metadata_loaded_ = true;
}

void AsciiGrid::parse_cell_values_() {   

    if (not this->is_metadata_loaded_) {
        this->parse_metadata_();
    }
    assert(this->ncols_ != 0);    
    assert(this->nrows_ != 0);
    
    // reserve memory
    this->cell_values_.clear();
    this->cell_values_.reserve(this->nrows_ * this->ncols_);
    
    // track position (for error reporting)
    unsigned long x = 0;
    unsigned long y = 0;
    
    // process optional nodata value
    std::string token;
    this->src_ >> token;
    if (this->src_.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading data");
    }    
    if (textutils::lower(token) == "nodata_value") {
        this->src_ >> this->nodata_value_;
        if (this->src_.eof()) {
            throw AsciiGridFormatEofError("unexpected EOF while reading data");
        }
        if (this->src_.fail()) {
            throw AsciiGridFormatValueError("invalid value specified for NODATA_value");
        }           
    } else {
        std::istringstream s(token);
        long v = 0;
        s >> v;
        if (s.fail()) {
            std::ostringstream msg;
            msg << "invalid value specified for cell (0,0) in file character position " << this->src_.tellg();
            throw AsciiGridFormatValueError(msg.str());
        }
        this->cell_values_.push_back(v);        
        x+= 1;
        if (x > this->ncols_) {
            x = 0;
            y += 1;
        }
    }
    
    long cell_value;
    this->src_ >> cell_value;
    while (not this->src_.eof()) {
        if (this->src_.fail()) {
            std::ostringstream msg;
            msg << "invalid value specified for cell (" << x << ", " << y << ") in file character position " << this->src_.tellg();
            throw AsciiGridFormatValueError(msg.str());
        }
        this->cell_values_.push_back(cell_value);
        x += 1;
        if (x >= this->ncols_) {
            x = 0;
            y += 1;
        }
        this->src_ >> cell_value;
    }
    this->is_cell_values_loaded_ = true;
    
}

} // asciigrid
} // gingko
