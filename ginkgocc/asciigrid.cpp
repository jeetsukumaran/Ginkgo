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
#include <iomanip>
#include "textutil.hpp"
#include "asciigrid.hpp"

namespace gingko {
namespace asciigrid {

void write_grid(const std::vector<long>& vals, 
        index_type ncols,
        index_type nrows, 
        std::ostream& out) {
    assert(vals.size() == ncols * nrows);        
    out << "ncols           " << ncols << std::endl;
    out << "nrows           " << nrows << std::endl;
    out << "xllcorner       0.0" << std::endl;
    out << "yllcorner       0.0" << std::endl;    
    out << "cellsize        1.0" << std::endl;
    out << "NODATA_value    -9999" << std::endl;
    index_type index = 0;
    for (index_type row = 0; row < nrows; ++row) {
        for (index_type col = 0; col < ncols; ++col, ++index) {
            out << vals[index] << " ";
        }
        out << std::endl;
    }
}        

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

void AsciiGrid::read_metadata_(std::string& metadata_name, index_type& metadata_value) {
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
        throw AsciiGridFormatValueError("value error while reading metadata (expecting index_type)");
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
    if (textutil::lower(name) != "ncols") {
        std::ostringstream msg;
        msg << "expecting \"ncols\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }    
    
    this->read_metadata_(name, this->nrows_);
    if (textutil::lower(name) != "nrows") {
        std::ostringstream msg;
        msg << "expecting \"nrows\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }
    
    this->read_metadata_(name, this->xllcorner_);
    if (textutil::lower(name) == "xllcorner") {
        // pass
    } else if (textutil::lower(name) == "xllcenter") {
        this->xllcenter_ = this->xllcorner_;
        this->xllcorner_ = 0;
    } else {
        std::ostringstream msg;
        msg << "expecting \"xllcorner\" or \"xllcenter\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }
    
    this->read_metadata_(name, this->yllcorner_);
    if (textutil::lower(name) == "yllcorner") {
        // pass
    } else if (textutil::lower(name) == "yllcenter") {
        this->yllcenter_ = this->yllcorner_;
        this->yllcorner_ = 0;
    } else {
        std::ostringstream msg;
        msg << "expecting \"yllcorner\" or \"yllcenter\" but found \"" << name << "\"";
        throw AsciiGridFormatTokenError(msg.str());
    }
    
    this->read_metadata_(name, this->cell_size_);
    if (textutil::lower(name) != "cellsize") {
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
    index_type x = 0;
    index_type y = 0;
    
    // process optional nodata value
    std::string token;
    this->src_ >> token;
    if (this->src_.eof()) {
        throw AsciiGridFormatEofError("unexpected EOF while reading data");
    }    
    if (textutil::lower(token) == "nodata_value") {
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
        if (x >= this->ncols_) {
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
    
//     out << std::setfill(' '); // reset    
//     std::cout << std::endl << "*** GRID DUMP ***" << std::endl;
//     index_type i = 0;
//     x = 0;
//     y = 0;
//     for (y = 0; y < this->nrows_; ++y) {
//         for (x = 0; x < this->ncols_; ++x, ++i) {
//             std::cout << std::setw(4) << this->cell_values_.at(i);
//         }
//         std::cout << std::endl;
//     }
}

} // asciigrid
} // gingko
