///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Biogeographical Evolution Simulator.
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

#if !defined(GINKGO_ASCIIGRID_H)
#define GINKGO_ASCIIGRID_H

#include <vector>
#include <fstream>
#include <istream>
#include <sstream>
#include <cassert>
#include <stdexcept>
#include <cstdlib>
#include <iomanip>
#include "textutil.hpp"
#include "asciigrid.hpp"

namespace ginkgo {
namespace asciigrid {

typedef unsigned long index_type;

/**
 * Writes out a vector of longs as an ASCII Grid file.
 * @param   vals    cell values
 * @param   ncols   number of cols
 * @param   nrows   number of rows
 * @param   out     destination
 */
template <typename T>
void write_grid(const std::vector<T>& vals,
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

/**
 * General i/o error.
 */
class AsciiGridIOError : public std::runtime_error {
    public:
        AsciiGridIOError(const char * msg) : std::runtime_error(msg) {}
        AsciiGridIOError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

/**
 * General format error.
 */
class AsciiGridFormatError : public std::runtime_error {
    public:
        AsciiGridFormatError(const char * msg) : std::runtime_error(msg) {}
        AsciiGridFormatError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

/**
 * EOF error.
 */
class AsciiGridFormatEofError : public AsciiGridFormatError {
    public:
        AsciiGridFormatEofError(const char * msg) : AsciiGridFormatError(msg) {}
        AsciiGridFormatEofError(const std::string& msg) : AsciiGridFormatError(msg) {}
};

/**
 * Unexpected token error.
 */
class AsciiGridFormatTokenError : public AsciiGridFormatError {
    public:
        AsciiGridFormatTokenError(const char * msg) : AsciiGridFormatError(msg) {}
        AsciiGridFormatTokenError(const std::string& msg) : AsciiGridFormatError(msg) {}
};

/**
 * Value error.
 */
class AsciiGridFormatValueError : public AsciiGridFormatError {
    public:
        AsciiGridFormatValueError(const char * msg) : AsciiGridFormatError(msg) {}
        AsciiGridFormatValueError(const std::string& msg) : AsciiGridFormatError(msg) {}
};


/**
 * Encapsulates parsing of an ESRI ASCII Grid format file into a vector
 * of index_types.
 *
 * (From the ArcWorkstation 8.3 Help File)-------------------------------------
 *
 *     The ASCII file must consist of header information containing a set of
 *     keywords, followed by cell values in row-major order. The file format
 *     is:
 *
 *         <NCOLS xxx>
 *         <NROWS xxx>
 *         <XLLCENTER xxx | XLLCORNER xxx>
 *         <YLLCENTER xxx | YLLCORNER xxx>
 *         <CELLSIZE xxx>
 *         {NODATA_VALUE xxx}
 *         row 1
 *         row 2
 *         .
 *         .
 *         .
 *         row n
 *
 *     where xxx is a number, and the keyword nodata_value is optional and
 *     defaults to -9999. Row 1 of the data is at the top of the grid, row 2
 *     is just under row 1 and so on.
 *
 *     For example:
 *
 *         ncols 480
 *         nrows 450
 *         xllcorner 378923
 *         yllcorner 4072345
 *         cellsize 30
 *         nodata_value -32768
 *         43 2 45 7 3 56 2 5 23 65 34 6 32 54 57 34 2 2 54 6
 *         35 45 65 34 2 6 78 4 2 6 89 3 2 7 45 23 5 8 4 1 62 ...
 *
 *     The nodata_value is the value in the ASCII file to be assigned to
 *     those cells whose true value is unknown. In the grid they will be
 *     assigned the keyword NODATA.
 *
 *     Cell values should be delimited by spaces. No carriage returns are
 *     necessary at the end of each row in the grid. The number of columns in
 *     the header is used to determine when a new row begins.
 *
 *     The number of cell values must be equal to the number of rows times
 *     the number of columns, or an error will be returned.
 * ----------------------------------------------------------------------------
 *
 * Note that: Grid values are stored as integers but can be read as floating
 * point values. *** @JS: Extended to a store grid values as floats ***
 * xllcorner and yllcorner are given as the EDGES of the grid,
 * NOT the centers of the edge cells. ARC/INFO supports other header strings
 * that allow the centers of the edge cells to be given using xllcenter and
 * yllcenter instead. The origin of the grid is the upper left and terminus
 * at the lower right.
 *
 */

template <class T=long>
class AsciiGrid {

    public:

        /**
         * Initializes metadata and binds to source stream.
         *
         * @param src   data source
         */
        AsciiGrid(std::istream& src);

        /**
         * Initializes metadata and binds to source file.
         *
         * @param fpath filepath of data source
         */
        AsciiGrid(const char * fpath);

        /**
         * Initializes metadata and binds to source file.
         *
         * @param fpath filepath of data source
         */
        AsciiGrid(const std::string& fpath);

        /**
         * Default no-op destructor.
         */
        ~AsciiGrid();

        /**
         * Returns number of columns in grid.
         * @return number of columns in grid
         */
        index_type get_ncols() {
            if (not this->is_metadata_loaded_) {
                this->parse_metadata_();
            }
            return this->ncols_;
        }

        /**
         * Returns number of rows in grid.
         * @return number of rows in grid
         */
        index_type get_nrows() {
            if (not this->is_metadata_loaded_) {
                this->parse_metadata_();
            }
            return this->nrows_;
        }

        /**
         * Returns values of cells loaded from grid.
         * @return vector of values of cells loaded from grid
         */
        std::vector<T> get_cell_values() {
            if (not this->is_cell_values_loaded_) {
                this->parse_cell_values_();
            }
            return this->cell_values_;
        }

        /**
         * Returns <code>true</code> if grid metadata specifies the same
         * dimensions as that required.
         * @param   x   required x dimension of grid
         * @param   y   required y dimension of grid
         * @return      <code>true</code> if grid dimensions are as specified
         */
        bool has_size(index_type x, index_type y) {
            if (this->get_ncols() == x && this->get_nrows() == y) {
                return true;
            } else {
                return false;
            }
        }

    private:

        /**
         * Initializes metadata.
         */
        void init_();

        /**
         * Read and parse a metadata row into into its components.
         *
         * @param metadata_name     store name of metadata
         * @param metadata_value    store value of metadata (as long)
         */
        void read_metadata_(std::string& metadata_name, long& metadata_value);

        /**
         * Read and parse a metadata row into into its components.
         *
         * @param metadata_name     store name of metadata
         * @param metadata_value    store value of metadata (as index_type)
         */
        void read_metadata_(std::string& metadata_name, index_type& metadata_value);

        /**
         * Read and parse a metadata row into into its components.
         *
         * @param metadata_name     store name of metadata
         * @param metadata_value    store value of metadata (as float)
         */
        void read_metadata_(std::string& metadata_name, float& metadata_value);

        /**
         * Load metadata from an input stream.
         */
         void parse_metadata_();

        /**
         * Load cell values from an input stream into given structure.
         *
         * @param cell_values   container to which extracted values will be
         *                      stored (existing values will be cleared)
         */
         void parse_cell_values_();

    private:
        /** Disabled copy constructor. */
        AsciiGrid(const AsciiGrid&);
        /** Disabled assignment operator. */
        const AsciiGrid& operator=(const AsciiGrid&);

    private:
        /** Input (file) stream. */
        std::ifstream       fsrc_;
        /** Input stream. */
        std::istream&       src_;
        /** Tracks whether or not metadata has been parsed. */
        bool                is_metadata_loaded_;
        /** Tracks whether or not metadata has been parsed. */
        bool                is_cell_values_loaded_;
        /** Number of columns in the grid. */
        index_type       ncols_;
        /** Number of rows in the grid. */
        index_type       nrows_;
        /** Geographical X-coordinate of the lower-left cell corner of the grid. */
        float               xllcorner_;
        /** Geographical Y-coordinate of the lower-left cell corner of the grid.  */
        float               yllcorner_;
        /** Geographical X-coordinate of the lower-left cell center of the grid. */
        float               xllcenter_;
        /** Geographical Y-coordinate of the lower-left cell center of the grid.  */
        float               yllcenter_;
        /** Size of cell.  */
        float               cell_size_;
        /** Value of flag indicating missing data. */
        long                nodata_value_;
        /** Cell values. */
        std::vector<T>   cell_values_;

};


template <class T>
AsciiGrid<T>::AsciiGrid(std::istream& src)
        : src_(src) {
    if (not this->src_) {
        throw AsciiGridIOError("invalid source stream");
    }
    this->init_();
}

template <class T>
AsciiGrid<T>::AsciiGrid(const char * fpath)
        : fsrc_(fpath),
          src_(fsrc_) {
    if (not this->src_) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw AsciiGridIOError(msg.str());
    }
    this->init_();
}

template <class T>
AsciiGrid<T>::AsciiGrid(const std::string& fpath)
        : fsrc_(fpath.c_str()),
          src_(fsrc_) {
    if (not this->src_) {
        std::ostringstream msg;
        msg << "invalid source \"" << fpath << "\"";
        throw AsciiGridIOError(msg.str());
    }
    this->init_();
}

template <class T>
AsciiGrid<T>::~AsciiGrid() { }

template <class T>
void AsciiGrid<T>::init_() {
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

template <class T>
void AsciiGrid<T>::read_metadata_(std::string& metadata_name, long& metadata_value) {
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

template <class T>
void AsciiGrid<T>::read_metadata_(std::string& metadata_name, index_type& metadata_value) {
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

template <class T>
void AsciiGrid<T>::read_metadata_(std::string& metadata_name, float& metadata_value) {
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

template <class T>
void AsciiGrid<T>::parse_metadata_() {
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

template <class T>
void AsciiGrid<T>::parse_cell_values_() {

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
        T v = 0;
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

    T cell_value;
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
} // ginkgo

#endif
