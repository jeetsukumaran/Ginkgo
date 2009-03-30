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

#if !defined(GINGKO_ASCIIGRID_H)
#define GINGKO_ASCIIGRID_H

#include <vector>
#include <fstream>
#include <istream>
#include <stdexcept>

namespace gingko {
namespace asciigrid {

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
 * of unsigned longs.
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
 *     ÅThe nodata_value is the value in the ASCII file to be assigned to
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
 * point values. xllcorner and yllcorner are given as the EDGES of the grid, 
 * NOT the centers of the edge cells. ARC/INFO supports other header strings 
 * that allow the centers of the edge cells to be given using xllcenter and 
 * yllcenter instead. The origin of the grid is the upper left and terminus 
 * at the lower right.
 *
 */
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
        unsigned long get_ncols() {
            if (not this->is_metadata_loaded_) {
                this->parse_metadata_();
            }
            return this->ncols_;
        }
            
        /**
         * Returns number of rows in grid.
         * @return number of rows in grid
         */
        unsigned long get_nrows() {
            if (not this->is_metadata_loaded_) {
                this->parse_metadata_();
            }        
            return this->nrows_;
        }
        
        /**
         * Returns values of cells loaded from grid.
         * @return vector of values of cells loaded from grid
         */
        std::vector<long> get_cell_values() {
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
        bool has_size(unsigned long x, unsigned long y) {
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
         * @param metadata_value    store value of metadata (as unsigned long)
         */       
        void read_metadata_(std::string& metadata_name, unsigned long& metadata_value);        
        
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
        unsigned long       ncols_;
        /** Number of rows in the grid. */
        unsigned long       nrows_;
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
        std::vector<long>   cell_values_;

};

} // asciigrid
} // gingko

#endif
