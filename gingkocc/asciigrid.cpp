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

#include "asciigrid.h"

namespace gingko {
namespace asciigrid {

AsciiGrid::AsciiGrid() {
    this->ncols_ = 0;
    this->nrows_ = 0;
    this->xllcorner_ = 0;
    this->yllcorner_ = 0;
    this->nodata_value_ = -99999;
}

AsciiGrid::~AsciiGrid() { }

void AsciiGrid::parse(std::istream& src) {
//     std::string token;
//     token = this->parse_value(src, "ncols", this->ncols);
//     token = this->parse_value(src, "nrows", this->nrows);
//     token = this->parse_value(src, "xllcorner", this->xllcorner);
//     token = this->parse_value(src, "yllcorner", this->yllcorner);
//     token = this->parse_value(src, "cellsize", this->cellsize);
//     token = this->parse_value(src, "NODATA_value", this->nodata_value);
}

} // asciigrid
} // gingko
