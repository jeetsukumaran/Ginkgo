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

#if !defined(GINGKO_DEFS_H)
#define GINGKO_DEFS_H

namespace gingko {

//! The maximum number of fitness factors in the system. A "fitness factor"
//! represents an environmental variable and a corresponding genotypic variable
//! that together contribute to the fitness of an individual organism.
const unsigned MAX_FITNESS_FACTORS = 10;

//! The number of neutral diploid locii to track
const unsigned NUM_NEUTRAL_DIPLOID_LOCII = 10;

//! The value type of a fitness factor (both environmental and genotypic.
typedef int                     FitnessFactorType;

//! The number of dimensions or "slots" reserved for fitness assessment. The
//! actual number may be less than this.
typedef FitnessFactorType       FitnessFactors[MAX_FITNESS_FACTORS];

//! The units for referencing cells on the landscape. We allow signed values
//! even though a negative coordinate is invalide so that negative values can 
//! be stored and used in calculations.
typedef long                    CellIndexType;

} // gingko namespace

#endif
