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

#include <fstream>
#include <stdexcept>
#include <climits>

#if !defined(GINKGO_DEFS_H)
#define GINKGO_DEFS_H

#define assert_m(cond, out, msg) if (!(cond)) { out << msg << std::endl; assert(cond); }

namespace ginkgo {

//! Organism memory manager default block size
const unsigned long ORGANISM_MEMORY_MANAGER_DEFAULT_BLOCK_SIZE = 1000;


//! The maximum number of fitness factors in the system. A "fitness factor"
//! represents an environmental variable and a corresponding genotypic variable
//! that together contribute to the fitness of an individual organism.
const unsigned MAX_FITNESS_TRAITS = 10;

//! The number of neutral diploid loci to track
const unsigned NUM_NEUTRAL_DIPLOID_loci = 10;

//! The value type of a fitness factor (both environmental and genotypic.
typedef float                   FitnessTraitType;

//! The number of dimensions or "slots" reserved for fitness assessment. The
//! actual number may be less than this.
typedef FitnessTraitType        FitnessTraits[MAX_FITNESS_TRAITS];

//! The global selection strength: the weighted least-squares distance of
//! a trait and its environmental optimum will be multiplied by this;
//! a value of 0 means no selection, while higher numbers increase the strength
//! of the selection.
const float DEFAULT_GLOBAL_SELECTION_STRENGTH = 1.0;

//! The units for referencing cells on the landscape.
#if defined(MINI)
    typedef unsigned char          CellIndexType;
    const CellIndexType MAX_LANDSCAPE_SIZE = static_cast<CellIndexType>(UCHAR_MAX);
#else
    typedef unsigned int           CellIndexType;
    const CellIndexType MAX_LANDSCAPE_SIZE = static_cast<CellIndexType>(UINT_MAX);
#endif

//! The carrying capacity type.
typedef unsigned int            PopulationCountType;

//! The movement economics type.
typedef int                     MovementCountType;

//! Numbers of generations.
typedef unsigned int            GenerationCountType;

/**
 * Landscape size too large error.
 */
class LandscapeSizeError : public std::runtime_error {
    public:
        LandscapeSizeError(const char * msg) : std::runtime_error(msg) {}
        LandscapeSizeError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

//! Landscape grid orientation: is (0,0) the upper-left? If not, it will be the lower-left
const bool GRID_ORIGIN_DEFAULT_UPPER_LEFT = true;

//! for vector of potential movement credits/capacity probabilities
const MovementCountType MAX_MOVEMENT_COUNT = 100;

} // ginkgo namespace

#endif
