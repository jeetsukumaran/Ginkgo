///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Phylogeographical Evolution Simulator.
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

#include "landscape.hpp"

#include <iostream>
#include <iomanip>
#include "textutil.hpp"
#include "ginkgo_defs.hpp"

using namespace ginkgo;

///////////////////////////////////////////////////////////////////////////////
// Landscape

// --- lifecycle and assignment ---

Landscape::Landscape()
    : size_x_(0),
      size_y_(0),
      size_(0),
      max_x_(0),
      max_y_(0),
      max_size_(0),
      origin_upper_left_(GRID_ORIGIN_DEFAULT_UPPER_LEFT),
      rng_(RandomNumberGenerator::get_instance()) {
    this->size_x_ = 0;
    this->size_y_ = 0;
    this->size_ = 0;
}

// clean up cells
Landscape::~Landscape() {
    for (std::vector<Cell*>::iterator cell = this->cells_.begin();
            cell != this->cells_.end();
            ++cell) {
        delete *cell;
    }
}

// --- initialization and set up ---

void Landscape::generate(CellIndexType size_x, CellIndexType size_y, unsigned num_fitness_traits) {
    this->size_x_ = size_x;
    this->size_y_ = size_y;
    this->max_x_ = this->size_x_ - 1;
    this->max_y_ = this->size_y_ - 1;
    this->size_ = size_x * size_y;
    this->max_size_ = this->size_ - 1;
    this->cells_.reserve(this->size_);
    for (CellIndexType i = 0; i < this->size_; ++i) {
        CellIndexType x = this->index_to_x(i);
        CellIndexType y = this->index_to_y(i);
        Cell* cell = new Cell(i, x, y, num_fitness_traits, *this);
        this->cells_.push_back(cell);
    }
}

// --- migration and movement ---

void Landscape::clear_migrants() {
    this->migrants_.clear();
}

void Landscape::process_migrants() {
    for (MigrationEvents::iterator m = this->migrants_.begin();
            m != this->migrants_.end();
            ++m) {
        m->first->set_fitness(-1);   // invalidate cached fitness
        this->cells_.at(m->second)->add_organism(m->first);
    }
    this->migrants_.clear();
}

// --- sampling for tree building ---

void Landscape::sample_organisms(Species * sp_ptr,
        PopulationCountType num_organisms_per_cell,
        const std::set<CellIndexType>& cell_indexes,
        std::vector<const Organism *>& samples) {

    std::set<CellIndexType>::const_iterator start;
    std::set<CellIndexType>::const_iterator end;
    std::set<CellIndexType> all_cells;
    if (cell_indexes.size() == 0) {
        for (CellIndexType i = 0; i < this->cells_.size(); ++i) {
            all_cells.insert(i);
        }
        start = all_cells.begin();
        end = all_cells.end();
    } else {
        start = cell_indexes.begin();
        end = cell_indexes.end();
    }
    for (std::set<CellIndexType>::const_iterator ci = start;
            ci != end;
            ++ci) {
        assert(static_cast<PopulationCountType>(*ci) < this->cells_.size());
        this->cells_[*ci]->sample_organisms(sp_ptr, samples, num_organisms_per_cell);
    }
}

void Landscape::count_organisms(Species * sp_ptr, std::vector<PopulationCountType>& counts) const {
    counts.reserve(this->cells_.size());
    for (std::vector<Cell *>::const_iterator ci = this->cells_.begin(); ci != this->cells_.end(); ++ci) {
        counts.push_back((*ci)->num_organisms(sp_ptr));
    }
}

LandscapeOrganismProvenanceProportions Landscape::get_organism_provenances(Species * sp_ptr) const {
    LandscapeOrganismProvenanceProportions landscape_migrant_freqs;
    landscape_migrant_freqs.reserve(this->cells_.size());
    for (std::vector<Cell *>::const_iterator ci = this->cells_.begin();
            ci != this->cells_.end();
            ++ci) {
        Cell * cip = *ci;
        PopulationCountType total_pop_size = cip->num_organisms(sp_ptr);
        if (total_pop_size == 0) {
            CellOrganismProvenanceProportions cell_migrant_freqs(this->cells_.size(), 0.0);
            landscape_migrant_freqs.push_back(cell_migrant_freqs);
        } else {
            std::vector<float> cell_migrant_freqs;
            cell_migrant_freqs.reserve(this->cells_.size());
            OrganismProvenances provenances = (*ci)->get_organism_provenances(sp_ptr);
            for (std::vector<Cell *>::const_iterator cj = this->cells_.begin();
                    cj != this->cells_.end();
                    ++cj) {
                PopulationCountType sub_count = provenances.get_count((*cj)->get_index());
                cell_migrant_freqs.push_back(static_cast<float>(sub_count)/total_pop_size);
            }
        }
    }
    return landscape_migrant_freqs;
}

// --- debug output ---

void Landscape::debug_dump_cell_xy(std::ostream& out) {
    CellIndexType i = 0;
    unsigned long max_width_x = static_cast<unsigned long>(log10(static_cast<double>(this->size_x_))) + 1;
    unsigned long max_width_y = static_cast<unsigned long>(log10(static_cast<double>(this->size_y_))) + 1;
    for (std::vector<Cell*>::iterator ci = this->cells_.begin();
            ci != this->cells_.end();
            ++ci, ++i) {
        CellIndexType x = this->index_to_x(i);
        CellIndexType y = this->index_to_y(i);
        if (x == 0) {
            out << std::endl;
            if ( (y > 0) && ((y % 5) == 0) ) {
                out << std::endl;
            }
        } else if ( (x % 5) == 0 ) {
            out << "  ";
        }
        out << std::setfill('0') << std::setw(max_width_x) << (*ci)->get_x();
        out << ",";
        out << std::setfill('0') << std::setw(max_width_y) << (*ci)->get_y();
        out << " ";
    }
    out << std::endl << std::endl;
    out << std::setfill(' '); // reset
}

void Landscape::debug_dump_cell_indexes(std::ostream& out) {
    CellIndexType i = 0;
    unsigned long max_width = static_cast<unsigned long>(log10(static_cast<double>(this->size_))) + 1;
    for (std::vector<Cell*>::iterator ci = this->cells_.begin();
            ci != this->cells_.end();
            ++ci, ++i) {
        CellIndexType x = this->index_to_x(i);
        CellIndexType y = this->index_to_y(i);
        if (x == 0) {
            out << std::endl;
            if ( (y > 0) && ((y % 5) == 0) ) {
                out << std::endl;
            }
        } else if ( (x % 5) == 0 ) {
            out << "  ";
        }
        out << std::setfill('0') << std::setw(max_width) << i;
        out << " ";
    }
    out << std::endl << std::endl;
    out << std::setfill(' '); // reset
}


void Landscape::debug_dump_carrying_capacity(std::ostream& out) {
    CellIndexType i = 0;
    for (std::vector<Cell*>::iterator ci = this->cells_.begin();
            ci != this->cells_.end();
            ++ci, ++i) {
        CellIndexType x = this->index_to_x(i);
        CellIndexType y = this->index_to_y(i);
        if (x == 0) {
            out << std::endl;
            if ( (y > 0) && ((y % 5) == 0) ) {
                out << std::endl;
            }
        } else if ( (x % 5) == 0 ) {
            out << "  ";
        }
        out << this->cells_[i]->get_carrying_capacity() << " ";
    }
    out << std::endl << std::endl;
}

// Landscape
///////////////////////////////////////////////////////////////////////////////
