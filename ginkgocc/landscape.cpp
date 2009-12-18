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

#include "landscape.hpp"

#include <iostream>
#include <iomanip>
#include "textutil.hpp"

using namespace ginkgo;

///////////////////////////////////////////////////////////////////////////////
// Landscape

// --- lifecycle and assignment ---

Landscape::Landscape(const SpeciesByLabel& species, RandomNumberGenerator& rng)
    : species_(species),
      rng_(rng) {
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
    this->size_ = size_x * size_y;
    this->cells_.reserve(this->size_);
    for (CellIndexType y = 0, index = 0; y < size_y; ++y) {
        for (CellIndexType x = 0; x < size_x; ++x, ++index) {
            Cell* cell = new Cell(index, x, y, num_fitness_traits, *this, this->species_, this->rng_);
            this->cells_.push_back(cell);
        }
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
        m->first.set_fitness(-1);   // invalidate cached fitness
        this->cells_.at(m->second)->insert_organism(m->first);
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

// --- debug output ---

void Landscape::debug_dump_cell_xy(std::ostream& out) {
    CellIndexType i = 0;
    unsigned long max_width_x = static_cast<unsigned long>(log10(static_cast<double>(this->size_x_))) + 1;
    unsigned long max_width_y = static_cast<unsigned long>(log10(static_cast<double>(this->size_y_))) + 1;
    for (std::vector<Cell*>::iterator ci = this->cells_.begin();
            ci != this->cells_.end();
            ++ci, ++i) {
        CellIndexType x = this->index_to_x(i);
        if (x == 0) {
            out << std::endl;
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
        if (x == 0) {
            out << std::endl;
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
        if (x == 0) {
            out << std::endl;
        }
        out << this->cells_[i]->get_carrying_capacity() << " ";
    }
    out << std::endl << std::endl;
}

// Landscape
///////////////////////////////////////////////////////////////////////////////
