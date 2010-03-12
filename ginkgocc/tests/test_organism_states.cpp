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

#include <iostream>
#include <cassert>
#include "../organism.hpp"
#include "../population.hpp"
#include "../species.hpp"

int main(int, char* []) {

    ginkgo::Organism m;
    m.init(ginkgo::Organism::Male);
    assert(!m.is_expired());
    assert(!m.is_female());
    assert(m.is_male());
    m.set_expired();
    assert(m.is_expired());
    assert(!m.is_female());
    assert(m.is_male());
    m.set_unexpired();
    assert(!m.is_expired());
    assert(!m.is_female());
    assert(m.is_male());

    ginkgo::Organism f;
    f.init(ginkgo::Organism::Female);
    assert(!f.is_expired());
    assert(f.is_female());
    assert(!f.is_male());
    f.set_expired();
    assert(f.is_expired());
    assert(f.is_female());
    assert(!f.is_male());
    f.set_unexpired();
    assert(!f.is_expired());
    assert(f.is_female());
    assert(!f.is_male());

    return 0;
}
