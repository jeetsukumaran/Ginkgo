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
#include <string>
#include <map>
#include <utility>
#include <ctime>
#include <fstream>
#include <iomanip>

#include "asciigrid.hpp"
#include "world.hpp"
#include "filesys.hpp"

using namespace gingko;

// constructor
World::World() 
    : species_(),
      rng_(),
      landscape_(species_, rng_),
      is_log_to_screen_(true) {
    this->current_generation_ = 0;    
}

// constructor
World::World(unsigned long seed) 
    : species_(),
      rng_(seed),
      landscape_(species_, rng_),
      is_log_to_screen_(true) {
    this->current_generation_ = 0;    
}    

// clean up species pool
World::~World() {
    for (SpeciesByLabel::iterator sp = this->species_.begin();
            sp != this->species_.end();
            ++sp) {
        assert (sp->second != NULL);            
        delete sp->second;
        sp->second = NULL;
        this->species_.erase(sp);
    }            
}

// --- initialization and set up ---

// Creates a new landscape.
void World::generate_landscape(CellIndexType size_x, CellIndexType size_y) {
    this->landscape_.generate(size_x, size_y, this->num_fitness_factors_);
}

// Adds a new species definition to this world.
Species& World::new_species(const std::string& label) {
    Species* sp = new Species(label, 
                              this->num_fitness_factors_, 
                              this->rng_);
    this->species_.insert(std::make_pair(std::string(label), sp));
    std::vector<int> default_movement_costs(this->landscape_.size(), 1);
    sp->set_movement_costs(default_movement_costs);
    return *sp;
}

// Populates the cell at (x,y) with organisms of the given species.
void World::seed_population(CellIndexType x, CellIndexType y, const std::string& species_label, unsigned long size) {
    assert(this->species_.find(species_label) != this->species_.end());
    this->landscape_.at(x, y).generate_new_organisms(this->species_[species_label], size);
}

// Populates the cell cell_index with organisms of the given species.
void World::seed_population(CellIndexType cell_index, const std::string& species_label, unsigned long size) {
    assert(this->species_.find(species_label) != this->species_.end());
    this->landscape_.at(cell_index).generate_new_organisms(this->species_[species_label], size);
}

// --- event handlers ---

WorldSettings& World::add_world_settings(unsigned long generation, const WorldSettings& world_settings) {
    this->world_settings_[generation] = world_settings;
    return this->world_settings_[generation];
}

// --- simulation cycles ---

void World::cycle() {

// Results in inflated population: the migrants get distributed after the 
// competition phase, resulting in a artificially (>> carrying capacity) 
// boosted population when entering the next generation's reproduction phase.
// This leads to a standing population at the end of each generation sometimes
// an order or more of magnitude above the carrying capacity of the cell.
//     for (CellIndexType i = this->landscape_.size()-1; i >= 0; --i) {
//         this->landscape_[i].reproduction(); 
//         this->landscape_[i].migration();
//         this->landscape_[i].survival();
//         this->landscape_[i].competition();
//     }
//     this->landscape_.process_migrants();
    std::ostringstream gen;
    gen << "Generation " << this->current_generation_ << " begun.";
    this->log_info(gen.str());
    this->log_info("Reproduction/migration phase.");
    for (CellIndexType i = this->landscape_.size()-1; i >= 0; --i) {
        this->landscape_[i].reproduction(); 
        this->landscape_[i].migration();
    }
    this->log_info("Processing migrants.");
    this->landscape_.process_migrants();
    this->log_info("Survival/competition phase.");
    for (CellIndexType i = this->landscape_.size()-1; i >= 0; --i) {    
        this->landscape_[i].survival();
        this->landscape_[i].competition();        
    }    
    this->log_info("Generation ended.");
    ++this->current_generation_;
}

void World::run() {    
    this->open_logs();
    this->log_extrasim_info("Starting simulation.");
    while (this->current_generation_ < this->generations_to_run_) {
        std::map<unsigned long, WorldSettings>::iterator wi = this->world_settings_.find(this->current_generation_);
        if (wi != this->world_settings_.end()) {
            if (wi->second.carrying_capacity.size() != 0) {
                this->log_info("Setting carrying capacity: \"" + wi->second.carrying_capacity + "\".");
                asciigrid::AsciiGrid grid(wi->second.carrying_capacity);
                this->landscape_.set_carrying_capacities(grid.get_cell_values());
            }
            if (wi->second.environments.size() != 0) {
                for (std::map<unsigned, std::string>::iterator ei = wi->second.environments.begin();
                     ei != wi->second.environments.end();
                     ++ei) {
                    std::ostringstream msg;
                    msg << "Setting environmental variable " <<  ei->first+1 <<  ": \"" <<  ei->second <<  "\"";
                    this->log_info(msg.str());
                    asciigrid::AsciiGrid grid(ei->second);
                    this->landscape_.set_environment(ei->first, grid.get_cell_values());                    
                }
            }            
            if (wi->second.movement_costs.size() != 0) {

            }
            if (wi->second.samples.size() != 0) {
                // process
            }
        }
        this->cycle();        
    }
    this->log_extrasim_info("Ending simulation.");
}

// --- logging and output ---

void World::open_ofstream(std::ofstream& out, const std::string& fpath) {
    std::string full_fpath = filesys::compose_path(this->output_dir_, fpath);
    out.open(full_fpath.c_str());
    if (not out) {
        throw WorldIOError("cannot open log file \"" + full_fpath + "\" for output");
    }
}

void World::open_logs() {    
    if (this->label_.size() == 0) {
        this->label_ = "world";
    }
    if (not this->infos_.is_open()) {
        this->open_ofstream(this->infos_, this->label_ + ".gingko.out.log");
    }
    if (not this->errs_.is_open()) {
        this->open_ofstream(this->errs_, this->label_ + ".gingko.err.log");
    }    
}

void World::close_logs() {
    if (this->infos_.is_open()) {
        this->infos_.close();
    }
    if (this->errs_.is_open()) {
        this->errs_.close();
    }    
}

std::string World::get_timestamp() {    
    time_t rawtime;    
    time ( &rawtime );
    struct tm * timeinfo = localtime ( &rawtime );
    char buffer[80];
    strftime (buffer,80,"%Y-%m-%d %H:%M:%S",timeinfo);
    return std::string(buffer);
}

std::string World::get_time_gen_stamp() {    
    std::ostringstream outs;
    outs << "[";
    outs << this->get_timestamp();
    outs << "]";    
    outs << " ";
    outs << std::setw(8) << std::setfill('0');
    outs << this->current_generation_;   
    return outs.str();
}

void World::log_extrasim_info(const std::string& message) {
    assert(this->infos_);
    std::ostringstream outs;
    outs << "[";
    outs << this->get_timestamp();
    outs << "]";
    outs << "         ";    
    if (this->is_log_to_screen_) {
        std::cout << outs.str() << " " << message << std::endl;
    }
    this->infos_ << outs.str() << " " << message << std::endl;
}

void World::log_info(const std::string& message) {
    assert(this->infos_);
    std::string ts = this->get_time_gen_stamp();
    if (this->is_log_to_screen_) {
        std::cout << ts << " " << message << std::endl;
    }
    this->infos_ << ts << " " << message << std::endl;
}

void World::log_error(const std::string& message) {
    assert(this->errs_);
    assert(this->infos_);    
    std::string ts = this->get_time_gen_stamp();
    if (this->is_log_to_screen_) {
        std::cerr << ts << " ERROR: " << message << std::endl;
    }
    this->errs_ << ts << " ERROR: " << message << std::endl;
    this->infos_ << ts << " ERROR: " << message << std::endl;
}


