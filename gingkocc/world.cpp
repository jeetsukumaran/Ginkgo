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
#include "tree.hpp"
#include "convert.hpp"

using namespace gingko;

// constructor
World::World() 
    : species_(),
      rng_(),
      landscape_(species_, rng_),
      num_fitness_factors_(1),
      fitness_factor_grain_(1),
      generations_to_run_(0),
      current_generation_(0),
      log_frequency_(10),
      is_log_to_screen_(true),
      is_produce_final_output_(true) {
    this->current_generation_ = 0;    
}

// constructor
World::World(unsigned long seed) 
    : species_(),
      rng_(seed),
      landscape_(species_, rng_),
      num_fitness_factors_(1),
      fitness_factor_grain_(1),
      generations_to_run_(0),
      current_generation_(0),
      log_frequency_(10),      
      is_log_to_screen_(true),
      is_produce_final_output_(true) {
    this->current_generation_ = 0;    
}    

// clean up species pool
World::~World() {
    for (SpeciesByLabel::iterator sp = this->species_.begin();
            sp != this->species_.end();
            ++sp) {
        assert (sp->second != NULL);            
        delete sp->second;
//         sp->second = NULL;
//         this->species_.erase(sp);
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
                              this->fitness_factor_grain_,
                              this->rng_);
    this->species_.insert(std::make_pair(std::string(label), sp));
    std::vector<long> default_movement_costs(this->landscape_.size(), 1);
    sp->set_movement_costs(default_movement_costs);
    return *sp;
}

// Populates the cell at (x,y) with organisms of the given species.
// void World::seed_population(CellIndexType x, CellIndexType y, const std::string& species_label, unsigned long size) {
//     assert(this->species_.find(species_label) != this->species_.end());
//     // this->landscape_.at(x, y).generate_new_organisms(this->species_[species_label], size);
// }

// Populates the cell cell_index with organisms of the given species.
void World::generate_seed_population(CellIndexType cell_index, 
        Species * species_ptr, 
        unsigned long pop_size,
        unsigned long ancestral_pop_size,
        unsigned long ancestral_generations) {
    std::ostringstream pre_msg;
    pre_msg << "[Generation " << this->current_generation_ << "] ";
    pre_msg << "Bootstrapping seed population of species " << species_ptr->get_label() << ": ";   
    pre_msg << ancestral_pop_size << " individuals for " << ancestral_generations << " generations.";
    this->log_info(pre_msg.str());
        
    this->landscape_.at(cell_index).generate_new_population(species_ptr,
            pop_size,
            ancestral_pop_size,
            ancestral_generations);
            
    std::ostringstream post_msg;
    post_msg << "[Generation " << this->current_generation_ << "] ";
    post_msg << "Seeding population ";
    post_msg << "of species " << species_ptr->get_label()  << " ";
    post_msg << "in (" << this->landscape_.index_to_x(cell_index) <<  "," << this->landscape_.index_to_y(cell_index) << "): ";
    post_msg << pop_size << " individuals drawn from an ancestral population of " << ancestral_pop_size << " ";
    post_msg << "after " << ancestral_generations << " generations.";
    this->log_info(post_msg.str());
}

// --- event handlers ---

void World::add_world_settings(unsigned long generation, const WorldSettings& world_settings) {
    this->world_settings_[generation] = world_settings;
}

void World::add_dispersal_event(unsigned long generation, const DispersalEvent& dispersal_event) {
    this->dispersal_events_.insert(std::make_pair(generation, dispersal_event));
}


void World::add_tree_sampling(unsigned long generation, const SamplingRegime& sampling_regime) {
    this->tree_samples_.insert(std::make_pair(generation, sampling_regime));
}

void World::add_occurrence_sampling(unsigned long generation, Species * species_ptr) {
    this->occurrence_samples_.insert(std::make_pair(generation, species_ptr));
}

void World::add_seed_population(CellIndexType cell_index, 
        Species * species_ptr,
        unsigned long pop_size,
        unsigned long ancestral_pop_size,
        unsigned long ancestral_generations) {
    this->seed_populations_.push_back( SeedPopulation(cell_index, species_ptr, pop_size, ancestral_pop_size, ancestral_generations) );           
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

    if ( this->current_generation_ % this->log_frequency_ == 0) {
        std::ostringstream gen;
        gen << "Generation " << this->current_generation_ << " life-cycle running.";
        this->log_info(gen.str());
    }
//     this->log_detail("Reproduction/migration phase.");
    for (CellIndexType i = 0; i < this->landscape_.size(); ++i) {
        this->landscape_[i].reproduction(); 
        this->landscape_[i].migration();
    }
//     this->log_detail("Processing migrants.");
    this->landscape_.process_migrants();
//     this->log_detail("Survival/competition phase.");
    for (CellIndexType i = 0; i < this->landscape_.size(); ++i) {    
        this->landscape_[i].survival();
        this->landscape_[i].competition();        
    }    
//     this->log_detail("Generation life-cycle completed.");
    ++this->current_generation_;
}

void World::run() {    
    this->open_logs();
    this->log_configuration();
        
    this->log_info("Starting simulation.");
    
    // startup    
    for (std::vector<SeedPopulation>::iterator spi = this->seed_populations_.begin();
            spi != this->seed_populations_.end();
            ++spi) {
        SeedPopulation& sp = *spi;
        this->generate_seed_population(sp.cell_index, 
            sp.species_ptr, 
            sp.pop_size, 
            sp.ancestral_pop_size,
            sp.ancestral_generations);
    }            
    
    while (this->current_generation_ <= this->generations_to_run_) {
        
        // clear organism labels
        for (std::map<std::string, Species *>::iterator spi = this->species_.begin(); 
                spi != this->species_.end(); 
                ++spi) {
            (spi->second)->clear_organism_labels();
        }
        
        // clear output filename stems
        this->output_filenames_.clear();
        
        // build trees requested in this generation
        this->process_tree_samplings();
                
        // save occurrence data requested in this generation
        this->process_occurrence_samplings();
        
        // process world changes
        this->process_world_settings();
        
        // process dispersal events
        this->process_dispersal_events();        
        
        // run the life cycle
        this->cycle();        
    }
    
    if (this->is_produce_final_output_) {
        this->log_info("Saving final set of occurrences and trees for all species.");
        std::set<CellIndexType> cell_indexes;
        for (CellIndexType i = 0; i < this->landscape_.size(); ++i) {
            cell_indexes.insert(i);        
        }
        for (std::map<std::string, Species *>::iterator spi = this->species_.begin(); 
                spi != this->species_.end(); 
                ++spi) {
            this->save_occurrences(spi->second);
            this->save_trees(spi->second, 0, cell_indexes, "COMPLETE");
        }
    }        
    
    this->log_info("Ending simulation.");
}

void World::process_world_settings() {
    std::map<unsigned long, WorldSettings>::iterator wi = this->world_settings_.find(this->current_generation_);    
    if (wi == this->world_settings_.end()) {
        return;
    }    
    if (wi->second.carrying_capacity.size() != 0) {
        this->log_info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Setting carrying capacity: \"" + wi->second.carrying_capacity + "\".");
        asciigrid::AsciiGrid grid(wi->second.carrying_capacity);
        this->landscape_.set_carrying_capacities(grid.get_cell_values());
    }
    if (wi->second.environments.size() != 0) {
        for (std::map<unsigned, std::string>::iterator ei = wi->second.environments.begin();
                 ei != wi->second.environments.end();
                 ++ei) {
            std::ostringstream msg;
            msg << "[Generation " << this->current_generation_ << "] Setting environmental variable " <<  ei->first+1 <<  ": \"" <<  ei->second <<  "\"";
            this->log_info(msg.str());
            asciigrid::AsciiGrid grid(ei->second);
            this->landscape_.set_environment(ei->first, grid.get_cell_values());                    
        }
    }            
    if (wi->second.movement_costs.size() != 0) {
        for (std::map<Species *, std::string>::iterator mi = wi->second.movement_costs.begin();
                 mi != wi->second.movement_costs.end();
                 ++mi) {
            std::ostringstream msg;
            msg << "[Generation " << this->current_generation_ << "] Setting movement costs for species " <<  mi->first->get_label() <<  ": \"" <<  mi->second <<  "\"";
            this->log_info(msg.str());
            asciigrid::AsciiGrid grid(mi->second);
            this->set_species_movement_costs(mi->first, grid.get_cell_values());                    
        }
    }
}

void World::process_dispersal_events() {
    typedef std::multimap<unsigned long, DispersalEvent> gen_disp_t;
    typedef std::pair<gen_disp_t::iterator, gen_disp_t::iterator> gen_disp_iter_pair_t;    
    gen_disp_iter_pair_t this_gen_dispersals = this->dispersal_events_.equal_range(this->current_generation_);
    for (gen_disp_t::iterator di = this_gen_dispersals.first; di != this_gen_dispersals.second; ++di) {            
        DispersalEvent& de = di->second;
        std::ostringstream msg;
        msg << "[Generation " << this->current_generation_ << "] Dispersal";
        if (de.species_ptr != NULL) {
            msg << " of " << de.species_ptr->get_label();
        }
        msg << " from (" << this->landscape_.index_to_x(de.source) << "," << this->landscape_.index_to_y(de.source) << ")";
        msg << " to (" << this->landscape_.index_to_x(de.destination) << "," << this->landscape_.index_to_y(de.destination) << ")";
        msg << " with probability " << de.probability << ": ";
        
        std::vector<const Organism *> organisms;
        this->landscape_[de.source].sample_organisms(de.species_ptr, organisms, de.num_organisms);
        this->landscape_.clear_migrants();
        unsigned long num_males = 0;
        unsigned long num_females = 0;
        for (std::vector<const Organism *>::iterator oi = organisms.begin();  oi != organisms.end(); ++oi) {
            if (this->rng().uniform_01() > de.probability) {
                Organism* og = const_cast<Organism*>(*oi);
                if (og->is_male()) {
                    ++num_males;
                } else {
                    ++num_females;
                }
                this->landscape_.add_migrant(*og, de.destination);
                og->set_expired(true);
            }                
        }
        msg << " " << num_females << " females and " << num_males << " males.";
        this->landscape_[de.source].purge_expired_organisms();
        this->landscape_.process_migrants();
        this->log_info(msg.str());              
    }
}

void World::process_tree_samplings() {
    typedef std::multimap<unsigned long, SamplingRegime> gen_sample_t;
    typedef std::pair<gen_sample_t::iterator, gen_sample_t::iterator> gen_sample_iter_pair_t;    
    gen_sample_iter_pair_t this_gen_samples = this->tree_samples_.equal_range(this->current_generation_);
    for (gen_sample_t::iterator i = this_gen_samples.first; i != this_gen_samples.second; ++i) {
        this->save_trees(i->second.species_ptr, i->second.num_organisms_per_cell, i->second.cell_indexes, i->second.label);
    }
}

void World::process_occurrence_samplings() {
    typedef std::multimap<unsigned long, Species *> gen_sample_t;
    typedef std::pair<gen_sample_t::iterator, gen_sample_t::iterator> gen_sample_iter_pair_t;    
    gen_sample_iter_pair_t this_gen_samples = this->occurrence_samples_.equal_range(this->current_generation_);
    for (gen_sample_t::iterator i = this_gen_samples.first; i != this_gen_samples.second; ++i) {
        this->save_occurrences(i->second);
    }
}  

// --- logging and output ---

void World::save_occurrences(Species * species_ptr ) {
    this->log_info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Saving occurrence data for species " + species_ptr->get_label() + ".");
    std::vector<long> counts;
    this->landscape_.count_organisms(species_ptr, counts);    
    std::ofstream occs;
    this->open_ofstream(occs,
        this->compose_output_filename(species_ptr->get_label(), "occurrences", "grd"));  
    asciigrid::write_grid(counts, this->landscape_.size_x(), this->landscape_.size_y(), occs);
}

void World::write_nexus_header(Species * sp_ptr,
        const std::vector<const Organism *>& organisms,
        std::ostream& out,
        bool add_diploid_allele_extensions) {
    out << "#NEXUS\n\n";
    out << "BEGIN TAXA;\n";
    if (add_diploid_allele_extensions) {
        out << "    DIMENSIONS NTAX=" << organisms.size()*2 << ";\n";
    } else {
        out << "    DIMENSIONS NTAX=" << organisms.size() << ";\n";
    }        
    out << "    TAXLABELS\n";
    for (std::vector<const Organism *>::const_iterator oi = organisms.begin();
            oi != organisms.end();
            ++oi) {
        if (add_diploid_allele_extensions) {
            out << "        " << sp_ptr->get_organism_label(**oi) << "_a1" << "\n";
            out << "        " << sp_ptr->get_organism_label(**oi) << "_a2" << "\n";
        } else {
            out << "        " << sp_ptr->get_organism_label(**oi) << "\n";
        }            
    }
    out << "    ;\n";
    out << "END;\n\n";
}        

void World::write_tree(Tree& tree, const std::string& species_label, unsigned long num_taxa, std::ostream& out) {
    try {
        tree.write_newick_tree(out);
    } catch (const TreeStructureMissingRootError& e) {
        std::ostringstream msg;
        msg << "tree for sample of organisms of species " << species_label;
        msg << " in generation " << this->current_generation_;
        msg << " could not be built due to missing root node";
        msg << " (organism sample size = " << num_taxa << ")";
        this->log_error(msg.str());
    } catch (const TreeStructureMultipleRootError& e) {
        std::ostringstream msg;
        msg << "tree for sample of organisms of species " << species_label;
        msg << " in generation " << this->current_generation_;
        msg << " could not be built due to multiple root nodes";
        msg << " (organism sample size = " << num_taxa << ")";
        this->log_error(msg.str());    
    }
}

void World::write_haploid_tree(Species * sp_ptr,
                const std::vector<const Organism *>& organisms,
                std::ostream& out) {
    Tree tree(&this->landscape_);
    
    // build tree
    for (std::vector<const Organism *>::const_iterator oi = organisms.begin();
            oi != organisms.end();
            ++oi) {
        tree.add_leaf((*oi)->get_haploid_node(), &sp_ptr->get_organism_label(**oi));
    }
    
    // write tree
    this->write_nexus_header(sp_ptr, organisms, out);
    out << "BEGIN TREES;\n";
    out << "    TREE HaploidLocus = [&R] ";
    this->write_tree(tree, sp_ptr->get_label(), organisms.size(), out);
    out << ";\n";
    out << "END;\n\n";
}

void World::write_diploid2_trees(Species * sp_ptr,
                const std::vector<const Organism *>& organisms,
                std::ostream& out) {
                
    this->write_nexus_header(sp_ptr, organisms, out, true);
    out << std::setfill(' '); // reset
    out << "BEGIN TREES;\n";
    for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_LOCII; ++i) {
        Tree tree(&this->landscape_);
        std::string allele1;
        std::string allele2;
        for (std::vector<const Organism *>::const_iterator oi = organisms.begin();
                oi != organisms.end();
                ++oi) {
            allele1 = sp_ptr->get_organism_label(**oi) + "_a1";
            allele2 = sp_ptr->get_organism_label(**oi) + "_a2";
            tree.add_leaf((*oi)->get_diploid_node1(i), &allele1);
            tree.add_leaf((*oi)->get_diploid_node2(i), &allele2);            
        }
        out << "    TREE DiploidLocus" << std::setw(2) << std::setfill('0') << i+1 << " = [&R] "; 
        this->write_tree(tree, sp_ptr->get_label(), organisms.size(), out);
        out << ";\n";
    }
    out << "END;\n\n";
    out << std::setfill(' '); // reset
}

void World::write_diploid1_trees(Species * sp_ptr,
        const std::vector<const Organism *>& organisms,
        std::ostream& out) {

    this->write_nexus_header(sp_ptr, organisms, out);
    out << std::setfill(' '); // reset
    out << "BEGIN TREES;\n";
    
    for (unsigned i = 0; i < NUM_NEUTRAL_DIPLOID_LOCII; ++i) {
        Tree tree(&this->landscape_);
        for (std::vector<const Organism *>::const_iterator oi = organisms.begin();
                oi != organisms.end();
                ++oi) {
            tree.add_leaf((*oi)->get_diploid_random_node(i, this->rng_), &sp_ptr->get_organism_label(**oi));          
        }
        out << "    TREE DiploidLocus" << std::setw(2) << std::setfill('0') << i+1 << " = [&R] "; 
        this->write_tree(tree, sp_ptr->get_label(), organisms.size(), out);
        out << ";\n";
    }
    out << "END;\n\n";        
    out << std::setfill(' '); // reset        
}        

void World::save_trees(Species * sp_ptr, 
                unsigned long num_organisms_per_cell, 
                const std::set<CellIndexType>& cell_indexes,
                const std::string& label) {
                
    std::ostringstream msg;
    msg << "[Generation " << this->current_generation_ << "] ";
    msg << "Sampling organisms of species " << sp_ptr->get_label() << " (";
    if (num_organisms_per_cell == 0) {
        msg << "all individuals per cell,";
    } else {
        msg << num_organisms_per_cell << " individuals per cell,";
    }
    if (cell_indexes.size() == 0) {
        msg << " all cells).";
    } else {
        msg << " from " << cell_indexes.size() << " cells).";
    }
    this->log_info(msg.str());
    
    std::vector<const Organism *> organisms;
    this->landscape_.sample_organisms(sp_ptr, num_organisms_per_cell, cell_indexes, organisms);
    std::string num_samples = convert::to_scalar<std::string>(organisms.size());
    
    if (organisms.size() == 0) {
        this->log_error("no organisms found in sample: aborting tree building");
        return;
    }
    
    this->log_info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Building trees for subsampled diploid locii alleles for " + num_samples + " organisms (" + num_samples + " leaves per tree).");  
    std::ofstream combined_trees;            
    this->open_ofstream(combined_trees,
        this->compose_output_filename(sp_ptr->get_label(), label, "diploid1.tre"));
    this->write_diploid1_trees(sp_ptr, organisms, combined_trees);          

    this->log_info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Building tree of haploid locus alleles for " + num_samples + " organisms (" + num_samples + " leaves per tree).");    
    std::ofstream haploid_trees;
        
    this->open_ofstream(haploid_trees, 
        this->compose_output_filename(sp_ptr->get_label(), label, "haploid.tre"));    
    this->write_haploid_tree(sp_ptr, organisms, haploid_trees);
    
    this->log_info("[Generation " + convert::to_scalar<std::string>(this->current_generation_) + "] Building trees for full diploid locii complement for " + num_samples + " organisms (" + convert::to_scalar<std::string>(organisms.size()*2) + " leaves per tree).");  
    std::ofstream diploid_trees;
    this->open_ofstream(diploid_trees,
        this->compose_output_filename(sp_ptr->get_label(), label, "diploid2.tre"));  
    this->write_diploid2_trees(sp_ptr, organisms, diploid_trees);
 
}

void World::open_ofstream(std::ofstream& out, const std::string& fpath) {
    std::string full_fpath = filesys::compose_path(this->output_dir_, fpath);
    out.open(full_fpath.c_str());
    if (not out) {
        throw WorldIOError("cannot open file \"" + full_fpath + "\" for output");
    }
}

std::string World::get_output_filename_stem() {
    if (this->output_filename_stem_.size() == 0) {
        if (this->label_.size() == 0) {
            this->output_filename_stem_ += "gingkorun";
        } else {
            this->output_filename_stem_ += this->label_;        
        }
        if (this->replicate_id_.size() > 0) {
            this->output_filename_stem_ += this->replicate_id_;
        }
    }
    return this->output_filename_stem_;           
}

std::string World::compose_output_filename(const std::string& species_label,
        const std::string& additional,
        const std::string& extension) {
    std::ostringstream f;
    f << this->get_output_filename_stem();
    f << "_G" << this->current_generation_;
    f << "_" << species_label;
    if (additional.size() > 0) {
        f << "_" << additional;
    }
    unsigned index = 0;
    std::string candidate_name = f.str() + "." + extension;
    std::set<std::string>::iterator fname_found = this->output_filenames_.find(candidate_name);
    while (fname_found != this->output_filenames_.end()) {
        ++index;
        candidate_name = f.str() + "-" + convert::to_scalar<std::string>(index) + "." + extension;
        fname_found = this->output_filenames_.find(candidate_name);
    }
    this->output_filenames_.insert(candidate_name);
    return candidate_name;
}        


void World::log_configuration() {
    std::ofstream out;
    this->open_ofstream(out, this->get_output_filename_stem() + ".conf.log");
    out <<  "GINGKO CONFIGURATION LOG " << this->get_timestamp() << std::endl;
        
    out << std::endl;
    out << "*** WORLD ***" << std::endl;
    out << "Label: " << this->label_ << std::endl;
    out << "Random seed: " << this->rng_.get_seed() << std::endl;
    out << "Fitness factors: " << this->num_fitness_factors_ << std::endl;
    out << "Fitness grain: " << this->fitness_factor_grain_ << std::endl;
    out << "Generations to run: " << this->generations_to_run_ << std::endl;
    out << "Output directory: " << this->output_dir_ << std::endl;
    out << "Replicate ID: " << this->replicate_id_ << std::endl; 

    out << std::endl;
    out << "*** LANDSCAPE ***" << std::endl;
    out << "Rows (X-dimension): " << this->landscape_.size_x() << std::endl;
    out << "Columns (Y-dimension): " << this->landscape_.size_y() << std::endl; 
    out << "Structure: " << std::endl;
    this->landscape_.debug_dump_structure(out);
    out << std::endl;
    out << "Default Carrying Capacity: " << std::endl;
    this->landscape_.debug_dump_carrying_capacity(out);

    out << std::endl;
    out << "*** LINEAGES ***" << std::endl;
    out << "(" << this->species_.size() << " lineages specified)" << std::endl;
    unsigned i = 0;
    for (SpeciesByLabel::iterator spi = this->species_.begin(); spi != this->species_.end(); ++spi) {
        i += 1;
        Species& lineage = *spi->second;
        out << std::endl;
        out << "\"" << lineage.get_label() << "\"" << std::endl;
        out << "     Fitness factors: " << lineage.get_num_fitness_factors() << std::endl;
        out << "     Fitness grain: " << lineage.get_fitness_factor_grain() << std::endl;
        out << "     Selection weights: ";
        std::vector<float> sw = lineage.get_selection_weights();
        for (std::vector<float>::iterator swi = sw.begin(); swi != sw.end(); ++swi) {
            if ((swi - sw.begin()) > 0) {
                out << " ";
            }
            out << *swi;
        }
        out << std::endl;
        out << "     Initial genotypic fitness factors: ";
        std::vector<FitnessFactorType> g = lineage.get_default_genotypic_fitness_factors();
        for (std::vector<FitnessFactorType>::iterator gi = g.begin(); gi != g.end(); ++gi) {
            if ((gi - g.begin()) > 0) {
                out << " ";
            }
            out << *gi;
        }           
        out << std::endl;
        out << "     Genotypic fitness mutation rate: " << lineage.get_mutation_rate() << std::endl;
        out << "     Genotypic fitness maximum mutation size: " << lineage.get_max_mutation_size() << std::endl;
        out << "     Fecundity: " <<  lineage.get_mean_reproductive_rate() << std::endl;
        out << "     Movement probability: " << lineage.get_movement_probability() << std::endl;
        out << "     Movement capacity: " << lineage.get_movement_capacity() << std::endl;
    }

    out << std::endl;
    out << "*** ENVIRONMENTS ***" << std::endl;
    out << "(" << this->world_settings_.size() << " generations with environmental changes specified)" << std::endl;
    
    for (std::map<unsigned long, WorldSettings>::iterator wi = this->world_settings_.begin(); wi != this->world_settings_.end(); ++wi) {
        out << "\n[ENVIRONMENT: GENERATION " << wi->first << "]" << std::endl;
        if (wi->second.carrying_capacity.size() != 0) {
            out << "    Carrying-capacity: \"" << wi->second.carrying_capacity << "\"" << std::endl;
        }
        if (wi->second.environments.size() != 0) {
            for (std::map<unsigned, std::string>::iterator ei = wi->second.environments.begin();
                     ei != wi->second.environments.end();
                     ++ei) {
                out << "    Environment factor #" << (ei->first + 1) << ": \"" << ei->second << "\"" << std::endl;                                                 
            }
        }            
        if (wi->second.movement_costs.size() != 0) {
            for (std::map<Species *, std::string>::iterator mi = wi->second.movement_costs.begin();
                     mi != wi->second.movement_costs.end();
                     ++mi) {
                out << "    Movement costs for lineage \"" << mi->first->get_label() << "\": \"" << mi->second << "\"" << std::endl;                   
            }
        }    
    }

    out << std::endl;
    out << "*** DISPERSALS ***" << std::endl;
    out << "(" << this->dispersal_events_.size() << " dispersal events specified)" << std::endl;    
    out << std::endl;    
    for (std::map<unsigned long, DispersalEvent>::iterator di = this->dispersal_events_.begin(); di != this->dispersal_events_.end(); ++di) {
        DispersalEvent& de = di->second;
        out << "    GENERATION " << di->first << ": " << "Dispersal";
        if (de.species_ptr != NULL) {
            out << " of " << de.species_ptr->get_label();
        }
        out << " from (" << this->landscape_.index_to_x(de.source) << "," << this->landscape_.index_to_y(de.source) << ")";
        out << " to (" << this->landscape_.index_to_x(de.destination) << "," << this->landscape_.index_to_y(de.destination) << ")";
        out << " with probability " << de.probability << "." << std::endl;
    }    
    
    out << std::endl;
    out << "*** OCCURRENCE SAMPLES ***";
    out << std::endl;    
    for (std::multimap<unsigned long, Species *>::iterator oci = this->occurrence_samples_.begin(); 
            oci != this->occurrence_samples_.end(); 
            ++oci) {
        out << "    GENERATION " << oci->first << ": " << oci->second->get_label() << std::endl;     
    }
    
    out << std::endl;
    out << "*** TREE SAMPLES ***";
    out << std::endl;    
    for (std::multimap<unsigned long, SamplingRegime>::iterator tci = this->tree_samples_.begin(); 
            tci != this->tree_samples_.end(); 
            ++tci) {
        SamplingRegime& sr = tci->second;            
        out << "    GENERATION " << tci->first << ":";
        if (sr.label.size() > 0) {
            out << " [sample \"" << sr.label << "\"]";
        }
        out << " Lineage \"" << sr.species_ptr->get_label() << "\"";
        if (sr.num_organisms_per_cell == 0) {
            out << ", all individuals";
        } else {
            out << ", " << sr.num_organisms_per_cell << " individuals";
        }
        if (sr.cell_indexes.size() == 0) {
            out << " from each cell, all cells";
        } else {
            out << " from each cell, from following cells: ";
            for (std::set<CellIndexType>::iterator ci = sr.cell_indexes.begin(); ci != sr.cell_indexes.end(); ++ci) {
                out << this->landscape_.index_to_x(*ci) << "," << this->landscape_.index_to_y(*ci) << " ";
            }
        }
        out << std::endl;
    }      
    out.close();
}

void World::open_logs() {    
    if (not this->infos_.is_open()) {
        this->open_ofstream(this->infos_, this->get_output_filename_stem() + ".out.log");
    }
    if (not this->errs_.is_open()) {
        this->open_ofstream(this->errs_, this->get_output_filename_stem() + ".err.log");
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
    return "[" + std::string(buffer) + "]";
}

void World::log_detail(const std::string& message) {
    assert(this->infos_);
    this->infos_ << this->get_timestamp();
    this->infos_ << " [Generation ";
    this->infos_ << this->current_generation_;
    this->infos_ << "] ";
    this->infos_ << message << std::endl;
}

void World::log_info(const std::string& message) {
    assert(this->infos_);
    std::string ts = this->get_timestamp();
    if (this->is_log_to_screen_) {
        std::cout << ts << " " << message << std::endl;
    }
    this->infos_ << ts << " " << message << std::endl;
}

void World::log_error(const std::string& message) {
    assert(this->errs_);
    assert(this->infos_);    
    std::string ts = this->get_timestamp();
    if (this->is_log_to_screen_) {
        std::cerr << ts << " ERROR: " << message << std::endl;
    }
    this->errs_ << ts << " ERROR: " << message << std::endl;
    this->infos_ << ts << " ERROR: " << message << std::endl;
}


