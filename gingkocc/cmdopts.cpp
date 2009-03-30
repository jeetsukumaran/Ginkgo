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

#include "cmdopts.hpp"
#include "textutils.hpp"
#include "filesys.hpp"

namespace gingko {

///////////////////////////////////////////////////////////////////////////////
// OptionArg

OptionArg::OptionArg(const char * help, const char * meta_var)
    : is_switch_(false),
      is_set_(false) {
    if (help != NULL) {
        this->help_ = help;
    }
    if (meta_var != NULL) {
        this->meta_var_ = meta_var;
    }
}
  
OptionArg::~OptionArg() {}
  

std::ostream& OptionArg::write_help(std::ostream& out) const {            
    std::string help_str;
    help_str += "  ";
    if (this->short_flag_.size() > 0) {               
        help_str += this->short_flag_;
        if (not this->is_switch_) {
            help_str += " ";
            if (this->meta_var_.size() == 0) {
                help_str += "VALUE";
            } else {
                help_str += this->meta_var_;
            }                        
        }
        if (this->long_flag_.size() > 0) {
            help_str += ", ";
        }
    }
    if (this->long_flag_.size() > 0) {            
        help_str += this->long_flag_;
        if (not this->is_switch_) {
            help_str += "=";
            if (this->meta_var_.size() == 0) {
                help_str += "VALUE";
            } else {
                help_str += this->meta_var_;
            } 
        }
    }            
    if (this->help_.size() > 0) {
        if (help_str.size() > CMDOPTS_OPTION_COL_WIDTH-2) {
            help_str += "\n";
        } else {
            while (help_str.size() < CMDOPTS_OPTION_COL_WIDTH) {
                help_str += " ";
            }                
        }
        std::string help_msg = this->help_;
        std::string::size_type defval = help_msg.find("%default");
        std::string replace_val = this->current_value_as_string();
        while (defval != std::string::npos) {
            help_msg.replace(defval, 8, replace_val.c_str());
            defval = help_msg.find("%default");
        }
        help_str += help_msg;
        std::string help_desc = textutils::textwrap(help_str, CMDOPTS_LINE_WIDTH, 0, CMDOPTS_OPTION_COL_WIDTH);
        help_str = help_desc;
    }                            
    out << help_str; 
    return out;
}

///////////////////////////////////////////////////////////////////////////////
// Specializations of TypedOptionArg

template <>
void TypedOptionArg<std::string>::process_value_string(const std::string& val_str) {
    *this->store_ = val_str;
    this->set_is_set(true);
}

///////////////////////////////////////////////////////////////////////////////
// OptionParser

OptionParser::OptionParser(const char * version,
        const char * description, 
        const char * usage)
    : show_help_(false),
      show_version_(false) {
    if (usage != NULL) {
        this->usage_.assign(usage);
    } else {
        this->usage_ = "%prog [options] [args]";
    }
    if (description != NULL) {
        this->description_.assign(description);
    }
    if (version != NULL) {
        this->version_.assign(version);
    }       
    this->version_option_ = this->add_switch(&this->show_version_, NULL, "--version", "show program's version number and exit");    
    this->help_option_ = this->add_switch(&this->show_help_, "-h", "--help",  "show this help message and exit");
}
    
OptionParser::~OptionParser() {
    for (std::vector<OptionArg *>::iterator oap = this->option_args_.begin();
            oap != this->option_args_.end();
            ++oap) {
        delete *oap;                    
    }                    
}

std::ostream& OptionParser::write_help(std::ostream& out, const char * progname) const {
    if (this->usage_.size() != 0) {
        if (progname == NULL) {
            progname = "PROGRAM";
        }
        std::string usage = "Usage: " + this->usage_;
        std::string::size_type pos = usage.find("%prog");
        while (pos != std::string::npos) {
            usage.replace(pos, 5, progname); 
            pos = usage.find("%prog");
        }
        out << usage << std::endl << std::endl;
    }
    if (this->description_.size() != 0) {
        out << textutils::textwrap(this->description_, CMDOPTS_LINE_WIDTH) << std::endl << std::endl;
    }
    out << "Options:" << std::endl;
    for (std::vector<OptionArg *>::const_iterator oa = this->option_args_.begin();
            oa != this->option_args_.end();
            ++oa) {
        (*oa)->write_help(out);
        out << std::endl;
    }
    return out;
}

void OptionParser::parse(int argc, char * argv[]) {

    for (int i = 1; i < argc; ++i) { 
        if (argv[i][0] == '-') {
            std::string arg_name;
            std::string arg_value;

            if (strncmp(argv[i], "--", 2) == 0) {
                bool parsing_name = true;
                for (char *a = argv[i]; *a; ++a) {
                    if (parsing_name) {
                        if (*a == '=') {
                            parsing_name = false;
                        } else {
                            arg_name += *a;
                        }
                    } else {
                        arg_value += *a;
                    }
                }                        
            } else if (argv[i][0] == '-') {
                std::string arg(argv[i]);
                if (arg.size() < 2) {
                    std::cerr << "unrecognized or incomplete option \"" << arg << "\"" << std::endl;
                    exit(1);
                }
                if (arg.size() == 2) {
                    arg_name = arg;
                } else {
                    arg_name = arg.substr(0, 2);
                    arg_value = arg.substr(2, arg.size());
                }
            }
            
            std::map< std::string, OptionArg * >::iterator oai = this->key_opt_map_.find(arg_name);
            if ( oai == this->key_opt_map_.end() ) {
                std::cerr << "unrecognized option \"" << arg_name << "\"" << std::endl;
                exit(1);
            }          
            
            OptionArg& oa = *(oai->second);
            
            if (not oa.is_switch()) {
                if (arg_value.size() == 0) {
                    if (i == argc-1) {
                        std::cerr << "expecting value for option \"" << arg_name << "\"" << std::endl;
                        exit(1);
                    } else {
                        arg_value = argv[i+1];
                        i += 1;
                    }                            
                }
                try {
                    oa.process_value_string(arg_value);                    
                } catch(OptionValueTypeError& e) {
                    std::cerr << "Invalid value passed to option " << arg_name << ": ";
                    std::cerr << "\"" << arg_value << "\"" << std::endl;                    
                    exit(1);
                }
            } else {
                TypedOptionArg<bool>* bool_opt = static_cast< TypedOptionArg<bool> *>(&oa);
                bool_opt->process_value(true);            
            }
            
        } else {
            this->pos_args_.push_back(argv[i]);
        }
        
            
        // help option specified
        if (this->show_help_) {
            this->write_help(std::cout, filesys::get_path_leaf(argv[0]).c_str());
            exit(0);
        }
        
        // show version
        if (this->show_version_) {
            std::cout << this->version_ << std::endl;
            exit(0);
        }          
    }
}

OptionArg* OptionParser::get_option_ptr(const char * flag) {
    std::map< std::string, OptionArg * >::iterator oai = this->key_opt_map_.find(flag);
    assert (oai != this->key_opt_map_.end() );
    return oai->second;
}

OptionArg& OptionParser::get_option(const char * flag) {
    return *(get_option_ptr(flag));
}

bool OptionParser::is_set(const char * flag) {
    OptionArg& oa = this->get_option(flag);
    return oa.is_set();
}

} // namespace gingko
