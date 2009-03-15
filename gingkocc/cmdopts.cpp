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

#include "cmdopts.h"

using namespace gingko;

///////////////////////////////////////////////////////////////////////////////
// OptionArg

OptionArg::OptionArg()
    : is_switch_(false),
      is_set_(false) 
  {}
  
OptionArg::~OptionArg() {}
  
std::string OptionArg::textwrap(const std::string& source, 
        unsigned line_width,
        unsigned first_line_indent, 
        unsigned subsequent_line_indent) const {
    std::ostringstream wrapped;
    unsigned col_count = 1;
    unsigned line_count = 1;
    for (std::string::const_iterator s = source.begin();
            s != source.end();
            ++s, ++col_count) {
        if (*s == '\n' or col_count > line_width) {
            wrapped << '\n';
            col_count = 0;
            line_count += 1;
            continue;
        }        
        if (col_count == 1 and line_count == 1 and first_line_indent > 0) {
            for (unsigned i = 0; i < first_line_indent; ++i) {
                wrapped << ' ';
            }
            col_count += first_line_indent;
        } else if (col_count == 1 and line_count > 1) {
            for (unsigned i = 0; i < subsequent_line_indent; ++i) {
                wrapped << ' ';
            }
            col_count += subsequent_line_indent;                    
        }
        wrapped << *s;
    }
    return wrapped.str();                         
} 

std::ostream& OptionArg::write_help(std::ostream& out) const {            
    std::string help_str;        
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
        if (help_str.size() > CMDOPTS_OPTION_COL_WIDTH) {
            help_str += "\n";
        } else {
            while (help_str.size() < CMDOPTS_OPTION_COL_WIDTH) {
                help_str += " ";
            }                
        }
        help_str += this->help_;
        std::string help_desc = this->textwrap(help_str, 
                                               CMDOPTS_LINE_WIDTH, 
                                               0, 
                                               CMDOPTS_OPTION_COL_WIDTH);
        help_str = help_desc;
    }                            
    out << help_str; 
    return out;
}

///////////////////////////////////////////////////////////////////////////////
// OptionParser

OptionParser::OptionParser() {
    this->help_option_ = dynamic_cast<BooleanOptionArg *>(this->add_option("-h", "--help", OptionArg::BOOLEAN,
                                           "show this message and exit"));
}
    
OptionParser::~OptionParser() {
    for (std::vector<OptionArg *>::iterator oap = this->option_args_.begin();
            oap != this->option_args_.end();
            ++oap) {
        delete *oap;                    
    }                    
}

OptionArg * OptionParser::add_option(const char * short_flag,
                                    const char * long_flag,
                                    OptionArg::option_arg_type val_type,
                                    const char * help,
                                    const char * meta_var) {
    OptionArg * oa;                              
    if (val_type == OptionArg::STRING) {            
        oa = new StringOptionArg();
    } else if (val_type == OptionArg::INTEGER) {      
        oa = new IntegerOptionArg();
    } else if (val_type == OptionArg::REAL) { 
        oa = new RealOptionArg();            
    } else if (val_type == OptionArg::BOOLEAN) { 
        oa = new BooleanOptionArg();            
    }
    assert ( oa );
    assert( short_flag != NULL or long_flag != NULL);
    if (short_flag != NULL) {
        oa->set_short_flag(short_flag);
    }
    if (long_flag != NULL) {
        oa->set_long_flag(long_flag);
    }   
    
    if (help != NULL) {
        oa->set_help(help);
    }   
    if (meta_var != NULL) {
        oa->set_meta_var(meta_var);
    } else if (long_flag != NULL) {
        oa->set_meta_var(long_flag);
    } else if (short_flag != NULL) {
        oa->set_meta_var(short_flag);
    }
    
//     if (val_type == OptionArg::STRING) {            
//         StringOptionArg * str_opt = dynamic_cast<StringOptionArg *>(oa);
//         str_opt->set_value(*static_cast<const std::string *>(default_value));
//     } else if (val_type == OptionArg::INTEGER) {            
//         IntegerOptionArg * int_opt = dynamic_cast<IntegerOptionArg *>(oa);
//         int_opt->set_value(*static_cast<const int *>(default_value));
//     } else if (val_type == OptionArg::REAL) {            
//         RealOptionArg * double_opt = dynamic_cast<RealOptionArg *>(oa);
//         double_opt->set_value(*static_cast<const double *>(default_value));
//     } else if (val_type == OptionArg::BOOLEAN) {
//         BooleanOptionArg * bool_opt = dynamic_cast<BooleanOptionArg *>(oa);
//         bool_opt->set_value(*static_cast<const bool *>(default_value));
//     }
    
    this->option_args_.push_back(oa);
    if (short_flag) {
        assert(short_flag[0] == '-' and short_flag[1] != 0 and short_flag[1] != '-');
        assert(this->key_opt_map_.find(short_flag) == this->key_opt_map_.end());
        this->key_opt_map_.insert(std::make_pair(short_flag, oa));
    }        
    if (long_flag) {
        assert(long_flag[0] == '-' and long_flag[1] =='-' and long_flag[2] != '-');
        assert(this->key_opt_map_.find(long_flag) == this->key_opt_map_.end());
        this->key_opt_map_.insert(std::make_pair(long_flag, oa));
    }        
    return oa;
}

std::ostream& OptionParser::write_help(std::ostream& out) const {
    for (std::vector<OptionArg *>::const_iterator oa = this->option_args_.begin();
            oa != this->option_args_.end();
            ++oa) {
        (*oa)->write_help(out);
        out << std::endl;
    }
    return out;
}

void OptionParser::parse(int argc, char * argv[]) {

    for (int i = 0; i < argc; ++i) { 
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
                std::cerr << "unrecognized command \"" << arg_name << "\"" << std::endl;
                exit(1);
            }
            
            if (oai->second == this->help_option_) {
                this->write_help(std::cerr);
                exit(1);
            }
            
            OptionArg& oa = *(oai->second);
            
            if (not oa.get_val_type() == OptionArg::BOOLEAN) {
                if (arg_value.size() == 0) {
                    if (i == argc-1) {
                        std::cerr << "expecting value for option \"" << arg_name << "\"" << std::endl;
                        exit(1);
                    } else {
                        arg_value = argv[i+1];
                        i += 1;
                    }                            
                }
                if (oa.get_val_type() == OptionArg::STRING) {            
                    StringOptionArg * str_opt = dynamic_cast<StringOptionArg *>(&oa);
                    str_opt->set_value(arg_value);
                } else if (oa.get_val_type() == OptionArg::INTEGER) {            
                    IntegerOptionArg * int_opt = dynamic_cast<IntegerOptionArg *>(&oa);
                    int_opt->set_value(atol(arg_value.c_str()));
                } else if (oa.get_val_type() == OptionArg::REAL) {            
                    RealOptionArg * double_opt = dynamic_cast<RealOptionArg *>(&oa);
                    double_opt->set_value(atof(arg_value.c_str()));
                }                        
            } else {
                BooleanOptionArg * bool_opt = static_cast<BooleanOptionArg *>(&oa);
                bool_opt->set_value(true);            
            }
            
        } else {
            this->pos_args_.push_back(argv[i]);
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

void OptionParser::store_value(const char * flag, void * data) {
    OptionArg * oa = this->get_option_ptr(flag);
    if (oa->get_val_type() == OptionArg::STRING) {            
        StringOptionArg * str_opt = dynamic_cast<StringOptionArg *>(this->get_option_ptr(flag));
        *(static_cast<std::string *>(data)) = str_opt->get_value();
    } else if (oa->get_val_type() == OptionArg::INTEGER) {            
        IntegerOptionArg * int_opt = dynamic_cast<IntegerOptionArg *>(this->get_option_ptr(flag));
        *(static_cast<int *>(data)) = int_opt->get_value();
    } else if (oa->get_val_type() == OptionArg::REAL) {            
        RealOptionArg * double_opt = dynamic_cast<RealOptionArg *>(this->get_option_ptr(flag));
        *(static_cast<double *>(data)) = double_opt->get_value();
    } else if (oa->get_val_type() == OptionArg::BOOLEAN) {            
        BooleanOptionArg * bool_opt = dynamic_cast<BooleanOptionArg *>(this->get_option_ptr(flag));
        *(static_cast<bool *>(data)) = bool_opt->get_value();
    } else {
        assert(0);
    }
}
