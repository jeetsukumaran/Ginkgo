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

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <utility>
#include <stdexcept>

namespace gingko {

const unsigned CMDOPTS_LINE_WIDTH = 78;
const unsigned CMDOPTS_OPTION_COL_WIDTH = 24;


class OptionValueTypeError : public std::runtime_error {
    public:
        OptionValueTypeError(const char * msg) : std::runtime_error(msg) {}
};

//! Wraps up a single option, including tracking flags and writing help.
class OptionArg {

    public:
    
        //! For validation and type casting.
        enum option_arg_type {
            STRING,
            INTEGER,
            REAL,
            BOOLEAN
        };    
        
        //! Default constructor.
        OptionArg(const char * help=NULL, const char * meta_var=NULL);      
        
        //! Default destructor.
        virtual ~OptionArg();
        
        //! Wraps a line of text to a given screen width.
        std::string textwrap(const std::string& source, 
                         unsigned line_width=CMDOPTS_LINE_WIDTH,
                         unsigned first_line_indent=0, 
                         unsigned subsequent_line_indent=0) const;                         
        
        //! Composes and writes out help entry for this option.
        std::ostream& write_help(std::ostream& out) const;
        
        std::string get_short_flag() const {
            return this->short_flag_;
        }
        
        void set_short_flag(const std::string s) {
            this->short_flag_ = s;
        }
        
        std::string get_long_flag() const {
            return this->long_flag_;
        }
        
        void set_long_flag(const std::string s) {
            this->long_flag_ = s;
        }
        
        std::string get_help() const {
            return this->help_;
        }
        
        void set_help(const char* s) {
            this->help_ = s;
        }
                
        std::string get_meta_var() const {
            return this->meta_var_;
        }
        
        void set_meta_var(const char*  s) {
            this->meta_var_ = s;
        }
        
        bool is_switch() const {
            return this->is_switch_;
        }
        
        void set_is_switch(bool val) {
            this->is_switch_ = val;
        }
        
        bool& is_set() {
            return this->is_set_;
        }
        
        virtual void process_value_string(const std::string& val_str) = 0;
  
    private:
        std::string     short_flag_;
        std::string     long_flag_;
        std::string     help_;
        std::string     meta_var_;
        bool            is_switch_;        
        bool            is_set_;

};

template <typename T>
class TypedOptionArg : public OptionArg {

    public:
        TypedOptionArg(void * store,
                       const char * help=NULL,
                       const char * meta_var=NULL,
                       void * default_value=NULL)
            : OptionArg(help, meta_var),
              default_value_() {
            if (store != NULL) {
                this->set_store(store);
            }
            if (default_value != NULL) {
                this->set_default_value(default_value);
                *this->store_ = this->default_value_;
            }
        }   
        
        virtual ~TypedOptionArg() {}
        
        T * get_store() {
            return this->store_;
        }
                
        void set_store(T * store) {
            this->store_ = store;
        }       

        void set_store(void * store) {
            this->store_ = static_cast<T *>(store);
        }       
        
        void process_value(const T& val) {
            *this->store_ = val;
        }      
        
        virtual void process_value_string(const std::string& val_str) {
            std::istringstream istr(val_str);
            T temp;
            istr >> temp;
            if (!istr.fail() and istr.eof()) {
                *this->store_ = temp;
                this->is_set() = true;
            } else {
                std::string msg;
                msg = "failed to convert \"" + val_str + "\"";
                throw OptionValueTypeError(msg.c_str());
            }
        }

        void set_default_value(void * val) {
            if (val != NULL) {
                this->default_value_ = *(static_cast<T *>(val));
            }                
        }

    private:
        T *     store_;
        T       default_value_;
};

///////////////////////////////////////////////////////////////////////////////
//! General option parser.
class OptionParser {

    public:
        
        OptionParser();            
        ~OptionParser();
        
        //! Add an option to those supported by the program. Must supply at
        //! least one of the following: short_flag or long_flag.
        //! Short flags start with a dash and are followed by one character
        //! (e.g., "-f").
        //! Long flags start with two dashes and are followed by one or more
        //! characters (e.g., "--filename").
        template <typename T>
        OptionArg * add_option(void * store,
                               const char * short_flag=NULL,
                               const char * long_flag=NULL,
                               const char * help=NULL,
                               const char * meta_var=NULL,
                               void * default_value=NULL) {                               
            OptionArg * oa;                              
            oa = new TypedOptionArg<T>(store, help, meta_var, default_value);
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
        
        OptionArg * add_switch(void * store,
                               const char * short_flag=NULL,
                               const char * long_flag=NULL,
                               const char * help=NULL,
                               const char * meta_var=NULL,
                               void * default_value=NULL) {
            OptionArg * switch_arg = this->add_option<bool>(store, short_flag, long_flag, help, meta_var, default_value);
            switch_arg->set_is_switch(true);
            return switch_arg;
        }                        
               
        //! Client must call this, passing in arguments from main().               
        void parse(int argc, char * argv[]);                
        
        //! Checks to see if a particular flag is set (i.e., specified by the 
        //! user.
        bool is_set(const char * flag);                
        
        //! Copies value into given data field.
        void store_value(const char * flag, void * data);

    private:
        std::ostream& write_help(std::ostream& out) const;               
        OptionArg* get_option_ptr(const char * flag);        
        OptionArg& get_option(const char * flag);     

    private:
        bool                                        show_help_;
        OptionArg *                                 help_option_;
        std::vector<OptionArg *>                    option_args_;
        std::vector< std::string >                  pos_args_;
        std::map< std::string, OptionArg * >        key_opt_map_;
};

} // namespace gingko