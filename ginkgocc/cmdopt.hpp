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

#include <cassert>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <utility>
#include <stdexcept>
#include <iomanip>
#include "textutil.hpp"


#if !defined(GINGKO_CMDOPTS_H)
#define GINGKO_CMDOPTS_H

namespace gingko {

const unsigned CMDOPTS_LINE_WIDTH = 78;
const unsigned CMDOPTS_OPTION_COL_WIDTH = 24;

/**
 * Thrown when types cannot be converted.
 */
class OptionValueTypeError : public std::runtime_error {
    public:
        OptionValueTypeError(const char * msg) : std::runtime_error(msg) {}
};

/**
 * Base class for an optional argument on the command line.
 */
class OptionArg {

    public:
    
        /**
         * Constructs argument handler with basic help info.
         *
         * @param help          help message
         * @param meta_var      string to display as value in help
         */
        OptionArg(const char * help=NULL, const char * meta_var=NULL);      
        
        /** Destructor. */
        virtual ~OptionArg();
        
        /**
         * Writes help message for this option out to given stream.
         * @param   out     stream to write message to
         * @return          output stream
         */
        std::ostream& write_help(std::ostream& out) const;

        /**
         * Sets the short flag bound to this argument.
         * @param   flag     short flag, including leading dash
         */
        void set_short_flag(const std::string flag) {
            this->short_flag_ = flag;
        }
        
        /**
         * Sets the long flag bound to this argument.
         * @param   flag     long flag, including leading dashes
         */
        void set_long_flag(const std::string flag) {
            this->long_flag_ = flag;
        }
        
        void set_meta_var(const char*  s) {
            this->meta_var_ = textutil::upper(s);
        }
        
        /**
         * Returns <code>true</code> if this option is a switch (controls 
         * a true/false option as indicated by its presence, and hence takes
         * no arguments).
         * @return      <code>true</code> if no arguments needed
         */        
        bool is_switch() const {
            return this->is_switch_;
        }
        
        /**
         * Set to <code>true</code> if this option is a switch (controls 
         * a true/false option as indicated by its presence, and hence takes
         * no arguments).
         * @param val      <code>true</code> if no arguments needed
         */           
        void set_is_switch(bool val) {
            this->is_switch_ = val;
        }
        
        /**
         * Returns <code>true</code> if this option was set by user.
         * @return <code>true</code> if this option was set by user=
         */            
        bool is_set() {
            return this->is_set_;
        }
        
        /**
         * Set to <code>true</code> if this option was set by user.
         * @param set <code>true</code> if this option was set by user
         */            
        void set_is_set(bool set) {
            this->is_set_ = set;
        }        
        
        /**
         * Processes string given as argument to this option.
         *
         * @param val_str       value passed to this option by user (as string)
         */
        virtual void process_value_string(const std::string& val_str) = 0;
        
        /**
         * Returns current value of argument as string.
         * @return current value of argument as string
         */
        virtual std::string current_value_as_string() const = 0;
        
  
    private:
    
        /** short flag string (single dash followed by single character) */
        std::string     short_flag_;
        /** long flag string (two dashes followed by one or more characters) */
        std::string     long_flag_;
        /** help message */
        std::string     help_;
        /** string to display as placeholder for value */
        std::string     meta_var_;
        /** <code>true</code> if this option takes no arguments */        
        bool            is_switch_;        
        /** <code>true</code> if value for this option was set by user */
        bool            is_set_;

};

/**
 * Templated derivation of OptionArg to handle various types of arguments.
 */
template <typename T>
class TypedOptionArg : public OptionArg {

    public:
    
        /**
         * Constructs argument handler with pointer to data store, flag and 
         * help info.
         *
         * @param  store        pointer to data store that will hold argument
         * @param  short_flag   short flag, including leading dash
         * @param  long_flag    long flag, including leading dashes    
         * @param  help         help message
         * @param  meta_var     string to display as value in help
         */
        TypedOptionArg(void * store,
                       const char * short_flag,
                       const char * long_flag,
                       const char * help=NULL,
                       const char * meta_var=NULL)
                : OptionArg(help, meta_var) {
            if (store != NULL) {
                this->set_store(store);
            }
            assert( short_flag != NULL or long_flag != NULL);
            if (short_flag != NULL) {
                this->set_short_flag(short_flag);
            }
            if (long_flag != NULL) {
                this->set_long_flag(long_flag);
            }             
        }   
        
        /** Destructor. */
        virtual ~TypedOptionArg() {}
        
        /**
         * Returns pointer to data store.
         * @return  pointer to data store
         */
        T * get_store() {
            return this->store_;
        }
             
        /**
         * Sets pointer to data store (explicit typing).
         * @param store  pointer to data store
         */             
        void set_store(T * store) {
            this->store_ = store;
        }       

        /**
         * Sets pointer to data store (void typing).
         * @param store  pointer to data store
         */  
        void set_store(void * store) {
            this->store_ = static_cast<T *>(store);
        }       
        
        /**
         * Handles typed value passed to this option.
         * @param val   value for this option
         */
        void process_value(const T& val) {
            *this->store_ = val;
        }      
        
        /**
         * Handles value passed to this argument as tring
         *
         * @param val_str   string representation of the value for this option
         */        
        virtual void process_value_string(const std::string& val_str) {
            std::istringstream istr(val_str);
            T temp;
            istr >> std::noskipws; 
            istr >> temp;
            if (!istr.fail() and istr.eof()) {
                *this->store_ = temp;
                this->set_is_set(true);
            } else {
                std::string msg;
                msg = "failed to convert \"" + val_str + "\"";
                throw OptionValueTypeError(msg.c_str());
            }
        }
        
        /**
         * Returns current value of argument as string.
         * @return current value of argument as string
         */
        virtual std::string current_value_as_string() const {
            assert(this->store_ != NULL);
            std::ostringstream ostr;
            ostr << *this->store_;
            return ostr.str();
        }
                

    private:
    
        /** pointer to storage for value for this option */
        T *     store_;
};


///////////////////////////////////////////////////////////////////////////////
//! General option parser.
class OptionParser {

    public:
    
        typedef std::vector< std::string >  PosArgs;
        
        /** Constructor. */
        OptionParser(const char * version=NULL, const char * description=NULL, const char * usage=NULL);    
        
        /** Destructor: frees memory associated with OptionArg objects. */
        ~OptionParser();
        
        /**
         * Add an option to those supported by the program. 
         * Must supply at least one of the following: short_flag or long_flag.
         * Short flags start with a dash and are followed by one character
         * (e.g., "-f").
         * Long flags start with two dashes and are followed by one or more
         * characters (e.g., "--filename").
         *
         * @param  store        pointer to data store that will hold argument
         * @param  short_flag   short flag, including leading dash
         * @param  long_flag    long flag, including leading dashes    
         * @param  help         help message
         * @param  meta_var     string to display as value in help
         */
        template <typename T>
        OptionArg * add_option(void * store,
                               const char * short_flag=NULL,
                               const char * long_flag=NULL,
                               const char * help=NULL,
                               const char * meta_var=NULL) {                               
            OptionArg * oa;                              
            oa = new TypedOptionArg<T>(store, short_flag, long_flag, help, meta_var);
            assert ( oa );
   
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
            
            if (meta_var != NULL) {
                oa->set_meta_var(meta_var);
            } else if (long_flag != NULL) {
                oa->set_meta_var(long_flag + 2);
            } else if (short_flag != NULL) {
                oa->set_meta_var(short_flag + 1);
            }
            
            return oa;
        }
        
        /**
         * Add boolean (switch) option to those supported by the program. 
         * Must supply least one of the following: short_flag or long_flag.
         * Short flags start with a dash and are followed by one character
         * (e.g., "-f").
         * Long flags start with two dashes and are followed by one or more
         * characters (e.g., "--filename").
         *
         * @param  store        pointer to data store that will hold argument
         * @param  short_flag   short flag, including leading dash
         * @param  long_flag    long flag, including leading dashes    
         * @param  help         help message
         * @param  meta_var     string to display as value in help
         */        
        OptionArg * add_switch(void * store,
                               const char * short_flag=NULL,
                               const char * long_flag=NULL,
                               const char * help=NULL,
                               const char * meta_var=NULL) {
            OptionArg * switch_arg = this->add_option<bool>(store, short_flag, long_flag, help, meta_var);
            switch_arg->set_is_switch(true);
            return switch_arg;
        }
        
        /**
         * Returns copy of the current program usage string.
         * @return  copy of the current program usage string
         */
         std::string get_usage() const {
            return this->usage_;
         }
         
        /**
         * Sets the current program usage string.
         * @param usage  program usage string
         */
         void set_usage(const char * usage) {
            this->usage_ = usage;
         }     
         
        /**
         * Returns copy of the current program description string.
         * @return  copy of the current program description string
         */
         std::string get_description() const {
            return this->description_;
         }
         
        /**
         * Sets the current program description string.
         * @param description  program description string
         */
         void set_description(const char * description) {
            this->description_ = description;
         }
         
        /**
         * Returns copy of the current program version string.
         * @return  copy of the current program version string
         */
         std::string get_version() const {
            return this->version_;
         }
         
        /**
         * Sets the current program version string.
         * @param version  program version string
         */
         void set_version(const char * version) {
            this->version_ = version;
         }          
               
        /**
         * Client must call this, passing in arguments from main().
         *
         * @param argc      number of argument strings in the command line
         * @param argv      array of argument strings
         */
        void parse(int argc, char * argv[]);                
        
        /**
         * Checks to see if a particular OptionArg is set.
         * @param   flag    flag string for option (including leading dashes)
         * @return          <code>true</code> if option is set
         */
        bool is_set(const char * flag);                
        
        /**
         * Copies value reference by pointer to data to OptionArg mapped to 
         * given flag.
         *
         * @param flag      flag for option to be set
         * @param data      pointer to data to set option
         */
        void store_value(const char * flag, void * data);
        
        /**
         * Returns copy of vector of positional (non-flagged) arguments.
         * @return  vector of positional (non-flagged) arguments
         */
        std::vector< std::string > get_args() {
            return this->pos_args_;
        }
        
        /**
         * Writes out help for all options to given output stream.
         * @param out   output stream to which to write help message
         * @return      same output stream
         */
        std::ostream& write_help(std::ostream& out) const;
        
        /**
         * Writes usage info to the given output stream.
         * @param out   output stream to which to write usage info
         * @return      same output stream
         */        
        std::ostream& write_usage(std::ostream& out) const;
        
        /**
         * Writes description to the given output stream.
         * @param out   output stream to which to write description
         * @return      same output stream
         */        
        std::ostream& write_description(std::ostream& out) const;
        
        /**
         * Writes version to the given output stream.
         * @param out   output stream to which to write version
         * @return      same output stream
         */        
        std::ostream& write_version(std::ostream& out) const;         

    private:
        /**
         * Searches for and returns pointer to OptionArg object correspondng to
         * given flag.
         * @param flag  flag for OptionArg to return
         */
        OptionArg* get_option_ptr(const char * flag);
        /**
         * Searches for and returns reference OptionArg object correspondng to
         * given flag.
         * @param flag  flag for OptionArg to return
         */        
        OptionArg& get_option(const char * flag);     

    private:
        /** stores value of help option arg switch */
        bool                                        show_help_;
        /** stores pointer to help option arg */
        OptionArg *                                 help_option_;
        /** stores value of help option arg switch */
        bool                                        show_version_;
        /** stores pointer to help option arg */
        OptionArg *                                 version_option_;        
        /** usage string */
        std::string                                 usage_;
        /** program description */
        std::string                                 description_;        
        /** program version */
        std::string                                 version_;        
        /** collection of OptionArg objects */
        std::vector<OptionArg *>                    option_args_;
        /** positional arguments on the command line */
        PosArgs                                     pos_args_;
        /** map of flag strings to option arguments */
        std::map< std::string, OptionArg * >        key_opt_map_;
        /** the program file name */
        std::string                                 prog_filename_;
};

} // namespace gingko

#endif
