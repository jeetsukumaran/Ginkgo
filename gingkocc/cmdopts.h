
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <utility>

namespace gingko {

const unsigned CMDOPTS_LINE_WIDTH = 78;
const unsigned CMDOPTS_OPTION_COL_WIDTH = 24;

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
        OptionArg();      
        
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
        
        option_arg_type get_val_type() const {
            return this->val_type_;
        }
        
        void set_val_type(option_arg_type val) {
            this->val_type_ = val;
        }
                
    private:
        std::string     short_flag_;
        std::string     long_flag_;
        std::string     help_;
        std::string     meta_var_;
        bool            is_switch_;        
        bool            is_set_;
        option_arg_type val_type_;

};

class StringOptionArg : public OptionArg {

    public:
        StringOptionArg() {
            this->set_val_type(OptionArg::STRING);
            this->set_is_switch(false);
        }        
        virtual ~StringOptionArg() {}    
        
        std::string get_default_value() const {
            return this->default_value_;
        }
       
        void set_default_value(const std::string& value) {
            this->default_value_ = value;
        }
        
        std::string get_value() const {
            return this->value_;
        }
       
        void set_value(const std::string& value) {
            this->value_ = value;
        }         
               
    private:
        std::string     default_value_;
        std::string     value_;
};

//! Specialization for integer option.
class IntegerOptionArg : public OptionArg {
    public:
        IntegerOptionArg() {
            this->set_val_type(OptionArg::INTEGER);
            this->set_is_switch(false);
        }           
        virtual ~IntegerOptionArg() {}    
        
        long get_default_value() const {
            return this->default_value_;
        }
       
        void set_default_value(const long& value) {
            this->default_value_ = value;
        }
        
        long get_value() const {
            return this->value_;
        }
       
        void set_value(const long& value) {
            this->value_ = value;
        }         
               
    private:
        long     default_value_;
        long     value_;
};

//! Specialization for floating point option.
class RealOptionArg : public OptionArg {
    public:
        RealOptionArg() {
            this->set_val_type(OptionArg::REAL);
            this->set_is_switch(false);
        }           
        virtual ~RealOptionArg() {}    
        
        double get_default_value() const {
            return this->default_value_;
        }
       
        void set_default_value(const double& value) {
            this->default_value_ = value;
        }
        
        double get_value() const {
            return this->value_;
        }
       
        void set_value(const double& value) {
            this->value_ = value;
        }         
               
    private:
        double     default_value_;
        double     value_;
};

//! Specialization for boolean option.
class BooleanOptionArg : public OptionArg {
    public:
        BooleanOptionArg() {
            this->set_val_type(OptionArg::BOOLEAN);
            this->set_is_switch(true);
        }           
        virtual ~BooleanOptionArg() {}    
        
        bool get_default_value() const {
            return this->default_value_;
        }
       
        void set_default_value(const bool& value) {
            this->default_value_ = value;
        }
        
        bool get_value() const {
            return this->value_;
        }
       
        void set_value(const bool& value) {
            this->value_ = value;
        }         
               
    private:
        bool     default_value_;
        bool     value_;
};

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
        OptionArg * add_option(const char * short_flag=NULL,
                               const char * long_flag=NULL,
                               OptionArg::option_arg_type val_type=OptionArg::STRING,
                               const char * help=NULL,
                               const char * meta_var=NULL);        
               
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
        BooleanOptionArg *                          help_option_;
        std::vector<OptionArg *>                    option_args_;
        std::vector< std::string >                  pos_args_;
        std::map< std::string, OptionArg * >        key_opt_map_;
};

} // namespace gingko