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

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <ctime>

#if !defined(LOGGING_H)
#define LOGGING_H

namespace ginkgo {

///////////////////////////////////////////////////////////////////////////////
// Exceptions

class LoggerIOError : public std::runtime_error {
    public:
        LoggerIOError(const char * msg) : std::runtime_error(msg) {}
        LoggerIOError(const std::string& msg) : std::runtime_error(msg.c_str()) {}
};

// Exceptions
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// LogStreamManager

class LogStreamManager {

    public:
        ~LogStreamManager();
        std::ofstream& create(const std::string& path);
        std::ofstream& create(const char * path);

    private:
        std::vector<std::ofstream *>    ofstream_ptrs_;

    ///////////////////////////////////////////////////////////////////////////
    // Singleton infrastructure

    public:
        static LogStreamManager& get_instance() {
            return LogStreamManager::instance_;
        }

    private:
        static LogStreamManager instance_;
        LogStreamManager() {}
        LogStreamManager(const LogStreamManager &);
        LogStreamManager & operator=(const LogStreamManager &);

};

// LogStreamManager
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// LogHandler

class LogHandler {

    public:

        enum LoggingLevel {
            LOG_DEBUG = 0,
            LOG_INFO = 1,
            LOG_WARNING = 2,
            LOG_ERROR = 3,
            LOG_CRITICAL = 4
        };

        /**
         * Constructs a log handler that writes all messages of "log_level"
         * or higher to output stream given by "dest".
         * @param dest      destintation stream to which to write log messages
         * @param log_level minimum level of message to handle
         */
        LogHandler(std::ostream& out, LoggingLevel logging_level=LogHandler::LOG_DEBUG);

        /**
         * Copy constructor.
         */
        LogHandler(const LogHandler& handler);

        /**
         * Assigment.
         */
        const LogHandler& operator=(const LogHandler& handler);

        /**
         * Destructor.
         */
        ~LogHandler();

        /**
         * Returns the minimum level of messages reported by this handler.
         * @returns logging level
         */
        LoggingLevel get_logging_level();

        /**
         * Sets the minimum level of messages reported by this handler.
         * @param @logging_level logging level
         */
        void set_logging_level(LoggingLevel logging_level);

        /**
         * Writes to destination if level is high enough.
         * @param message   message to write (const char *)
         */
        void emit(LoggingLevel level, const char * message);

        /**
         * Writes to destination if level is high enough.
         * @param message   message to write (std::string)
         */
        void emit(LoggingLevel level, const std::string& message);

    private:
        std::ostream*   dest_;
        std::string     path_;
        LoggingLevel    logging_level_;

};

// LogHandler
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Logger

class Logger {

    public:
        /**
         * Constructor: gets reference to LogStreamManager singleton object.
         */
        Logger();

        /**
         * No-op destructor.
         */
        ~Logger();

        /**
         * Creates and adds a new handler.
         * @param path          destination file path to which to write logging messages
         * @param logging_level minimum level of message that this handler will log
         * @returns             reference to logging handler just added.
         */
        LogHandler& add_handler(const std::string& path, LogHandler::LoggingLevel logging_level=LogHandler::LOG_DEBUG);

        /**
         * Creates and adds a new handler.
         * @param path          destination file path to which to write logging messages
         * @param logging_level minimum level of message that this handler will log
         * @returns             reference to logging handler just added.
         */
        LogHandler& add_handler(const char * path, LogHandler::LoggingLevel logging_level=LogHandler::LOG_DEBUG);

        /**
         * Creates and adds a new handler.
         * @param dest          output stream to which to write logging messages
         * @param logging_level minimum level of message that this handler will log
         * @returns             reference to logging handler just added.
         */
        LogHandler& add_handler(std::ostream& dest, LogHandler::LoggingLevel logging_level=LogHandler::LOG_DEBUG);

        /**
         * Sets the log time to NOW.
         */
        void reset_elapsed_time();

        /**
         * Sets the log time to NOW, and simultaneous sets the string
         * leader to be printed in the time definition.
         */
        void reset_elapsed_time(const std::string& leader) {
            this->reset_elapsed_time(leader.c_str());
        }

        /**
         * Sets the log time to NOW, and simultaneous sets the string
         * leader to be printed in the time definition.
         */
        void reset_elapsed_time(const char * leader) {
            this->set_elapsed_time_leader(leader);
            this->reset_elapsed_time();
        }

        /**
         * Sets the token preceding the elapsed time in the timestamp.
         */
        void set_elapsed_time_leader(const std::string& leader) {
            this->elapsed_time_leader_ = leader;
        }

        /**
         * Sets the token preceding the elapsed time in the timestamp.
         */
        void set_elapsed_time_leader(const char * leader) {
            this->elapsed_time_leader_ = leader;
        }

        /**
         * Writes a message at the "DEBUG" level.
         *
         * @param message      message to write
         */
        void debug(const std::string& message);

        /**
         * Writes a message at the "DEBUG" level.
         *
         * @param message      message to write
         */
        void debug(const char * message);

        /**
         * Writes a message at the "INFO" level.
         *
         * @param message      message to write
         */
        void info(const std::string& message);

        /**
         * Writes a message at the "INFO" level.
         *
         * @param message      message to write
         */
        void info(const char * message);

        /**
         * Writes a message at the "WARNING" level.
         *
         * @param message      message to write
         */
        void warning(const std::string& message);

        /**
         * Writes a message at the "WARNING" level.
         *
         * @param message      message to write
         */
        void warning(const char * message);

        /**
         * Writes a message at the "ERROR" level.
         *
         * @param message      message to write
         */
        void error(const std::string& message);

        /**
         * Writes a message at the "ERROR" level.
         *
         * @param message      message to write
         */
        void error(const char * message);

        /**
         * Writes a message at the "CRITICAL" level.
         *
         * @param message      message to write
         */
        void critical(const std::string& message);

        /**
         * Writes a message at the "CRITICAL" level.
         *
         * @param message      message to write
         */
        void critical(const char * message);

        /**
         * Writes a message at the specified level.
         *
         * @param logging_level message severity level
         * @param message       message to write
         */
        void log(LogHandler::LoggingLevel logging_level, const std::string& message);

        /**
         * Writes a message at the specified level.
         *
         * @param logging_level message severity level
         * @param message       message to write
         */
        void log(LogHandler::LoggingLevel logging_level, const char * message);

        /**
         * Composes and returns pretty time stamp.
         */
        void set_timestamps();

        /**
         * Displays elapsed time in log.
         */
        void show_elapsed_time() {
            this->is_show_elapsed_time_ = true;
        }

        /**
         * Hides elapsed time in log.
         */
        void hide_elapsed_time() {
            this->is_show_elapsed_time_ = false;
        }

        /**
         * Returns timestamp.
         */
        const char * get_current_timestamp() {
            this->set_timestamps();
            return this->current_timestamp_;
        }

        /**
         * Return total logging time hours.
         */
        double get_hours_since_launched();

    private:
        LogStreamManager&               log_stream_manager_;
        std::vector<LogHandler>         log_handlers_;
        time_t                          launch_time_;
        time_t                          registered_start_time_;
        char                            current_timestamp_[80];
        double                          elapsed_hours_;
        std::string                     elapsed_time_leader_;
        bool                            is_show_elapsed_time_;

}; // logger

// Logger
///////////////////////////////////////////////////////////////////////////////

}

#endif
