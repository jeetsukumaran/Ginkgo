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

#include <cassert>
#include <iostream>
#include <ctime>
#include <sstream>
#include <iomanip>
#include "logging.hpp"

using namespace ginkgo;

///////////////////////////////////////////////////////////////////////////////
// LogStreamManager

// singleton instance
LogStreamManager LogStreamManager::instance_;

// free allocate stream pointerss
LogStreamManager::~LogStreamManager() {
    for (std::vector<std::ofstream *>::iterator ofi = this->ofstream_ptrs_.begin();
            ofi != this->ofstream_ptrs_.end();
            ++ofi) {
        assert(*ofi != NULL);
        if ((*ofi)->is_open()) {
            (*ofi)->close();
        }
        delete (*ofi);
    }
}

// creates (opens) a new stream
std::ofstream& LogStreamManager::create(const std::string& path) {
    return this->create(path.c_str());
}

// creates (opens) a new stream
std::ofstream& LogStreamManager::create(const char * path) {
    std::ofstream * out = new std::ofstream(path);
    if (!*out) {
        throw LoggerIOError(path);
    }
    this->ofstream_ptrs_.push_back(out);
    return *(this->ofstream_ptrs_.back());
}

// LogStreamManager
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// LogHandler

// constructor with stream
LogHandler::LogHandler(std::ostream& dest, LogHandler::LoggingLevel logging_level)
        : dest_(&dest),
          logging_level_(logging_level) {
}

// copy constructor
LogHandler::LogHandler(const LogHandler& handler) {
    *this = handler;
}

// assignment
const LogHandler& LogHandler::operator=(const LogHandler& handler) {
    this->path_ = handler.path_;
    this->dest_ = handler.dest_;
    this->logging_level_ = handler.logging_level_;
    return *this;
}

// destructor
LogHandler::~LogHandler() {}

// returns the logging level
LogHandler::LoggingLevel LogHandler::get_logging_level() {
    return this->logging_level_;
}

// sets the logging level
void  LogHandler::set_logging_level(LogHandler::LoggingLevel logging_level) {
    this->logging_level_ = logging_level;
}

// writes to log
void LogHandler::emit(LoggingLevel level, const char * message) {
    if (level >= this->logging_level_) {
        assert(this->dest_);
        (*(this->dest_)) << message << std::endl;
    }
}

// writes to log
void LogHandler::emit(LoggingLevel level, const std::string& message) {
    this->emit(level, message.c_str());
}

// LogHandler
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Logger

// Constructor: gets reference to LogStreamManager singleton object.
Logger::Logger():
            log_stream_manager_(LogStreamManager::get_instance()),
            elapsed_time_leader_("T+"),
            is_show_elapsed_time_(true) {
    time (&this->launch_time_);
    this->reset_elapsed_time();
}

// Destructor: no-op
Logger::~Logger() {

}

// adds a new handler
LogHandler& Logger::add_handler(const std::string& path, LogHandler::LoggingLevel logging_level) {
    return this->add_handler(path.c_str(), logging_level);
}

// adds a new handler
LogHandler& Logger::add_handler(const char * path, LogHandler::LoggingLevel logging_level) {
    std::ostream& out = this->log_stream_manager_.create(path);
    return this->add_handler(out, logging_level);
}

// adds a new handler
LogHandler& Logger::add_handler(std::ostream& dest, LogHandler::LoggingLevel logging_level) {
    this->log_handlers_.push_back(LogHandler(dest, logging_level));
    return this->log_handlers_.back();
}

// resets the log start time
void Logger::reset_elapsed_time() {
    time (&this->registered_start_time_);
}

// writes a debug message
void Logger::debug(const std::string& message) {
    this->log(LogHandler::LOG_DEBUG, message.c_str());
}

// writes a debug message
void Logger::debug(const char * message) {
    this->log(LogHandler::LOG_DEBUG, message);
}

// writes a info message
void Logger::info(const std::string& message) {
    this->log(LogHandler::LOG_INFO, message.c_str());
}

// writes a info message
void Logger::info(const char * message) {
    this->log(LogHandler::LOG_INFO, message);
}

// writes a warning message
void Logger::warning(const std::string& message) {
    this->log(LogHandler::LOG_WARNING, message.c_str());
}

// writes a warning message
void Logger::warning(const char * message) {
    this->log(LogHandler::LOG_WARNING, message);
}

// writes a error message
void Logger::error(const std::string& message) {
    this->log(LogHandler::LOG_ERROR, message.c_str());
}

// writes a error message
void Logger::error(const char * message) {
    this->log(LogHandler::LOG_ERROR, message);
}

// writes a critical message
void Logger::critical(const std::string& message) {
    this->log(LogHandler::LOG_CRITICAL, message.c_str());
}

// writes a critical message
void Logger::critical(const char * message) {
    this->log(LogHandler::LOG_CRITICAL, message);
}

// general log-level message
void Logger::log(LogHandler::LoggingLevel logging_level, const std::string& message) {
    this->log(logging_level, message.c_str());
}

// general log-level message
void Logger::log(LogHandler::LoggingLevel logging_level, const char * message) {
    this->set_timestamps();
    for (std::vector<LogHandler>::iterator logh = this->log_handlers_.begin();
            logh != this->log_handlers_.end();
            ++logh) {
        std::ostringstream msg;
        msg << "[" << this->current_timestamp_;
        if (this->is_show_elapsed_time_) {
            msg << " -- " << this->elapsed_time_leader_ << std::fixed << std::setprecision(4) << this->elapsed_hours_ << "h";
        }
        msg << "]";
        msg << " " << message;
        (*logh).emit(logging_level, msg.str());
    }
}

void Logger::set_timestamps() {
    time_t rawtime;
    time ( &rawtime );
    struct tm * timeinfo = localtime ( &rawtime );
    strftime(this->current_timestamp_, 80, "%Y-%m-%d %H:%M:%S", timeinfo);
    double elapsed_seconds = difftime(rawtime, this->registered_start_time_);
    this->elapsed_hours_ = elapsed_seconds / 3600;
}

double Logger::get_hours_since_launched() {
    time_t rawtime;
    time ( &rawtime );
    double elapsed_seconds = difftime(rawtime, this->launch_time_);
    return elapsed_seconds / 3600;
}

// Logger
///////////////////////////////////////////////////////////////////////////////
