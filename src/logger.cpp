#include "logger.hpp"

Logger::Logger(const std::string& path, bool append)
    : ofs_(path, append ? std::ios::app : std::ios::out)
{
    if (!ofs_) {
        std::lock_guard<std::mutex> lk(mu());
        std::cerr << "WARN: could not open log file '" << path
                  << "'; logging to stderr\n";
    }
}

