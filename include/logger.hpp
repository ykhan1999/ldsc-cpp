#pragma once
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>

class Logger {
public:
    explicit Logger(const std::string& path, bool append=false);
    template <class... Ts>
    void log(const Ts&... parts) {
        std::ostringstream oss;
        (oss << ... << parts);
        oss << '\n';
        std::lock_guard<std::mutex> lk(mu());
        if (ofs_) ofs_ << oss.str();
        else      std::cerr << oss.str();
    }
private:
    static std::mutex& mu() { static std::mutex m; return m; }
    std::ofstream ofs_;
};

