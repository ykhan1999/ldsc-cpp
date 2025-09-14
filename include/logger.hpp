// ldsc-cpp — fast LD Score Regression tools (CLI)
// Copyright (C) 2025  Yousef Khan
// Project: https://github.com/ykhan1999/ldsc-cpp
// Contact (electronic): yousefkhan125@gmail.com
// Contact (paper mail): 6304 Holland Meadow Ln, Gaithersburg, MD 20882
//
// Acknowledgment: This project is an independent C++ reimplementation inspired by
// the original LDSC (ldsc) work by Bulik-Sullivan et al. (Price Lab) and related
// contributors. All credit for the original method and reference implementation
// belongs to the LDSC authors. See the CREDITS and NOTICE files for details.
//
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.


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

