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
#include <vector>

inline std::vector<char*> shift_argv(int argc, char** argv, int shift) {
    std::vector<char*> out;
    for (int i=0;i<argc-shift;++i) out.push_back(argv[i+shift]);
    return out;
}

