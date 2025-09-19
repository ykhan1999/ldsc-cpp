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

#include "subcommand.hpp"
#include "munge.hpp"
#include "ph2.hpp"
#include "ldsc.hpp"
#include <iostream>
#include <string>

// -------- Safe defaults for metadata (avoid #ifdef in stream expressions) ----
#ifndef PROJECT_NAME_STR
#  define PROJECT_NAME_STR "ldsc-cpp"
#endif
#ifndef PROJECT_VERSION_STR
#  define PROJECT_VERSION_STR "0.4.0"
#endif
#ifndef PROJECT_LICENSE_STR
#  define PROJECT_LICENSE_STR "GPL-3.0-or-later"
#endif
#ifndef PROJECT_URL_STR
#  define PROJECT_URL_STR "https://github.com/ykhan1999/ldsc-cpp"
#endif
#ifndef PROJECT_AUTHOR_STR
#  define PROJECT_AUTHOR_STR "Yousef Khan"
#endif

// -------- GNU license texts ---------------------------------------------------
static const char* kWarranty =
"This program comes with ABSOLUTELY NO WARRANTY.\n"
"For details see the \"Disclaimer of Warranty\" section of the GNU GPL v3:\n"
"https://www.gnu.org/licenses/gpl-3.0.html#warranty\n";

static const char* kConditions =
"This is free software; you can redistribute it and/or modify it under\n"
"the terms of the GNU General Public License as published by the Free\n"
"Software Foundation; version 3 or (at your option) any later version.\n"
"Full text: https://www.gnu.org/licenses/gpl-3.0.html\n";

static const char* kCredits =
"Credits / Acknowledgments\n"
"  ldsc-cpp is an independent C++ reimplementation inspired by the\n"
"  original LDSC (ldsc) work by Bulik-Sullivan et al. (Price Lab)\n"
"  and collaborators. All credit for the original methodology and\n"
"  reference implementation belongs to the LDSC authors.\n";

// -------- UX helpers ----------------------------------------------------------
static int usage() {
  std::cerr <<
    "Usage:\n"
    "  ldsc munge [munge-args...]\n"
    "  ldsc ph2   [partitioned-h2-args...]\n"
    "  ldsc ldsc  [ldsc-args...]\n"
    "\n"
    "Global options:\n"
    "  --version     Show version (with license and LDSC credits)\n"
    "  --license     Show GPL conditions\n"
    "  --credits     Show method/project acknowledgments (LDSC)\n"
    "  --warranty    Show warranty disclaimer\n";
  return 2;
}

static void interactive_banner() {
  std::cout
    << PROJECT_NAME_STR << " " << PROJECT_VERSION_STR << " — " << PROJECT_LICENSE_STR << "\n"
    << PROJECT_URL_STR << "\n"
    << "Acknowledgment: inspired by the original LDSC (ldsc) work by Bulik-Sullivan\n"
       "et al. (Price Lab) and collaborators.\n"
       "This program comes with ABSOLUTELY NO WARRANTY; for details type `show w`.\n"
       "This is free software; you may redistribute it under certain conditions;\n"
       "type `show c` for details. Type `show credits` for acknowledgments.\n\n"
       "Type `quit` to exit.\n";

  std::string line;
  while (true) {
    std::cout << "> ";
    if (!std::getline(std::cin, line)) break;
    if (line == "show w")            std::cout << kWarranty;
    else if (line == "show c")       std::cout << kConditions;
    else if (line == "show credits") std::cout << kCredits;
    else if (line == "quit" || line == "exit") break;
    else if (!line.empty())          std::cout << "Unknown command. Try `show w`, `show c`, `show credits`, or `quit`.\n";
  }
}

int main(int argc, char** argv) {
  // No args → interactive banner (GNU guidance for interactive mode)
  if (argc < 2) {
    interactive_banner();
    return 0;
  }

  std::string cmd = argv[1];

  // Global flags that short-circuit execution
  if (cmd == "--version" || cmd == "-V") {
    std::cout
      << PROJECT_NAME_STR << " " << PROJECT_VERSION_STR << " — " << PROJECT_LICENSE_STR << "\n"
      << PROJECT_URL_STR << "\n"
      << "Acknowledgment: inspired by the original LDSC (ldsc) work by Bulik-Sullivan\n"
         "et al. (Price Lab) and collaborators.\n";
    return 0;
  }
  if (cmd == "--license")  { std::cout << kConditions; return 0; }
  if (cmd == "--credits")  { std::cout << kCredits;    return 0; }
  if (cmd == "--warranty") { std::cout << kWarranty;   return 0; }

  // Subcommands
  if (cmd == "munge") {
    auto sub = shift_argv(argc, argv, 2);
    return run_munge(static_cast<int>(sub.size()), sub.data());
  } else if (cmd == "ph2") {
    auto sub = shift_argv(argc, argv, 2);
    return run_ph2(static_cast<int>(sub.size()), sub.data());
  } else if (cmd == "ldsc") {
    auto sub = shift_argv(argc, argv, 2);
    return run_ldsc(static_cast<int>(sub.size()),sub.data());
  } else {
    std::cerr << "Unknown subcommand or option: " << cmd << "\n";
    return usage();
  }
}

