#include "subcommand.hpp"
#include "munge.hpp"
#include "ph2.hpp"
#include <iostream>
#include <string>

static int usage() {
    std::cerr <<
      "Usage:\n"
      "  ldsc munge [munge-args...]\n"
      "  ldsc ph2   [partitioned-h2-args...]\n";
    return 2;
}

int main(int argc, char** argv) {
    if (argc < 2) return usage();
    std::string cmd = argv[1];

    if (cmd == "munge") {
        auto sub = shift_argv(argc, argv, 2);
        return run_munge(static_cast<int>(sub.size()), sub.data());
    } else if (cmd == "ph2") {
        auto sub = shift_argv(argc, argv, 2);
        return run_ph2(static_cast<int>(sub.size()), sub.data());
    } else {
        std::cerr << "Unknown subcommand: " << cmd << "\n";
        return usage();
    }
}

