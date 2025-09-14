#pragma once
#include <vector>

inline std::vector<char*> shift_argv(int argc, char** argv, int shift) {
    std::vector<char*> out;
    for (int i=0;i<argc-shift;++i) out.push_back(argv[i+shift]);
    return out;
}

