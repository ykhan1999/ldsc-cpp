
# ldsc-cpp

> **Acknowledgement:** This project is an independent C++ reimplementation inspired by the original **LDSC (LD Score Regression)** developed by Brendan K. Bulik-Sullivan, Hilary K. Finucane, Po-Ru Loh, Alkes Price, and collaborators. The original LDSC is licensed under GPL-3.0 and is available at the Price Lab’s repository. This project follows GPL-3.0 and owes substantial conceptual credit to those authors.

What is this?
-------------
ldsc-cpp provides two subcommands:
- ldsc munge  — QC/convert GWAS summary statistics to a compact LDSC-style format.
- ldsc ldsc    — options for either rg or partitioned heritability via two-step weighted regression on LD scores

It outputs human-readable logs and a *.summary.txt with headline estimates.

Getting started: Prebuilt binaries
-------------------------------
Linux x86_64: https://github.com/ykhan1999/ldsc-cpp/releases/download/v0.6.0/ldsc_linux-x86_64_glibc2.17.tar.gz

MacOS (all architectures): https://github.com/ykhan1999/ldsc-cpp/releases/download/v0.6.0/ldsc_macos-universal.tar.gz

Installation Instructions:

Linux:
```bash
wget https://github.com/ykhan1999/ldsc-cpp/releases/download/v0.6.0/ldsc_linux-x86_64_glibc2.17.tar.gz
tar xzf ldsc_linux-x86_64_glibc2.17.tar.gz
chmod +x ldsc
./ldsc --version
# optional sanity
ldd ./ldsc
```

MacOS:
```bash
wget https://github.com/ykhan1999/ldsc-cpp/releases/download/v0.6.0/ldsc_macos-universal.tar.gz
tar xzf ldsc_macos-universal.tar.gz
chmod +x ldsc
xattr -d com.apple.quarantine ldsc
./ldsc --version
```
Usage
-------------
Please see https://github.com/ykhan1999/ldsc-cpp/wiki/Quick-start

Building from Source
-----
Prerequisites: CMake >= 3.20, a C++17 compiler, and (recommended) Ninja.

Configure & build (RelWithDebInfo):
  cmake --preset default
  cmake --build --preset default --parallel

Binary path:
  build/default/ldsc    (or ldsc.exe on Windows)

OpenMP:
  - Enabled automatically if found.
  - macOS: brew install libomp
  - Windows: MSVC works out-of-the-box; MinGW needs -fopenmp (handled by CMake).

Release + package:
  cmake --preset release
  cmake --build --preset release --parallel
  cpack --preset release       # creates .zip/.tar.gz in build/release

Notes & differences vs. the original LDSC implementation
--------------------------------------------------------
- Accepts per-chr directories OR {} patterns (e.g., ld/{}.l2.ldscore).
- If per-chr *.l2.M_5_50 files exist, their totals are summed and broadcast as M; otherwise M defaults to SNP count.

License
-------
This project is licensed GPL-3.0-only. See LICENSE.
A NOTICE file documents inspiration and upstream attribution.

Citation
--------
If you use this tool, please cite the original LD Score Regression papers by the LDSC authors in addition to this repository.

Final acknowledgement
---------------------
Deep credit and thanks to the creators and maintainers of the original LDSC. This reimplementation stands
on their shoulders; any mistakes are ours.
