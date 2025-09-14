# ldsc-cpp

Portable C++17 reimplementation of LDSC utilities:

- `ldsc munge` – QC and convert GWAS summary stats into LDSC format.
- `ldsc ph2` – partitioned heritability estimation via weighted regression on LD scores.

## Build

Prerequisites: CMake ≥ 3.20, C++17 compiler, Ninja (recommended).

```bash
cmake --preset default
cmake --build --preset default --parallel

