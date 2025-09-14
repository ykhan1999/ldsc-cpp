ldsc-cpp

Acknowledgement (credit):
This project is an independent C++ reimplementation inspired by the original LDSC (LD Score Regression)
developed by Brendan K. Bulik-Sullivan, Hilary K. Finucane, Po-Ru Loh, Alkes Price, and collaborators.
The original LDSC is licensed under GPL-3.0 and is available at the Price Lab’s repository. This project
follows GPL-3.0 and owes substantial conceptual credit to those authors.

What is this?
-------------
ldsc-cpp provides two subcommands:
- ldsc munge  — QC/convert GWAS summary statistics to a compact LDSC-style format.
- ldsc ph2    — partitioned heritability via two-step weighted regression on LD scores
                (supports per-chromosome or flat LD score files).

It outputs human-readable logs and a *.summary.txt with headline estimates.

Build
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

Prebuilt binaries (placeholder)
-------------------------------
Prebuilt archives for Linux, macOS (Intel/Apple), and Windows will be posted on the project’s
Releases page in the future.

Quick start with sample (public) GWAS data
------------------------------------------
The following is a generic walkthrough using publicly available resources. Replace paths with your filenames.

1) Prepare reference files (once)
   - w_hm3.snplist — HapMap3 SNP list.
   - eur_w_ld_chr/ — Per-chromosome LD score directory with files like:
       eur_w_ld_chr/
         1.l2.ldscore
         2.l2.ldscore
         ...
         22.l2.ldscore
         (optional) 1.l2.M_5_50 ... 22.l2.M_5_50

2) Munge GWAS
   Assume you’ve downloaded a tab-delimited GWAS file (e.g., from the GWAS Catalog) to gwas_raw.tsv
   with columns: SNP, A1, A2, P, and an effect column (e.g., BETA).

   ./ldsc munge      --sumstats gwas_raw.tsv      --out META_EUR      --id SNP      --p P      --A1 A1 --A2 A2      --eff_col BETA --eff_type beta      --merge-alleles w_hm3.snplist

   Result:
     META_EUR.sumstats  (columns: SNP N Z [A1 A2 FRQ])
     META_EUR.log

   Notes:
     - If your file has case/control counts or info/maf columns, you can pass:
         --N_cas_col N_CASES --N_con_col N_CONT
         --info INFO --maf MAF --keep-maf
     - If there’s no per-SNP N, you can provide a constant:
         --N 100000

3) Partitioned heritability

   Per-chromosome reference/weights:
     ./ldsc ph2        --out META_EUR_h2        --h2 META_EUR.sumstats        --ref-ld-chr eur_w_ld_chr        --w-ld-chr  eur_w_ld_chr

   Flat single-file mode (if you merged LD scores beforehand):
     ./ldsc ph2        --out META_EUR_h2        --h2 META_EUR.sumstats        --ref-ld ref.ldscore.tsv        --w-ld  wld.ldscore.tsv

   Result:
     META_EUR_h2.summary.txt  — top-line estimates (intercept, h2, SEs, QC metrics)
     META_EUR_h2.log          — detailed call/QC/merge info

Quick reference
---------------
ldsc munge (common flags):
  --sumstats <file>                (TSV with headers)
  --out <prefix>
  --id <col>                       SNP ID column (e.g., SNP)
  --p <col>                        p-value column
  --A1 <col> --A2 <col>            effect/other allele (omit with --no-alleles)
  --eff_col <col> --eff_type {Zscore|OR|beta|log_odds}
Optional QC/enrichment:
  --merge-alleles w_hm3.snplist
  --info <col> --info-min 0.9
  --maf <col> --maf-min 0.01 --keep-maf
Sample size:
  --N <scalar> OR --N_cas <scalar> --N_con <scalar>
  (or per-row --N_col, --N_cas_col, --N_con_col)

Outputs:
  <out>.sumstats
  <out>.log

ldsc ph2 (common flags):
  --out <prefix>
  --h2 <file.sumstats>
Exactly one of:
  --ref-ld-chr <dir | pattern-with-{}>     OR
  --ref-ld <flat-file>
Exactly one of:
  --w-ld-chr <dir | pattern-with-{}>       OR
  --w-ld <flat-file>
Optional:
  --M <comma-list | file>    (per-annotation M; otherwise auto-broadcasts M from *.l2.M_5_50 if present,
                              or uses SNP count)
  --n-blocks <int>           (jackknife blocks, default 200)
  --no-intercept             (fix intercept to 1)
  --intercept-h2 <value>     (freeze intercept to provided value)

Outputs:
  <out>.summary.txt
  <out>.log

Notes & differences vs. the original LDSC implementation
--------------------------------------------------------
- Two-step weighting with chi^2 cap at 30 for stability.
- Accepts per-chr directories OR {} patterns (e.g., ld/{}.l2.ldscore).
- If per-chr *.l2.M_5_50 files exist, their totals are summed and broadcast as M; otherwise M defaults to SNP count.
- OpenMP is optional; results should be statistically consistent single- vs multi-threaded.

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

