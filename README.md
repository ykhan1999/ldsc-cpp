
# ldsc-cpp

> **Acknowledgement:** This project is an independent C++ reimplementation inspired by the original **LDSC (LD Score Regression)** developed by Brendan K. Bulik-Sullivan, Hilary K. Finucane, Po-Ru Loh, Alkes Price, and collaborators. The original LDSC is licensed under GPL-3.0 and is available at the Price Lab’s repository. This project follows GPL-3.0 and owes substantial conceptual credit to those authors.

What is this?
-------------
ldsc-cpp provides two subcommands:
- ldsc munge  — QC/convert GWAS summary statistics to a compact LDSC-style format.
- ldsc ph2    — partitioned heritability via two-step weighted regression on LD scores
                (supports per-chromosome or flat LD score files).

It outputs human-readable logs and a *.summary.txt with headline estimates.

Getting started: Prebuilt binaries
-------------------------------
Prebuilt archives for Linux, macOS (Intel/Apple), and Windows will be posted on the project’s
Releases page in the future.

Getting started: Building from Source
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

Quick start with sample GWAS data
------------------------------------------

1) Prepare reference files (once)
   - LD score files, which can be flat files or organized per-chromosome. Sample LD score files computed precomputed from 1000 Genomes European ancestry samples:
		```bash
		wget https://y3782016.eero.online/eur_w_ld_chr.tar.gz
		tar -xzvf eur_w_ld_chr.tar.gz
		```
   - HapMap3 SNP list
	   - This is a reference of high-quality SNPs to include for LDSC computations, which will be specific to your LD score panel. A sample file of well-imputed HapMap3 SNPs in European-ancestry-like individuals is below:
		```bash
		wget https://y3782016.eero.online/w_hm3.snplist
		```


2) Munge GWAS
   Using sample data from GWAS Catalog:
	```bash
	wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006901/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz
	gzip -d Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz
   ./ldsc munge \
   --sumstats Meta-analysis_Wood_et_al+UKBiobank_2018.txt \
   --out meta_height \
   --id SNP \
   --p P \
   --A1 Tested_Allele \
   --A2 Other_Allele \
   --eff_col BETA \
   --eff_type beta \
   --merge-alleles w_hm3.snplist \
     --N 700000
	```
	
   Result:
     META_height.sumstats  (columns: SNP N Z A1 A2)
     META_height.log

   Notes:
     - If your file has case/control counts or info/maf columns, you can pass:
		```bash
		--N_cas_col N_CASES --N_con_col N_CONT
		--info INFO --maf MAF --keep-maf
		```
     - If there’s no per-SNP N, you can provide a constant:
		```bash
		--N 700000
		```

3) Partitioned heritability



   Per-chromosome reference/weights:
      ```bash
		./ldsc ph2 --out \
		meta_height_h2 \
		--h2 meta_height.sumstats \
		--ref-ld-chr eur_w_ld_chr \
		--w-ld-chr eur_w_ld_chr
      ```

   Flat single-file mode (if you have merged LD scores):
      ```bash
		./ldsc ph2 --out \
		meta_height_h2 \
		--h2 meta_height.sumstats \
		--ref-ld ref.ldscore.tsv \
		--w-ld  wld.ldscore.tsv
      ```

   Result:
     meta_height_h2.summary.txt — top-line estimates (intercept, h2, SEs, QC metrics)
     meta_height_h2.log          — detailed call/QC/merge info

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


