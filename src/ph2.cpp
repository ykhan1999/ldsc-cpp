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


// ph2.cpp — LDSC-style h2 via two-step GLS with wLD weights, χ² cap=30 for weighting,
// zero-variance annotation pruning, observed-scale h2 from M_5_50 and N_bar.
// Supports flat (--ref-ld/--w-ld) and per-chr (--ref-ld-chr/--w-ld-chr) inputs.
//
// Expected columns:
//   sumstats: SNP Z N
//   ref/w LD score tables: must include SNP; numeric columns auto-detected.
//
// Output:
//   <out>.ph2.summary.txt  (Intercept + each annotation slope + SE)
// Logs include: mean χ², λGC, max χ², and "Total Observed scale h2: v (se)"
//
// (C) 2025

#include "ph2.hpp"
#include "logger.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <functional>

// ============ Args ============
struct Ph2Args {
    std::string out;            // --out
    std::string sumstats;       // --h2 (clean .sumstats)
    std::string ref_ld;         // --ref-ld (single flat file)
    std::string ref_ld_chr;     // --ref-ld-chr (dir with 1..22.l2.ldscore OR {} pattern)
    std::string w_ld;           // --w-ld (single flat file)
    std::string w_ld_chr;       // --w-ld-chr (dir or pattern)
    std::string M;              // --M (optional: comma-list or single-line file)
    bool overlap_annot = false; // --overlap-annot (writes .results too)
    bool print_coefficients = false; // --print-coefficients (kept for parity)
    bool no_intercept = false;  // --no-intercept (LDSC semantics: fixes intercept to 1)
    double intercept_h2 = std::numeric_limits<double>::quiet_NaN(); // --intercept-h2
    int    n_blocks    = 200;   // --n-blocks (jackknife blocks)
};

static Ph2Args parse_ph2(int argc, char** argv) {
    Ph2Args a;
    auto need = [&](int& i, const std::string& k){
        if (i+1>=argc) { std::cerr<<"Missing value for "<<k<<"\n"; std::exit(2); }
        return std::string(argv[++i]);
    };
    for (int i=0;i<argc;i++){
        std::string k = argv[i];
        if      (k=="--out")            a.out = need(i,k);
        else if (k=="--h2")             a.sumstats = need(i,k);
        else if (k=="--ref-ld")         a.ref_ld = need(i,k);
        else if (k=="--ref-ld-chr")     a.ref_ld_chr = need(i,k);
        else if (k=="--w-ld")           a.w_ld = need(i,k);
        else if (k=="--w-ld-chr")       a.w_ld_chr = need(i,k);
        else if (k=="--M")              a.M = need(i,k);
        else if (k=="--overlap-annot")  a.overlap_annot = true;
        else if (k=="--print-coefficients") a.print_coefficients = true;
        else if (k=="--no-intercept")   a.no_intercept = true;
        else if (k=="--intercept-h2")   a.intercept_h2 = std::stod(need(i,k));
        else if (k=="--n-blocks")       a.n_blocks = std::stoi(need(i,k));
        else if (k.rfind("--",0)==0) { std::cerr<<"Unknown flag in ph2: "<<k<<"\n"; std::exit(2); }
    }
    if (a.out.empty())        { std::cerr<<"ph2: --out required\n"; std::exit(2); }
    if (a.sumstats.empty())   { std::cerr<<"ph2: --h2 required (clean .sumstats)\n"; std::exit(2); }
    if ((a.ref_ld.empty() == a.ref_ld_chr.empty())) {
        std::cerr<<"ph2: specify exactly one of --ref-ld (flat) or --ref-ld-chr (per-chr)\n"; std::exit(2);
    }
    if ((a.w_ld.empty() == a.w_ld_chr.empty())) {
        std::cerr<<"ph2: specify exactly one of --w-ld (flat) or --w-ld-chr (per-chr)\n"; std::exit(2);
    }
    if (a.n_blocks <= 1) {
        std::cerr<<"ph2: --n-blocks must be > 1\n"; std::exit(2);
    }
    if (a.no_intercept) a.intercept_h2 = 1.0; // LDSC semantics: fix intercept to 1
    return a;
}

static void echo_call(Logger& log, const Ph2Args& a) {
    std::ostringstream ss;
    ss << "Call (ph2): --out " << a.out
       << " --h2 " << a.sumstats
       << (a.ref_ld.empty() ? " --ref-ld-chr "+a.ref_ld_chr : " --ref-ld "+a.ref_ld)
       << (a.w_ld.empty()   ? " --w-ld-chr "+a.w_ld_chr     : " --w-ld "+a.w_ld);
    if (!a.M.empty()) ss << " --M " << a.M;
    if (a.overlap_annot) ss << " --overlap-annot";
    if (a.print_coefficients) ss << " --print-coefficients";
    if (!std::isnan(a.intercept_h2)) ss << " --intercept-h2 " << a.intercept_h2;
    if (a.n_blocks!=200) ss << " --n-blocks " << a.n_blocks;
    log.log(ss.str());
}

// ============ Small utilities ============
static std::vector<std::string> split_ws(const std::string& s){
    std::vector<std::string> v; std::string t; std::istringstream iss(s);
    while (iss>>t) v.push_back(t);
    return v;
}

static std::vector<std::string> expand_chr_input(const std::string& s, const std::string& suffix) {
    std::vector<std::string> out;
    auto pos = s.find("{}");
    if (pos != std::string::npos) {
        out.reserve(22);
        for (int c=1;c<=22;++c){
            std::ostringstream one; one<<s.substr(0,pos)<<c<<s.substr(pos+2);
            out.push_back(one.str());
        }
        return out;
    }
    if (std::filesystem::is_directory(s)) {
        out.reserve(22);
        for (int c=1;c<=22;++c){
            std::ostringstream one; one<<s<<"/"<<c<<suffix;
            out.push_back(one.str());
        }
        return out;
    }
    out.push_back(s); // single file
    return out;
}

struct Table {
    std::vector<std::string> snp;
    std::vector<std::string> cols;
    std::vector<std::vector<double>> data; // data[j][i] = column j, row i
};
struct Sumstats {
    std::vector<std::string> snp;
    std::vector<double> Z, N;
};

static Sumstats read_sumstats(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open sumstats: "+path);
    std::string line; if (!std::getline(in,line)) throw std::runtime_error("Empty sumstats: "+path);
    auto hdr = split_ws(line);
    int iS=-1, iZ=-1, iN=-1;
    for (int i=0;i<(int)hdr.size();++i){
        if (hdr[i]=="SNP") iS=i; else if (hdr[i]=="Z") iZ=i; else if (hdr[i]=="N") iN=i;
    }
    if (iS<0||iZ<0||iN<0) throw std::runtime_error("sumstats must have SNP Z N columns");
    Sumstats ss;
    while (std::getline(in,line)){
        if (line.empty()) continue;
        auto f = split_ws(line);
        if ((int)f.size() <= std::max(iS,std::max(iZ,iN))) continue;
        ss.snp.push_back(f[iS]);
        ss.Z.push_back(std::strtod(f[iZ].c_str(), nullptr));
        ss.N.push_back(std::strtod(f[iN].c_str(), nullptr));
    }
    return ss;
}

using ColSelector = std::function<std::vector<int>(const std::vector<std::string>&)>;

static std::vector<int> pick_ref_cols(const std::vector<std::string>& hdr){
    static const std::unordered_set<std::string> skip = {
        "SNP","CHR","BP","CM","MAF","A1","A2","L2_MAF"
    };
    std::vector<int> idx;
    for (int i=0;i<(int)hdr.size();++i){
        if (!skip.count(hdr[i]) && hdr[i]!="SNP") idx.push_back(i);
    }
    if (idx.empty()) throw std::runtime_error("ref-LD has no usable annotation columns.");
    return idx;
}

static std::vector<int> pick_w_cols(const std::vector<std::string>& hdr){
    static const std::unordered_set<std::string> meta = {"SNP","CHR","BP","CM","MAF","A1","A2","L2_MAF"};
    const char* prefer[] = {"L2","LDSCORE","L2_no_mhc"};
    for (const char* p : prefer){
        for (int i=0;i<(int)hdr.size();++i) if (hdr[i]==p) return {i};
    }
    int best=-1;
    for (int i=0;i<(int)hdr.size();++i) if (!meta.count(hdr[i]) && hdr[i]!="SNP") best=i;
    if (best<0) throw std::runtime_error("w-LD has no usable numeric column.");
    return {best};
}

static int find_col(const std::vector<std::string>& hdr, const std::string& name){
    for (int i=0;i<(int)hdr.size();++i) if (hdr[i]==name) return i;
    return -1;
}

static Table read_ld_one_select(const std::string& path, const ColSelector& pick) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open LD file: "+path);
    std::string line; if (!std::getline(in,line)) throw std::runtime_error("Empty LD file: "+path);
    auto hdr = split_ws(line);
    int iSNP = find_col(hdr, "SNP");
    if (iSNP<0) throw std::runtime_error("LD file missing SNP header: "+path);

    auto idx = pick(hdr);

    Table t;
    t.cols.reserve(idx.size());
    for (int j : idx) t.cols.push_back(hdr[j]);
    t.data.resize(idx.size());

    std::string line2;
    while (std::getline(in,line2)){
        if (line2.empty()) continue;
        auto f = split_ws(line2);
        if ((int)f.size() <= iSNP) continue;
        t.snp.push_back(f[iSNP]);
        for (size_t k=0;k<idx.size();++k){
            int j = idx[k];
            double v = (j<(int)f.size()) ? std::strtod(f[j].c_str(), nullptr)
                                         : std::numeric_limits<double>::quiet_NaN();
            t.data[k].push_back(v);
        }
    }
    return t;
}

static Table read_ld_fileset_any(const std::string& single_or_dir_or_pattern,
                                 const std::string& suffix,
                                 const ColSelector& pick) {
    Table agg;
    auto files = expand_chr_input(single_or_dir_or_pattern, suffix);
    bool first = true;
    for (const auto& f: files) {
        auto t = read_ld_one_select(f, pick);
        if (first) {
            agg = std::move(t);
            first = false;
        } else {
            if (t.cols != agg.cols) throw std::runtime_error("LD column sets differ across files: "+f);
            agg.snp.insert(agg.snp.end(), t.snp.begin(), t.snp.end());
            for (size_t j=0;j<agg.data.size();++j){
                agg.data[j].insert(agg.data[j].end(), t.data[j].begin(), t.data[j].end());
            }
        }
    }
    return agg;
}

static Table read_ld_flat(const std::string& path, const ColSelector& pick){
    return read_ld_one_select(path, pick);
}

// ============ Auto M_5_50 ============
static long long read_M_total_from_dir_or_pattern(const std::string& s){
    if (!std::filesystem::is_directory(s)) return -1;
    long long total = 0;
    for (int c=1;c<=22;++c){
        std::ostringstream p; p<<s<<"/"<<c<<".l2.M_5_50";
        std::ifstream in(p.str());
        if (!in) throw std::runtime_error("Could not open M_5_50 file: "+p.str());
        long long m_chr = 0;
        if (!(in>>m_chr)) throw std::runtime_error("Could not parse integer from "+p.str());
        total += m_chr;
    }
    return total;
}

// ============ Merge ============
struct MergeResult {
    std::vector<std::string> snp;
    std::vector<double> y_chisq;   // Z^2
    std::vector<double> N;
    std::vector<std::string> annot_names;
    std::vector<std::vector<double>> annot; // k x n
    std::vector<double> wld; // regression weight LD (one column)
};

static MergeResult merge_sumstats_ld(const Sumstats& ss, const Table& refld, const Table& wld, Logger& log) {
    std::unordered_map<std::string,int> idx_ref, idx_w;
    idx_ref.reserve(refld.snp.size()*2);
    idx_w.reserve(wld.snp.size()*2);
    for (int i=0;i<(int)refld.snp.size();++i) idx_ref.emplace(refld.snp[i], i);
    for (int i=0;i<(int)wld.snp.size();++i)   idx_w.emplace(wld.snp[i], i);

    MergeResult r;
    r.annot_names = refld.cols;
    r.annot.resize(refld.cols.size());
    size_t dropped = 0;
    for (size_t i=0;i<ss.snp.size();++i){
        auto it1 = idx_ref.find(ss.snp[i]);
        auto it2 = idx_w.find(ss.snp[i]);
        if (it1==idx_ref.end() || it2==idx_w.end()) { dropped++; continue; }
        int ir = it1->second, iw = it2->second;
        r.snp.push_back(ss.snp[i]);
        double z = ss.Z[i];
        r.y_chisq.push_back(z*z);
        r.N.push_back(ss.N[i]);
        for (size_t j=0;j<refld.cols.size();++j) r.annot[j].push_back(refld.data[j][ir]);
        if (wld.cols.size()!=1) throw std::runtime_error("--w-ld/--w-ld-chr must yield exactly one weight column.");
        r.wld.push_back(wld.data[0][iw]);
    }
    if (r.snp.empty()) throw std::runtime_error("After merging sumstats with LD files, no SNPs remain.");
    log.log("After merging, ", r.snp.size(), " SNPs remain (dropped ", dropped, ").");
    return r;
}

// ============ Linear algebra / WLS / Jackknife ============
static void xtwx_xtwy(const std::vector<std::vector<double>>& X, const std::vector<double>& w,
                      const std::vector<double>& y, std::vector<std::vector<double>>& XtWX,
                      std::vector<double>& XtWy) {
    int p = (int)X.size();
    int n = (int)w.size();
    XtWX.assign(p, std::vector<double>(p, 0.0));
    XtWy.assign(p, 0.0);
    for (int i=0;i<n;++i){
        double wi = w[i];
        for (int a=0;a<p;++a){
            double xa = X[a][i];
            XtWy[a] += xa * wi * y[i];
            for (int b=a;b<p;++b){
                XtWX[a][b] += xa * wi * X[b][i];
            }
        }
    }
    for (int a=0;a<p;++a) for (int b=0;b<a;++b) XtWX[a][b] = XtWX[b][a];
}

static std::vector<double> solve_spd(std::vector<std::vector<double>> A, std::vector<double> b){
    int p = (int)A.size();
    for (int i=0;i<p;++i) A[i][i] += 1e-12;
    for (int i=0;i<p;++i){
        int piv=i; double best=std::fabs(A[i][i]);
        for (int r=i+1;r<p;++r){ double v=std::fabs(A[r][i]); if (v>best){ best=v; piv=r; } }
        if (best==0) throw std::runtime_error("Singular normal matrix.");
        if (piv!=i){ std::swap(A[piv],A[i]); std::swap(b[piv],b[i]); }
        double d = A[i][i];
        for (int j=i;j<p;++j) A[i][j] /= d; b[i] /= d;
        for (int r=0;r<p;++r){
            if (r==i) continue;
            double f = A[r][i];
            if (f==0) continue;
            for (int j=i;j<p;++j) A[r][j] -= f*A[i][j];
            b[r] -= f*b[i];
        }
    }
    return b;
}

struct WLSResult { std::vector<double> beta; };
static WLSResult wls(const std::vector<std::vector<double>>& X, const std::vector<double>& w,
                     const std::vector<double>& y) {
    std::vector<std::vector<double>> XtWX;
    std::vector<double> XtWy;
    xtwx_xtwy(X, w, y, XtWX, XtWy);
    WLSResult r;
    r.beta = solve_spd(XtWX, XtWy);
    return r;
}

static std::vector<int> block_boundaries(int n, int B){
    B = std::max(2, std::min(B, n));
    std::vector<int> cuts(B+1,0);
    for (int b=0;b<=B;++b) cuts[b] = (int)std::llround((long double)b * n / B);
    return cuts;
}

static void jackknife(const std::vector<std::vector<double>>& X, const std::vector<double>& w,
                      const std::vector<double>& y, int B,
                      std::vector<double>& beta, std::vector<double>& se) {
    int n = (int)w.size();
    int p = (int)X.size();
    auto full = wls(X, w, y);
    beta = full.beta;

    auto cuts = block_boundaries(n, B);
    std::vector<std::vector<double>> betas; betas.reserve(B);
    for (int b=0;b<B;++b){
        int lo = cuts[b], hi = cuts[b+1];
        int keep = n - (hi-lo);
        std::vector<double> wy; wy.reserve(keep);
        std::vector<double> ww; ww.reserve(keep);
        std::vector<std::vector<double>> WX(p, std::vector<double>()); 
        for (int a=0;a<p;++a) WX[a].reserve(keep);
        for (int i=0;i<n;++i){
            if (i>=lo && i<hi) continue;
            ww.push_back(w[i]); wy.push_back(y[i]);
            for (int a=0;a<p;++a) WX[a].push_back(X[a][i]);
        }
        auto fit = wls(WX, ww, wy);
        betas.push_back(std::move(fit.beta));
    }
    se.assign(p, 0.0);
    for (int a=0;a<p;++a){
        long double mean=0.0L;
        for (int b=0;b<B;++b) mean += betas[b][a];
        mean /= B;
        long double var=0.0L;
        for (int b=0;b<B;++b){ long double d = betas[b][a] - mean; var += d*d; }
        var *= (B-1.0L) / B;
        se[a] = std::sqrt((double)var);
    }
}

// ============ Helpers for LDSC weighting & h2 ============
static double vec_mean(const std::vector<double>& v){
    long double s=0.0L; for (double x: v) s += x;
    return (v.empty()? 0.0 : (double)(s / (long double)v.size()));
}

// Build X with (intercept?, annotations...)
static void build_design(const MergeResult& R, bool include_intercept,
                         std::vector<std::vector<double>>& X){
    int n = (int)R.snp.size();
    int k = (int)R.annot.size();
    if (include_intercept){
        X.assign(k+1, std::vector<double>(n, 1.0));
        for (int j=0;j<k;++j) X[j+1] = R.annot[j];
    } else {
        X.assign(k, std::vector<double>(n, 0.0));
        for (int j=0;j<k;++j) X[j] = R.annot[j];
    }
}

// compute observed-scale h2 and its SE from slopes & SEs: h2 = sum_j (M_j/N_bar)*beta_j
static void h2_from_betas(const std::vector<double>& beta, const std::vector<double>& se,
                          bool have_intercept, const std::vector<double>& M,
                          double N_bar, double& h2, double& h2_se){
    int offs = have_intercept ? 1 : 0;
    int k = (int)M.size();
    h2 = 0.0; long double var = 0.0L;
    for (int j=0;j<k;++j){
        double c = M[j]/N_bar;
        h2 += c * beta[j+offs];
        var += (long double)c*c * (long double)se[j+offs]* (long double)se[j+offs];
    }
    h2_se = std::sqrt((double)var);
}

// Build LDSC step-1 weights: w_i = 1 / wld_i
static void weights_step1(const MergeResult& R, std::vector<double>& w){
    int n = (int)R.snp.size();
    w.resize(n);
    for (int i=0;i<n;++i){
        double v = R.wld[i];
        if (!(v>0)) v = 1.0;
        w[i] = 1.0 / v;
    }
}

// Build LDSC step-2 weights using step-1 betas.
// pred_term_i = 1 + sum_j beta1_j * l_ij  (cap the implied χ² at 30 in weighting)
// w_i = 1 / ( wld_i * pred_term_i^2 )
static void weights_step2(const MergeResult& R,
                          const std::vector<double>& beta1, bool have_intercept,
                          std::vector<double>& w){
    int n = (int)R.snp.size();
    int k = (int)R.annot.size();
    int offs = have_intercept ? 1 : 0;
    w.resize(n);
    for (int i=0;i<n;++i){
        // linear predictor from annotations only (exclude intercept)
        long double lp = 0.0L;
        for (int j=0;j<k;++j){
            lp += (long double)beta1[j+offs] * (long double)R.annot[j][i];
        }
        // expected chi^2 ≈ intercept + lp; LDSC caps χ² at 30 for weighting.
        // We emulate by capping the total expectation to 30.
        long double exp_chi2 = (have_intercept ? (long double)beta1[0] : 1.0L) + lp;
        if (exp_chi2 < 1.0L) exp_chi2 = 1.0L;
        if (exp_chi2 > 30.0L) exp_chi2 = 30.0L;

        // weights use (1 + polygenic term) ≈ exp_chi2, so:
        long double denom_scale = exp_chi2; // equivalent to (1 + N h2 L2 / M) in 1-annot case
        double base = R.wld[i]; if (!(base>0)) base = 1.0;
        long double denom = (long double)base * denom_scale * denom_scale;
        if (!(denom>0)) denom = 1.0L;
        w[i] = 1.0 / (double)denom;
    }
}

// ============ ph2 main ============
int run_ph2(int argc, char** argv) {
    Ph2Args args = parse_ph2(argc, argv);
    Logger log(args.out + ".log", /*append=*/true);
    echo_call(log, args);

    try {
        // 1) Sumstats
        log.log("Reading sumstats from ", args.sumstats, " ...");
        Sumstats ss = read_sumstats(args.sumstats);
        log.log("Read ", (int)ss.snp.size(), " SNPs from sumstats.");

        // 2) LD files
        log.log("Reading reference LD ...");
        Table refld = args.ref_ld.empty()
            ? read_ld_fileset_any(args.ref_ld_chr, ".l2.ldscore", pick_ref_cols)
            : read_ld_flat(args.ref_ld, pick_ref_cols);
        log.log("Reference LD: ", refld.snp.size(), " rows, ", (int)refld.cols.size(), " columns.");

        log.log("Reading regression weights (w-LD) ...");
        Table wld = args.w_ld.empty()
            ? read_ld_fileset_any(args.w_ld_chr, ".l2.ldscore", pick_w_cols)
            : read_ld_flat(args.w_ld, pick_w_cols);
        if (wld.cols.size()!=1) throw std::runtime_error("--w-ld/--w-ld-chr must yield exactly one weight column.");

        // 3) Merge
        MergeResult R = merge_sumstats_ld(ss, refld, wld, log);

        // 4) M vector
        std::vector<double> M(R.annot.size(), 0.0);
        if (!args.M.empty()){
            std::vector<double> M_in;
            if (args.M.find(',')!=std::string::npos){
                std::istringstream iss(args.M); std::string tok;
                while (std::getline(iss,tok,',')) M_in.push_back(std::strtod(tok.c_str(), nullptr));
            } else {
                std::ifstream mf(args.M);
                if (!mf) throw std::runtime_error("Could not open --M: "+args.M);
                std::string line; if (!std::getline(mf,line)) throw std::runtime_error("Empty --M file");
                std::istringstream iss(line); double v; while (iss>>v) M_in.push_back(v);
            }
            if ((int)M_in.size() != (int)R.annot.size())
                throw std::runtime_error("--M length does not match # of LD columns after pruning.");
            M = std::move(M_in);
            log.log("Using user-supplied M vector with length ", (int)M.size(), ".");
        } else {
            long long M_total = read_M_total_from_dir_or_pattern(args.ref_ld_chr);
            if (M_total > 0) {
                for (double &x : M) x = (double)M_total;
                log.log("Auto M_5_50 total (sum over chr) = ", (long long)M_total, " ; broadcasting to ", (int)M.size(), " annotations.");
            } else {
                for (double &x : M) x = (double)R.snp.size();
                log.log("No per-chromosome M_5_50 found; using merged SNP count (", (int)R.snp.size(), ") as M proxy.");
            }
        }

        // 5) Zero-variance annotation pruning (after M prepared; we’ll keep M aligned)
        {
            // If we drop columns, mirror-drop in M as well
            std::vector<int> keep;
            for (size_t j=0;j<R.annot.size();++j){
                const auto& col = R.annot[j];
                double mu=0.0; for (double v: col) mu += v; mu /= std::max<size_t>(1,col.size());
                long double var=0.0L; for (double v: col){ long double d=v-mu; var += d*d; }
                if (var>0) keep.push_back((int)j);
            }
            if (keep.empty()) throw std::runtime_error("All LD Score columns have zero variance.");
            if ((int)keep.size() != (int)R.annot.size()) {
                Logger tmp(args.out + ".log", /*append=*/true);
                tmp.log("Removing ", (int)R.annot.size() - (int)keep.size(), " partitioned LD columns with zero variance.");
                std::vector<std::string> names2; names2.reserve(keep.size());
                std::vector<std::vector<double>> annot2; annot2.reserve(keep.size());
                std::vector<double> M2; M2.reserve(keep.size());
                for (int j: keep){ names2.push_back(R.annot_names[j]); annot2.push_back(std::move(R.annot[j])); M2.push_back(M[j]); }
                R.annot_names.swap(names2);
                R.annot.swap(annot2);
                M.swap(M2);
            }
        }

        // 6) Some QC prints (from raw y)
        {
            const auto& y0 = R.y_chisq;
            long double mean=0.0L; for (double v: y0) mean += v; mean/=std::max(1,(int)y0.size());
            std::vector<double> tmp=y0;
            std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end());
            double med = tmp[tmp.size()/2];
            double mx = *std::max_element(y0.begin(), y0.end());
            log.log("Mean chi^2 = ", (double)mean);
            log.log("Lambda GC = ", med/0.4549);
            log.log("Max chi^2 = ", mx);
        }

        // 7) Two-step estimator with χ² cap for weighting
        log.log("Using two-step estimator with cutoff at 30.");

        const bool include_intercept = std::isnan(args.intercept_h2); // default: free intercept
        const double N_bar = vec_mean(R.N);

        // STEP 1: weights = 1 / wLD
        std::vector<double> w1; weights_step1(R, w1);

        std::vector<std::vector<double>> X1;
        build_design(R, include_intercept, X1);

        // If intercept is fixed (no_intercept / intercept-h2), subtract it from y as per LDSC semantics.
        std::vector<double> y = R.y_chisq;
        if (!include_intercept) {
            // intercept_h2 given: subtract that constant from y
            for (double& v: y) v -= args.intercept_h2;
        }

        // STEP 1 regression (jackknife for SE won’t be used; we just need betas)
        std::vector<double> beta1, se1;
        jackknife(X1, w1, y, std::min(args.n_blocks, (int)R.snp.size()), beta1, se1);

        // STEP 2: reweight using step-1 betas with χ² cap=30
        std::vector<double> w2; weights_step2(R, beta1, include_intercept, w2);

        // STEP 2 regression (final): get betas & SEs
        std::vector<double> beta, se;
        jackknife(X1, w2, y, std::min(args.n_blocks, (int)R.snp.size()), beta, se);

        // 8) Convert slopes to observed-scale h2 (sums over annotations)
        double h2_obs=0.0, h2_se=0.0;
        h2_from_betas(beta, se, include_intercept, M, N_bar, h2_obs, h2_se);

        // 10) Log summary
        {
            const auto& y0 = R.y_chisq;
            long double mean=0.0L; for (double v: y0) mean += v; mean/=std::max(1,(int)y0.size());
            std::vector<double> tmp=y0;
            std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end());
            double med = tmp[tmp.size()/2];
            double mx = *std::max_element(y0.begin(), y0.end());

            log.log("Total Observed scale h2: ", h2_obs, " (", h2_se, ")");
            log.log("Lambda GC: ", med/0.4549);
            log.log("Mean Chi^2: ", (double)mean);
            if (include_intercept) log.log("Intercept: ", beta[0], " (", se[0], ")");
            log.log("Max chi^2: ", mx);
        }

        // write a compact summary file with top-line estimates
        {
            // Recompute quick QC stats from original chi^2 vector R.y_chisq
            const auto& y_all = R.y_chisq;
            long double mean_all = 0.0L; for (double v: y_all) mean_all += v;
            mean_all /= std::max(1,(int)y_all.size());
            std::vector<double> tmp_all = y_all;
            std::nth_element(tmp_all.begin(), tmp_all.begin()+tmp_all.size()/2, tmp_all.end());
            double med_all = tmp_all[tmp_all.size()/2];
            double mx_all = *std::max_element(y_all.begin(), y_all.end());

            // Figure out reported intercept value/SE (fixed vs free)
            const bool include_intercept = std::isnan(args.intercept_h2);
            double intercept_est = include_intercept ? beta[0] : args.intercept_h2;
            double intercept_se  = include_intercept ? se[0]   : 0.0;

            // Derive an M_5_50 total to print (if broadcast, M[0] is the total; else sum)
            long double M_total_print = 0.0L;
            if (!M.empty()) {
                bool all_same = true;
                for (size_t j=1;j<M.size();++j) if (M[j]!=M[0]) { all_same=false; break; }
                if (all_same) M_total_print = M[0];
                else { for (double v: M) M_total_print += v; }
            }

            std::ofstream s(args.out + ".summary.txt");
            if (!s) throw std::runtime_error("Could not open summary output file.");

            s << "Field\tValue\n";
            s << "SNPs_used\t" << R.snp.size() << "\n";
            s << "Annotations\t" << R.annot_names.size() << "\n";
            s << "N_bar\t" << N_bar << "\n";
            s << "M_5_50_total\t" << (double)M_total_print << "\n";
            s << "Intercept\t" << intercept_est << "\n";
            s << "Intercept_SE\t" << intercept_se << "\n";
            s << "h2_observed\t" << h2_obs << "\n";
            s << "h2_observed_SE\t" << h2_se << "\n";
            s << "Mean_Chi2\t" << (double)mean_all << "\n";
            s << "Lambda_GC\t" << (med_all/0.4549) << "\n";
            s << "Max_Chi2\t" << mx_all << "\n";
            s.close();

            log.log("Wrote summary to ", args.out, ".summary.txt");
        }

        //finish
        log.log("ph2: completed.");
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ph2 ERROR: " << e.what() << "\n";
        return 1;
    }
}
