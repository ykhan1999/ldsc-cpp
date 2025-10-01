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

// main.cpp — C++17 LDSC-style munge with explicit columns and --merge-alleles support.
// No external deps. Outputs TSV (uncompressed) "<out>.sumstats".

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// Avoid ODR clash with the global Logger used by ph2.cpp
#define Logger MungeLogger
using namespace std;

// ---------- Thread-safe logger ----------
struct Logger {
    explicit Logger(const std::string& path, bool append=false)
        : ofs(path, append ? std::ios::app : std::ios::out)
    {
        if (!ofs) {
            lock_guard<mutex> lk(mu());
            cerr << "WARN: could not open log file '" << path << "'; logging to stderr\n";
        }
    }
    template <class... Ts>
    void log(const Ts&... parts) {
        ostringstream oss;
        (oss << ... << parts);
        oss << '\n';
        lock_guard<mutex> lk(mu());
        if (ofs) ofs << oss.str(); else cerr << oss.str();
    }
  private:
    static mutex& mu(){ static mutex m; return m; }
    ofstream ofs;
};

// -------- Output syntax if help is called ----
#ifndef MUNGE_HELP
#   define MUNGE_HELP \
    "This is the main subcommand used to prepare sumstats for ldsc computations:\n" \
    "Usage: ./ldsc munge --[flag] [string]\n" \
    "   REQUIRED FLAGS\n" \
    "       --sumstats sumstats.txt (NOTE .gz file not supported)\n" \
    "       --id [SNP IDs column name]\n" \
    "       --p [p-value column name]\n" \
    "       --A1 [A1 column name]\n" \
    "       --A2 [A2 column name]\n" \
    "       --A1 [A1 column name]\n" \
    "       --eff-col [column with effect values]\n" \
    "       --eff-type [type of effect value; options: beta, OR, Zscore, log_odds]\n" \
    "       --out [output-prefix]\n" \
    "       Either of:\n" \
    "           --merge-alleles snplist.txt (known good snplist with cols SNP A1 A2)\n" \
    "           --no-alleles true\n" \
    "       Either of:\n" \
    "           --N-col [N column name]\n" \
    "           --N [raw N value]\n" \
    "   OPTIONAL FLAGS\n" \
    "       --maf [MAF column name] (filter by MAF specified with --maf-min)\n" \
    "       --maf-min [min MAF value]\n" \
    "       --info [INFO column name] (filter by INFO specified with --info-min)\n" \
    "       --info-min [min INFO value]\n" \
    "       --N-cas-col [column with N cases]\n" \
    "       --N-con-col [column with N controls]\n" \
    "       --N-cas [raw N cases value]\n" \
    "       --N-con [raw N controls value]\n"
#endif

// ---------- CLI ----------
struct Args {
    string sumstats;      // --sumstats
    string out;           // --out (prefix)
    string id_col;        // --id
    string P_col;         // --p
    string A1_col;        // --A1
    string A2_col;        // --A2
    string N_col;         // --N-col (optional)
    string N_cas_col;     // --N-cas-col
    string N_con_col;     // --N-con-col
    string stats_col;     // --eff-col
    string stats_type;    // --eff-type {Zscore, OR, beta, log_odds}
    string info_col;      // --info
    string maf_col;       // --maf
    string merge_alleles; // --merge-alleles (file with SNP A1 A2)
    double info_min = 0.9;  // --info-min
    double maf_min  = 0.01; // --maf-min
    bool   no_alleles = false; // --no-alleles
    bool   keep_maf   = false; // --keep-maf
    double N_flag = numeric_limits<double>::quiet_NaN(); // --N
    double N_cas  = numeric_limits<double>::quiet_NaN(); // --N-cas
    double N_con  = numeric_limits<double>::quiet_NaN(); // --N-con
};

static Args parse_args(int argc, char** argv) {
    if (argc <=1 ){
        std::cout << MUNGE_HELP;
        std::exit(2);
    }
    Args a;
    auto need = [&](int& i, const string& k){ if (i+1>=argc) { cerr<<"Missing value for "<<k<<"\n"; exit(2);} return string(argv[++i]); };
    for (int i=0;i<argc;i++){
        string k = argv[i];
        if      (k=="--sumstats") a.sumstats = need(i,k);
        else if (k=="--out")      a.out = need(i,k);
        else if (k=="--id")       a.id_col = need(i,k);
        else if (k=="--p")        a.P_col = need(i,k);
        else if (k=="--A1")       a.A1_col = need(i,k);
        else if (k=="--A2")       a.A2_col = need(i,k);
        else if (k=="--N-col")    a.N_col = need(i,k);
        else if (k=="--N-cas-col")  a.N_cas_col = need(i,k);
        else if (k=="--N-con-col")   a.N_con_col = need(i,k);
        else if (k=="--eff-col")  a.stats_col  = need(i,k);
        else if (k=="--eff-type") a.stats_type = need(i,k);
        else if (k=="--info")     a.info_col = need(i,k);
        else if (k=="--maf")      a.maf_col  = need(i,k);
        else if (k=="--merge-alleles") a.merge_alleles = need(i,k);
        else if (k=="--info-min") a.info_min = stod(need(i,k));
        else if (k=="--maf-min")  a.maf_min  = stod(need(i,k));
        else if (k=="--no-alleles") a.no_alleles = true;
        else if (k=="--keep-maf")   a.keep_maf   = true;
        else if (k=="--N")       a.N_flag = stod(need(i,k));
        else if (k=="--N-cas")   a.N_cas  = stod(need(i,k));
        else if (k=="--N-con")   a.N_con  = stod(need(i,k));
        else if (k=="--help") {
            std::cout << MUNGE_HELP;
            std::exit(2);
        }
        else {
            std::cerr<<"Unknown flag in munge: "<<k<<"\n"<<MUNGE_HELP;
            std::exit(2);
        }
    }
    // required
    if (a.sumstats.empty() || a.out.empty() || a.id_col.empty() || a.P_col.empty()) {
        cerr<<"Required: --sumstats --out --id --p\n"; exit(2);
    }
    if (!a.no_alleles && (a.A1_col.empty() || a.A2_col.empty())) {
        cerr<<"Provide --A1 and --A2 or use --no-alleles\n"; exit(2);
    }
    if (!a.merge_alleles.empty() && a.no_alleles) {
        cerr<<"--merge-alleles is incompatible with --no-alleles\n"; exit(2);
    }
    if (a.stats_col.empty() || a.stats_type.empty()) {
        cerr<<"Provide --eff_col and --eff_type {Zscore, OR, beta, log_odds}\n"; exit(2);
    }
    if (std::isnan(a.N_flag) && a.N_col.empty()) {
        cerr<<"Provide one of either --N [total N] or --N-col [column with variant-level N]\n"; exit(2);
    }
    for (auto& c: a.stats_type) c = tolower(static_cast<unsigned char>(c));
    if (a.stats_type!="zscore" && a.stats_type!="or" && a.stats_type!="beta" && a.stats_type!="log_odds") {
        cerr<<"--eff_type must be one of: Zscore, OR, beta, log_odds\n"; exit(2);
    }
    return a;
}

// ---------- Utils ----------
static vector<string> split_ws(const string& s){ vector<string> v; string t; istringstream iss(s); while(iss>>t) v.push_back(t); return v; }
static bool is_valid_base(char c){ c=toupper(static_cast<unsigned char>(c)); return c=='A'||c=='C'||c=='G'||c=='T'; }
static inline char comp_base(char b){ switch(toupper(static_cast<unsigned char>(b))){case 'A':return 'T';case 'T':return 'A';case 'C':return 'G';case 'G':return 'C';default:return '?';}}
static bool is_strand_ambiguous(char a, char b){
    a=toupper(static_cast<unsigned char>(a)); b=toupper(static_cast<unsigned char>(b));
    if (a==b) return true;
    return (a=='A'&&b=='T')||(a=='T'&&b=='A')||(a=='C'&&b=='G')||(a=='G'&&b=='C');
}
static bool keep_alleles(const string& A1, const string& A2){
    if (A1.size()!=1 || A2.size()!=1) return false;
    if (!is_valid_base(A1[0]) || !is_valid_base(A2[0])) return false;
    return !is_strand_ambiguous(A1[0],A2[0]);
}
static bool pass_info(double x, double minv){ if (std::isnan(x)) return true; if (x<0||x>2.0) return false; return x>=minv; }
static bool pass_maf(double& f, double minv){ if (std::isnan(f)) return true; if (f<0||f>1) return false; f = min(f,1.0-f); return f>minv; }

// ---- P-value sanitization to avoid infinities (p=0) and nonsense (>1) -------
static inline double sanitize_p(double p) {
    if (std::isnan(p)) return std::numeric_limits<double>::quiet_NaN();
    if (p <= 0.0) return 1e-323; // floor to avoid inf z/chi2
    if (p > 1.0)  return 1.0;    // cap
    return p;
}
struct PClampStats { uint64_t zero_or_neg = 0, over_one = 0; } g_pstats;
static inline double sanitize_p_counting(double p) {
    if (std::isnan(p)) return p;
    if (p <= 0.0) { ++g_pstats.zero_or_neg; return 1e-323; }
    if (p >  1.0) { ++g_pstats.over_one;   return 1.0; }
    return p;
}

// Acklam invnorm for P->Z (Z = Phi^{-1}(1 - p/2))
static double inv_norm_cdf(double p){
    if (p<=0) return -numeric_limits<double>::infinity();
    if (p>=1) return  numeric_limits<double>::infinity();
    static const double a1=-3.969683028665376e+01, a2=2.209460984245205e+02, a3=-2.759285104469687e+02,
                        a4=1.383577518672690e+02, a5=-3.066479806614716e+01, a6=2.506628277459239e+00;
    static const double b1=-5.447609879822406e+01, b2=1.615858368580409e+02, b3=-1.556989798598866e+02,
                        b4=6.680131188771972e+01, b5=-1.328068155288572e+01;
    static const double c1=-7.784894002430293e-03, c2=-3.223964580411365e-01, c3=-2.400758277161838e+00,
                        c4=-2.549732539343734e+00, c5= 4.374664141464968e+00, c6= 2.938163982698783e+00;
    static const double d1= 7.784695709041462e-03, d2= 3.224671290700398e-01, d3= 2.445134137142996e+00,
                        d4= 3.754408661907416e+00;
    double q,r,x;
    if (p<0.02425){ q=sqrt(-2*log(p));
        x=((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6; x/=((((d1*q+d2)*q+d3)*q+d4)*q+1);
    } else if (p>1-0.02425){ q=sqrt(-2*log(1-p));
        x=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1);
    } else { q=p-0.5; r=q*q;
        x=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1);
    }
    // one Halley step
    double u = 0.5*(1+erf(x/sqrt(2))) - p;
    double v = (1.0/sqrt(2*M_PI))*exp(-0.5*x*x);
    x = x - u/(v + x*u/2);
    return x;
}
static double p_to_z(double p){ return inv_norm_cdf(1.0 - 0.5*p); }

// ---------- Row model ----------
struct Row {
    string id, A1, A2;
    double P=numeric_limits<double>::quiet_NaN();
    double N=numeric_limits<double>::quiet_NaN();
    double N_cas=numeric_limits<double>::quiet_NaN();
    double N_con=numeric_limits<double>::quiet_NaN();
    double INFO=numeric_limits<double>::quiet_NaN();
    double FRQ=numeric_limits<double>::quiet_NaN();
    double SSTAT=numeric_limits<double>::quiet_NaN(); // signed stat as provided
    double Z=numeric_limits<double>::quiet_NaN();     // output
    int    z_sign = 1;                                // +1 normally; -1 if we swapped to match --merge-alleles
};

// ---------- Merge-alleles support ----------
// Map SNP -> pair(M1,M2) from the merge file (uppercase)
using MergeMap = unordered_map<string, pair<char,char>>;

static MergeMap read_merge_alleles(const string& path){
    if (path.empty()) return {};
    ifstream in(path);
    if (!in) throw runtime_error("Could not open --merge-alleles file: " + path);
    string line;
    if (!getline(in, line)) throw runtime_error("Empty --merge-alleles file");
    auto hdr = split_ws(line);
    // Expect exactly SNP, A1, A2 (in any order)
    int iS=-1,i1=-1,i2=-1;
    for (int i=0;i<(int)hdr.size();++i){
        string h = hdr[i];
        for (auto& c: h) c = toupper(static_cast<unsigned char>(c));
        if (h=="SNP") iS=i; else if (h=="A1") i1=i; else if (h=="A2") i2=i;
    }
    if (iS<0 || i1<0 || i2<0) throw runtime_error("--merge-alleles must have columns SNP, A1, A2");
    MergeMap mp; mp.reserve(1<<20);
    while (getline(in, line)){
        if (line.empty()) continue;
        auto f = split_ws(line);
        if ((int)f.size() <= max(iS, max(i1,i2))) continue;
        string snp = f[iS];
        string a1 = f[i1], a2 = f[i2];
        if (snp.empty() || a1.empty() || a2.empty()) continue;
        char A = toupper(static_cast<unsigned char>(a1[0]));
        char B = toupper(static_cast<unsigned char>(a2[0]));
        if (!is_valid_base(A) || !is_valid_base(B)) continue;
        mp.emplace(snp, make_pair(A,B));
    }
    return mp;
}

// True if observed alleles (x1,x2) are compatible with reference (r1,r2)
// allowing ref flip and/or strand flip (same as MATCH_ALLELES in LDSC).
static bool alleles_match(char x1, char x2, char r1, char r2){
    x1 = toupper(static_cast<unsigned char>(x1));
    x2 = toupper(static_cast<unsigned char>(x2));
    r1 = toupper(static_cast<unsigned char>(r1));
    r2 = toupper(static_cast<unsigned char>(r2));
    if (!is_valid_base(x1)||!is_valid_base(x2)||!is_valid_base(r1)||!is_valid_base(r2)) return false;

    // strand & ref match
    if (x1==r1 && x2==r2) return true;
    // ref match, strand flip
    if (x1==comp_base(r1) && x2==comp_base(r2)) return true;
    // ref flip, strand match
    if (x1==r2 && x2==r1) return true;
    // ref flip, strand flip
    if (x1==comp_base(r2) && x2==comp_base(r1)) return true;
    return false;
}

// ---------- IO helpers ----------
static vector<string> read_header(ifstream& ifs){
    string line; if(!getline(ifs,line)) throw runtime_error("Empty --sumstats file");
    return split_ws(line);
}

// ---------- Main parse/filter ----------
static vector<Row> parse_and_filter(const Args& a, Logger& log){
    // Load merge-alleles map first (if any)
    MergeMap merge;
    if (!a.merge_alleles.empty()){
        log.log("Reading --merge-alleles from ", a.merge_alleles);
        merge = read_merge_alleles(a.merge_alleles);
        log.log("Read ", merge.size(), " SNPs for allele merge.");
    }

    ifstream ifs(a.sumstats);
    if(!ifs) throw runtime_error("Could not open --sumstats");
    vector<string> hdr = read_header(ifs);
    unordered_map<string,int> idx; for(int i=0;i<(int)hdr.size();++i) idx[hdr[i]] = i;

    auto must_idx = [&](const string& name)->int{
        auto it = idx.find(name);
        if (it==idx.end()) throw runtime_error("Missing required column: "+name);
        return it->second;
    };

    // required column indices
    int id_i = must_idx(a.id_col);
    int p_i  = must_idx(a.P_col);
    int A1_i = -1, A2_i = -1;
    if (!a.no_alleles){ A1_i = must_idx(a.A1_col); A2_i = must_idx(a.A2_col); }

    // optional indices
    int N_i=-1, Nc_i=-1, Nn_i=-1, info_i=-1, maf_i=-1, sstat_i=-1;
    if (!a.N_col.empty()     && idx.count(a.N_col))     N_i = idx[a.N_col];
    if (!a.N_cas_col.empty() && idx.count(a.N_cas_col)) Nc_i = idx[a.N_cas_col];
    if (!a.N_con_col.empty() && idx.count(a.N_con_col)) Nn_i = idx[a.N_con_col];
    if (!a.info_col.empty()  && idx.count(a.info_col))  info_i = idx[a.info_col];
    if (!a.maf_col.empty()   && idx.count(a.maf_col))   maf_i  = idx[a.maf_col];
    sstat_i = must_idx(a.stats_col);

    vector<Row> rows; rows.reserve(1<<20);
    size_t total_read=0, kept=0;
    size_t drop_na=0, drop_p_nonnum=0, drop_info=0, drop_maf=0, drop_alleles=0,
        drop_merge_notin=0, drop_merge_mismatch=0,
        kept_merge_swap=0;

    string line;
    while (getline(ifs,line)){
        if (line.empty()) continue;
        ++total_read;

        vector<string> f = split_ws(line);
        if ((int)f.size() <= max({id_i,p_i,A1_i,A2_i,N_i,Nc_i,Nn_i,info_i,maf_i,sstat_i})) { drop_na++; continue; }

        Row r;
        r.id = f[id_i];
        if (r.id.empty()){ drop_na++; continue; }

        // If merging alleles: restrict to listed SNPs
        if (!merge.empty() && !merge.count(r.id)) { drop_merge_notin++; continue; }

        // P (sanitize instead of dropping out-of-bounds)
        {
            const string& s = f[p_i];
            if (s=="."||s=="NA"){ drop_na++; continue; }
            char* e=nullptr; double pv = strtod(s.c_str(), &e);
            if (!e || *e!='\0'){ drop_p_nonnum++; continue; } // non-numeric only
            r.P = sanitize_p_counting(pv);
        }

        // SSTAT
        {
            const string& s = f[sstat_i];
            if (s!="." && s!="NA"){
                char* e=nullptr; r.SSTAT = strtod(s.c_str(), &e);
                if (!e || *e!='\0') r.SSTAT = numeric_limits<double>::quiet_NaN();
            }
        }

        // Alleles
        if (!a.no_alleles){
            r.A1 = f[A1_i]; r.A2 = f[A2_i];
            if (r.A1.empty()||r.A2.empty()){ drop_na++; continue; }
            for (auto& c: r.A1) c = toupper(static_cast<unsigned char>(c));
            for (auto& c: r.A2) c = toupper(static_cast<unsigned char>(c));
            if (!keep_alleles(r.A1, r.A2)){ drop_alleles++; continue; }
            // If merging, enforce allele compatibility with reference file
            if (!merge.empty()){
                auto it = merge.find(r.id);
                // if absent we would have dropped earlier
                if (it!=merge.end()){
                    char R1 = it->second.first, R2 = it->second.second;
                    if (!alleles_match(r.A1[0], r.A2[0], R1, R2)){
                        drop_merge_mismatch++;
                        continue;
                    }
                    // Strict policy:
                    // - exact match: keep
                    // - reversed: swap alleles and mark Z to flip later
                    // - anything else (incl. strand complements): drop
                    if (r.A1[0]==R1 && r.A2[0]==R2) {
                        // keep
                    } else if (r.A1[0]==R2 && r.A2[0]==R1) {
                        std::swap(r.A1, r.A2);
                        r.z_sign = -1;
                        kept_merge_swap++;
                    } else {
                        drop_merge_mismatch++;
                        continue;
                    }
                }
            }
        }

        // N / cases/controls
        if (N_i>=0){
            const string& s = f[N_i];
            if (s!="." && s!="NA"){ char* e=nullptr; r.N = strtod(s.c_str(), &e); if (!e||*e!='\0') r.N = NAN; }
        }
        if (Nc_i>=0){
            const string& s = f[Nc_i];
            if (s!="." && s!="NA"){ char* e=nullptr; r.N_cas = strtod(s.c_str(), &e); if (!e||*e!='\0') r.N_cas = NAN; }
        }
        if (Nn_i>=0){
            const string& s = f[Nn_i];
            if (s!="." && s!="NA"){ char* e=nullptr; r.N_con = strtod(s.c_str(), &e); if (!e||*e!='\0') r.N_con = NAN; }
        }

        // INFO/MAF
        if (info_i>=0){
            const string& s = f[info_i];
            if (s!="." && s!="NA"){ char* e=nullptr; r.INFO = strtod(s.c_str(), &e); if (!e||*e!='\0') r.INFO = NAN; }
            if (!pass_info(r.INFO, a.info_min)){ drop_info++; continue; }
        }
        if (maf_i>=0){
            const string& s = f[maf_i];
            if (s!="." && s!="NA"){ char* e=nullptr; r.FRQ = strtod(s.c_str(), &e); if (!e||*e!='\0') r.FRQ = NAN; }
            double tmp = r.FRQ;
            if (!pass_maf(tmp, a.maf_min)){ drop_maf++; continue; }
            r.FRQ = tmp; // now MAF
        }

        rows.push_back(std::move(r));
        ++kept;
        if (kept % 1000000 == 0) cerr << ".";
    }
    cerr << " done\n";

    // Summary
    {
        ostringstream msg;
        msg << "Read " << total_read << " SNPs from --sumstats file.\n";
        if (!a.merge_alleles.empty()){
            msg << "Removed " << drop_merge_notin << " SNPs not in --merge-alleles.\n";
        }
        msg << "Removed " << drop_na           << " SNPs with missing values.\n";
        msg << "Removed " << drop_info         << " SNPs with INFO <= " << a.info_min << ".\n";
        msg << "Removed " << drop_maf          << " SNPs with MAF <= "  << a.maf_min  << ".\n";
        msg << "Removed " << drop_p_nonnum     << " SNPs with non-numeric p-values.\n";
        msg << "Removed " << drop_alleles      << " variants that were not SNPs or were strand-ambiguous.\n";
        if (!a.merge_alleles.empty()){
            msg << "Removed " << drop_merge_mismatch << " SNPs whose alleles did not match --merge-alleles.\n";
            msg << "Flipped alleles for " << kept_merge_swap << " SNPs to match --merge-alleles (Z sign reversed).\n";
        }
        msg << rows.size() << " SNPs remain.";
        log.log(msg.str());

        // Report p-value sanitization counts (informational)
        if (g_pstats.zero_or_neg || g_pstats.over_one) {
            log.log("Sanitized p-values: ", (unsigned long long)g_pstats.zero_or_neg,
                    " <= 0; ", (unsigned long long)g_pstats.over_one, " > 1");
        }
    }

    // Dedup by RSID
    {
        size_t before = rows.size();
        unordered_set<string> seen; seen.reserve(rows.size()*2);
        vector<Row> dedup; dedup.reserve(rows.size());
        for (auto& r: rows) if (seen.insert(r.id).second) dedup.push_back(std::move(r));
        rows.swap(dedup);
        log.log("Removed ", (before-rows.size()), " SNPs with duplicated rs numbers (", rows.size(), " SNPs remain).");
    }

    return rows;
}

// ---------- N handling & write out ----------
static void process_and_write(vector<Row>& rows, const Args& a, Logger& log){
    // Fill N if missing
    bool anyN=false;
    for (auto& r: rows){ if (!std::isnan(r.N)) { anyN=true; break; } }
    if (!anyN){
        if (!std::isnan(a.N_flag)){
            for (auto& r: rows) r.N = a.N_flag;
            log.log("Using constant N = ", a.N_flag);
        } else if (!std::isnan(a.N_cas) && !std::isnan(a.N_con)){
            double tot = a.N_cas + a.N_con;
            for (auto& r: rows) r.N = tot;
            log.log("Using N_cas = ", a.N_cas, "; N_con = ", a.N_con);
        } else {
            bool have_row_cc=false;
            for (auto& r: rows) if(!std::isnan(r.N_cas) && !std::isnan(r.N_con)){ have_row_cc=true; break; }
            if (have_row_cc){
                for (auto& r: rows) if (std::isnan(r.N) && !std::isnan(r.N_cas) && !std::isnan(r.N_con)) r.N = r.N_cas + r.N_con;
                log.log("Synthesized N from row-wise cases+controls.");
            } else {
                throw runtime_error("Cannot determine N (no --N/--N-cas/--N-con and no N/N_cases/N_cont columns).");
            }
        }
    }

    // Drop low N: default (90th percentile)/1.5
    vector<double> Ns; Ns.reserve(rows.size());
    for (auto& r: rows) if (!std::isnan(r.N)) Ns.push_back(r.N);
    if (!Ns.empty()){
        nth_element(Ns.begin(), Ns.begin()+(Ns.size()*9)/10, Ns.end());
        double n90 = Ns[(Ns.size()*9)/10];
        double nmin = n90/1.5;
        size_t before = rows.size();
        rows.erase(remove_if(rows.begin(), rows.end(), [&](const Row& r){ return std::isnan(r.N) || r.N < nmin; }), rows.end());
        log.log("Removed ", (before-rows.size()), " SNPs with N < ", nmin, " (", rows.size(), " SNPs remain).");
    }

    // Compute Z from sanitized P
    for (auto& r: rows) r.Z = p_to_z(r.P);

    // Orient Z by signed stat
    double null_signed = (a.stats_type=="or") ? 1.0 : 0.0;
    for (auto& r: rows){
        if (!std::isnan(r.SSTAT) && !std::isnan(r.Z)){
            if (r.SSTAT < null_signed) r.Z = -r.Z;
        }
        // If we swapped A1/A2 earlier to align with the reference, reverse Z sign.
        if (r.z_sign < 0 && !std::isnan(r.Z)) {
            r.Z = -r.Z;
        }
    }

    // Write output
    const string ofn = a.out + ".sumstats";
    ofstream ofs(ofn);
    if (!ofs) throw runtime_error("Could not open output: "+ofn);
    ofs.setf(std::ios::fixed); ofs<<setprecision(3);
    ofs << "SNP\tN\tZ";
    if (!a.no_alleles) ofs << "\tA1\tA2";
    if (a.keep_maf && !a.maf_col.empty()) ofs << "\tFRQ";
    ofs << "\n";
    size_t printed=0;
    for (auto& r: rows){
        ofs << r.id << '\t' << r.N << '\t' << r.Z;
        if (!a.no_alleles) ofs << '\t' << r.A1 << '\t' << r.A2;
        if (a.keep_maf && !a.maf_col.empty()) ofs << '\t' << (std::isnan(r.FRQ)?0.0:r.FRQ);
        ofs << '\n';
        printed++;
    }
    ofs.close();
    log.log("Writing summary statistics for ", printed, " SNPs to ", ofn, ".");

    // QC metrics (finite-only)
    vector<double> chisq; chisq.reserve(rows.size());
    for (auto& r: rows){
        if (std::isfinite(r.Z)) {
            double c = r.Z * r.Z;
            if (std::isfinite(c)) chisq.push_back(c);
        }
    }

    if (!chisq.empty()){
        double sum = 0.0; for (double c: chisq) sum += c;
        double mean_chi2 = sum / chisq.size();
        log.log("Mean chi^2 = ", std::round(mean_chi2*1000)/1000.0);

        vector<double> tmp = chisq;
        nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end());
        double med = tmp[tmp.size()/2];
        double lambda_gc = med / 0.4549; // median of chi2_1
        ostringstream oss; oss<<fixed<<setprecision(3)<<lambda_gc;
        log.log("Lambda GC = ", oss.str());

        double mx = *max_element(chisq.begin(), chisq.end());
        ostringstream oss2; oss2<<fixed<<setprecision(3)<<mx;
        log.log("Max chi^2 = ", oss2.str());

        // p < 5e-8 threshold for 1 df chi-square ≈ 29.716
        const double GW_CHI2 = 29.716;
        size_t gws = count_if(chisq.begin(), chisq.end(), [&](double x){ return x > GW_CHI2; });
        log.log(to_string(gws), " Genome-wide significant SNPs (some may have been removed by filtering).");
    } else {
        log.log("No finite chi^2 values to summarize.");
    }
}

int run_munge(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    try{
        Args args = parse_args(argc, argv);
        Logger log(args.out + ".log");

        // Echo call (compact)
        {
            ostringstream call;
            call << "Call: --sumstats " << args.sumstats << " --out " << args.out
                 << " --rsid " << args.id_col << " --p-value " << args.P_col
                 << " --signed-col " << args.stats_col << " --signed-type " << args.stats_type;
            if (!args.no_alleles) call << " --A1 " << args.A1_col << " --A2 " << args.A2_col;
            if (!args.N_col.empty()) call << " --N-col " << args.N_col;
            if (!args.N_cas_col.empty()) call << " --N_cases " << args.N_cas_col;
            if (!args.N_con_col.empty()) call << " --N_cont " << args.N_con_col;
            if (!args.info_col.empty()) call << " --info " << args.info_col << " --info-min " << args.info_min;
            if (!args.maf_col.empty())  call << " --maf "  << args.maf_col  << " --maf-min "  << args.maf_min;
            if (!args.merge_alleles.empty()) call << " --merge-alleles " << args.merge_alleles;
            if (args.keep_maf) call << " --keep-maf";
            if (!std::isnan(args.N_flag)) call << " --N " << args.N_flag;
            if (!std::isnan(args.N_cas))  call << " --N-cas " << args.N_cas;
            if (!std::isnan(args.N_con))  call << " --N-con " << args.N_con;
            log.log(call.str());
        }

        auto rows = parse_and_filter(args, log);
        if (rows.empty()) throw runtime_error("After applying filters, no SNPs remain.");
        process_and_write(rows, args, log);
        log.log("Conversion finished.");
        return 0;
    } catch (const exception& e){
        cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
#undef Logger

