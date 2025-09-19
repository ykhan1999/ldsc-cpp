// ldsc.cpp — LDSC (single-annotation) with proper two-step option
// - Univariate h2 with optional two-step (Step1 intercept on filtered SNPs; Step2 slope on all SNPs with intercept fixed)
// - Bivariate rg uses univariate fits (so benefits from two-step), bivariate intercept estimated as before
// - No partitioned outputs
// (C) 2025  GPL-3.0-or-later

#include "logger.hpp"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

int run_ldsc(int argc, char** argv);

// ------------------------- CLI -------------------------
struct LdscArgs {
    std::string out;        // --out
    std::string h2;         // --h2 <sumstats>
    std::string rg1, rg2;   // --rg f1,f2
    std::string ref_ld;     // --ref-ld
    std::string ref_ld_chr; // --ref-ld-chr
    std::string w_ld;       // --w-ld
    std::string w_ld_chr;   // --w-ld-chr
    std::string M;          // --M
    int n_blocks = 200;     // --n-blocks
    bool fix_int = false;   // --no-intercept (univariate a=1 fixed)
    bool fix_int12 = false; // --no-intercept12 (bivariate a12=0 fixed)
    int irwls_iter = 2;     // --irwls-iter
    bool two_step = false;  // --two-step [thresh]
    double two_step_thresh = 30.0;
};

static LdscArgs parse_ldsc(int argc, char** argv) {
    auto need = [&](int& i, const std::string& k){
        if (i+1 >= argc) { std::cerr<<"Missing value for "<<k<<"\n"; std::exit(2); }
        return std::string(argv[++i]);
    };

    LdscArgs a;
    for (int i=0; i<argc; ++i) {
        std::string k = argv[i];
        if      (k=="--out")            a.out = need(i,k);
        else if (k=="--h2")             a.h2 = need(i,k);
        else if (k=="--rg")             { auto s=need(i,k); auto p=s.find(','); if(p==std::string::npos){std::cerr<<"--rg needs f1,f2\n"; std::exit(2);} a.rg1=s.substr(0,p); a.rg2=s.substr(p+1); }
        else if (k=="--ref-ld")         a.ref_ld = need(i,k);
        else if (k=="--ref-ld-chr")     a.ref_ld_chr = need(i,k);
        else if (k=="--w-ld")           a.w_ld = need(i,k);
        else if (k=="--w-ld-chr")       a.w_ld_chr = need(i,k);
        else if (k=="--M")              a.M = need(i,k);
        else if (k=="--n-blocks")       a.n_blocks = std::stoi(need(i,k));
        else if (k=="--no-intercept")   a.fix_int = true;
        else if (k=="--no-intercept12") a.fix_int12 = true;
        else if (k=="--irwls-iter")     a.irwls_iter = std::max(1, std::stoi(need(i,k)));
        else if (k=="--two-step") {
            a.two_step = true;
            // Optional numeric value; default 30 if next token is another flag or end.
            if (i+1 < argc) {
                std::string nxt = argv[i+1];
                if (nxt.rfind("--", 0) != 0) { ++i; a.two_step_thresh = std::stod(nxt); }
            }
        }
        else if (k.rfind("--",0)==0) { std::cerr<<"Unknown flag in ldsc: "<<k<<"\n"; std::exit(2); }
    }

    if (a.out.empty()) { std::cerr<<"ldsc: --out required\n"; std::exit(2); }
    const bool want_h2 = !a.h2.empty();
    const bool want_rg = !a.rg1.empty() && !a.rg2.empty();
    if (want_h2 == want_rg) { std::cerr<<"ldsc: specify exactly one of --h2 or --rg\n"; std::exit(2); }
    if ((a.ref_ld.empty() == a.ref_ld_chr.empty())) {
        std::cerr<<"ldsc: specify exactly one of --ref-ld or --ref-ld-chr\n"; std::exit(2);
    }
    if ((a.w_ld.empty() == a.w_ld_chr.empty())) {
        std::cerr<<"ldsc: specify exactly one of --w-ld or --w-ld-chr\n"; std::exit(2);
    }
    if (a.n_blocks <= 1) { std::cerr<<"ldsc: --n-blocks must be > 1\n"; std::exit(2); }
    return a;
}

static void echo_call(Logger& log, const LdscArgs& a) {
    std::ostringstream ss;
    ss << "Call (ldsc): --out " << a.out;
    if (!a.h2.empty()) ss << " --h2 " << a.h2;
    if (!a.rg1.empty()) ss << " --rg " << a.rg1 << "," << a.rg2;
    ss << (a.ref_ld.empty() ? " --ref-ld-chr "+a.ref_ld_chr : " --ref-ld "+a.ref_ld);
    ss << (a.w_ld.empty()   ? " --w-ld-chr "+a.w_ld_chr     : " --w-ld "+a.w_ld);
    if (!a.M.empty()) ss << " --M " << a.M;
    if (a.fix_int) ss << " --no-intercept";
    if (a.fix_int12) ss << " --no-intercept12";
    if (a.n_blocks!=200) ss << " --n-blocks " << a.n_blocks;
    if (a.irwls_iter!=2) ss << " --irwls-iter " << a.irwls_iter;
    if (a.two_step) ss << " --two-step " << a.two_step_thresh;
    log.log(ss.str());
}

// --------------------- Small utilities ---------------------
static std::vector<std::string> split_ws(const std::string& s){
    std::vector<std::string> v; std::string t; std::istringstream iss(s);
    while (iss>>t) v.push_back(t);
    return v;
}

static int find_col(const std::vector<std::string>& hdr, const std::string& name){
    for (int i=0;i<(int)hdr.size();++i) if (hdr[i]==name) return i;
    return -1;
}

static std::vector<std::string> expand_chr_input(const std::string& s, const std::string& suffix) {
    std::vector<std::string> out;
    auto pos = s.find("{}");
    if (pos != std::string::npos) { out.reserve(22); for (int c=1;c<=22;++c){ std::ostringstream one; one<<s.substr(0,pos)<<c<<s.substr(pos+2); out.push_back(one.str()); } return out; }
    if (std::filesystem::is_directory(s)) { out.reserve(22); for (int c=1;c<=22;++c){ std::ostringstream one; one<<s<<"/"<<c<<suffix; out.push_back(one.str()); } return out; }
    out.push_back(s); return out;
}

struct Table {
    std::vector<std::string> snp;
    std::vector<std::string> cols;
    std::vector<std::vector<double>> data; // data[j][i]
};
struct Sumstats {
    std::vector<std::string> snp;
    std::vector<double> Z, N;
};

static Sumstats read_sumstats(const std::string& path, Logger& log) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open sumstats: "+path);
    std::string line; if (!std::getline(in,line)) throw std::runtime_error("Empty sumstats: "+path);
    auto hdr = split_ws(line);
    int iS=-1, iZ=-1, iN=-1;
    for (int i=0;i<(int)hdr.size();++i){ if (hdr[i]=="SNP") iS=i; else if (hdr[i]=="Z") iZ=i; else if (hdr[i]=="N") iN=i; }
    if (iS<0||iZ<0||iN<0) throw std::runtime_error("sumstats must have SNP Z N columns");
    Sumstats ss;
    size_t total=0, drop=0;
    std::string row;
    while (std::getline(in,row)){
        if (row.empty()) continue;
        ++total;
        auto f = split_ws(row);
        if ((int)f.size() <= std::max(iS,std::max(iZ,iN))) { ++drop; continue; }
        const std::string& sS=f[iS], &sZ=f[iZ], &sN=f[iN];
        if (sS.empty()) { ++drop; continue; }
        char* e1=nullptr; double z = std::strtod(sZ.c_str(), &e1);
        char* e2=nullptr; double n = std::strtod(sN.c_str(), &e2);
        if (!e1||*e1!='\0'||!std::isfinite(z)) { ++drop; continue; }
        if (!e2||*e2!='\0'||!std::isfinite(n)||(n<=0)) { ++drop; continue; }
        ss.snp.push_back(sS); ss.Z.push_back(z); ss.N.push_back(n);
    }
    log.log("Read ", (int)ss.snp.size(), " SNPs from ", path, " (parsed ", total, ", dropped ", drop, ").");
    return ss;
}

static std::vector<int> pick_ref_1col(const std::vector<std::string>& hdr){
    const char* prefer[] = {"L2","LDSCORE","L2_no_mhc"};
    for (const char* p : prefer)
        for (int i=0;i<(int)hdr.size();++i) if (hdr[i]==p) return {i};
    static const std::unordered_set<std::string> meta = {"SNP","CHR","BP","CM","MAF","A1","A2","L2_MAF"};
    for (int i=0;i<(int)hdr.size();++i) if (!meta.count(hdr[i]) && hdr[i]!="SNP") return {i};
    throw std::runtime_error("ref-LD has no usable numeric column.");
}
static std::vector<int> pick_w_1col(const std::vector<std::string>& hdr){ return pick_ref_1col(hdr); }

static Table read_ld_one_select(const std::string& path, const std::vector<int>& idx) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open LD file: "+path);
    std::string line; if (!std::getline(in,line)) throw std::runtime_error("Empty LD file: "+path);
    auto hdr = split_ws(line);
    int iSNP = find_col(hdr, "SNP");
    if (iSNP<0) throw std::runtime_error("LD file missing SNP header: "+path);

    Table t; t.cols.reserve(idx.size()); t.data.resize(idx.size());
    for (int j: idx) t.cols.push_back(hdr[j]);

    std::string row;
    while (std::getline(in,row)){
        if (row.empty()) continue;
        auto f = split_ws(row);
        if ((int)f.size() <= iSNP) continue;
        t.snp.push_back(f[iSNP]);
        for (size_t k=0;k<idx.size();++k){
            int j = idx[k];
            double v = (j<(int)f.size()) ? std::strtod(f[j].c_str(), nullptr) : std::numeric_limits<double>::quiet_NaN();
            t.data[k].push_back(v);
        }
    }
    return t;
}

static Table read_ld_fileset_1col(const std::string& single_or_dir_or_pattern, const std::string& suffix, bool for_weights){
    Table agg; bool first=true;
    auto files = expand_chr_input(single_or_dir_or_pattern, suffix);
    for (const auto& f: files) {
        std::ifstream in(f); if(!in) throw std::runtime_error("Could not open LD file: "+f);
        std::string line; if(!std::getline(in,line)) throw std::runtime_error("Empty LD file: "+f);
        auto hdr = split_ws(line);
        auto idx = for_weights? pick_w_1col(hdr) : pick_ref_1col(hdr);
        in.close();
        auto t = read_ld_one_select(f, idx);
        if (first) { agg = std::move(t); first=false; }
        else {
            if (t.cols != agg.cols) throw std::runtime_error("LD column sets differ across files: "+f);
            agg.snp.insert(agg.snp.end(), t.snp.begin(), t.snp.end());
            for (size_t j=0;j<agg.data.size();++j)
                agg.data[j].insert(agg.data[j].end(), t.data[j].begin(), t.data[j].end());
        }
    }
    if (agg.cols.size()!=1) throw std::runtime_error("Expected exactly one LD column.");
    return agg;
}

static long long detect_M_total(const std::string& ref_ld_chr, size_t ld_rows){
    if (std::filesystem::is_directory(ref_ld_chr)) {
        long long total=0;
        for (int c=1;c<=22;++c){
            std::ostringstream p; p<<ref_ld_chr<<"/"<<c<<".l2.M_5_50";
            std::ifstream in(p.str()); if (!in) return -1;
            long long m_chr=0; if (!(in>>m_chr)) return -1;
            total += m_chr;
        }
        return total;
    }
    return (long long)ld_rows;
}

struct Merge1 {
    std::vector<std::string> snp;
    std::vector<double> y_chi2;  // Z^2
    std::vector<double> N;
    std::vector<double> ld;      // ref LD
    std::vector<double> wld;     // weight LD
};
struct Merge2 {
    std::vector<std::string> snp;
    std::vector<double> Z1, Z2;
    std::vector<double> N1, N2;
    std::vector<double> ld, wld;
};

static Merge1 merge_univ(const Sumstats& ss, const Table& refld, const Table& wld) {
    std::unordered_map<std::string,int> ir, iw;
    ir.reserve(refld.snp.size()*2); iw.reserve(wld.snp.size()*2);
    for (int i=0;i<(int)refld.snp.size();++i) ir.emplace(refld.snp[i], i);
    for (int i=0;i<(int)wld.snp.size();++i)   iw.emplace(wld.snp[i], i);

    Merge1 m; size_t drop=0;
    for (size_t i=0;i<ss.snp.size();++i){
        auto itr = ir.find(ss.snp[i]);
        auto itw = iw.find(ss.snp[i]);
        if (itr==ir.end() || itw==iw.end()) { ++drop; continue; }
        double z=ss.Z[i], n=ss.N[i];
        double l = refld.data[0][itr->second];
        double wl = wld.data[0][itw->second];
        if (!std::isfinite(z)||!std::isfinite(n)||(n<=0)||!std::isfinite(l)||!std::isfinite(wl)) { ++drop; continue; }
        m.snp.push_back(ss.snp[i]);
        m.y_chi2.push_back(z*z);
        m.N.push_back(n);
        m.ld.push_back(l);
        m.wld.push_back(wl);
    }
    if (m.snp.empty()) throw std::runtime_error("merge_univ: no overlapping SNPs.");
    return m;
}

static Merge2 merge_bivar(const Sumstats& s1, const Sumstats& s2, const Table& refld, const Table& wld){
    std::unordered_map<std::string,int> i1,i2,ir,iw;
    i1.reserve(s1.snp.size()*2); i2.reserve(s2.snp.size()*2);
    ir.reserve(refld.snp.size()*2); iw.reserve(wld.snp.size()*2);
    for (int i=0;i<(int)s1.snp.size();++i) i1.emplace(s1.snp[i], i);
    for (int i=0;i<(int)s2.snp.size();++i) i2.emplace(s2.snp[i], i);
    for (int i=0;i<(int)refld.snp.size();++i) ir.emplace(refld.snp[i], i);
    for (int i=0;i<(int)wld.snp.size();++i)   iw.emplace(wld.snp[i], i);

    Merge2 m;
    for (const auto& kv : i1){
        const auto& snp = kv.first;
        auto a = i2.find(snp); auto b = ir.find(snp); auto c = iw.find(snp);
        if (a==i2.end()||b==ir.end()||c==iw.end()) continue;
        int i = kv.second, j=a->second, r=b->second, w=c->second;
        double z1=s1.Z[i], n1=s1.N[i], z2=s2.Z[j], n2=s2.N[j];
        double l = refld.data[0][r], wl = wld.data[0][w];
        if (!std::isfinite(z1)||!std::isfinite(z2)||!std::isfinite(n1)||!std::isfinite(n2)
            || !(n1>0) || !(n2>0) || !std::isfinite(l) || !std::isfinite(wl)) continue;
        m.snp.push_back(snp);
        m.Z1.push_back(z1); m.Z2.push_back(z2);
        m.N1.push_back(n1); m.N2.push_back(n2);
        m.ld.push_back(l);  m.wld.push_back(wl);
    }
    if (m.snp.empty()) throw std::runtime_error("merge_bivar: no overlapping SNPs.");
    return m;
}

static std::vector<int> block_boundaries(int n, int B){
    B = std::max(2, std::min(B, n));
    std::vector<int> cuts(B+1,0);
    for (int b=0;b<=B;++b) cuts[b] = (int)std::llround((long double)b * n / B);
    return cuts;
}
static std::vector<int> all_idx(int n){ std::vector<int> v(n); std::iota(v.begin(), v.end(), 0); return v; }

// X'WX utilities
static void xtwx_xtwy(const std::vector<std::vector<double>>& X, const std::vector<double>& w,
                      const std::vector<double>& y, std::vector<std::vector<double>>& XtWX,
                      std::vector<double>& XtWy) {
    int p=(int)X.size(), n=(int)w.size();
    XtWX.assign(p, std::vector<double>(p,0.0));
    XtWy.assign(p,0.0);
    for (int i=0;i<n;++i){
        double wi=w[i];
        for (int a=0;a<p;++a){
            double xa=X[a][i];
            XtWy[a]+=xa*wi*y[i];
            for (int b=a;b<p;++b) XtWX[a][b]+=xa*wi*X[b][i];
        }
    }
    for (int a=0;a<p;++a) for (int b=0;b<a;++b) XtWX[a][b]=XtWX[b][a];
}
static std::vector<double> solve_spd(std::vector<std::vector<double>> A, std::vector<double> b){
    int p=(int)A.size();
    for (int i=0;i<p;++i) A[i][i]+=1e-12;
    for (int i=0;i<p;++i){
        int piv=i; double best=std::fabs(A[i][i]);
        for (int r=i+1;r<p;++r){ double v=std::fabs(A[r][i]); if (v>best){best=v; piv=r;} }
        if (best==0) throw std::runtime_error("Singular normal matrix.");
        if (piv!=i){ std::swap(A[piv],A[i]); std::swap(b[piv],b[i]); }
        double d=A[i][i];
        for (int j=i;j<p;++j) A[i][j]/=d; b[i]/=d;
        for (int r=0;r<p;++r){
            if (r==i) continue;
            double f=A[r][i]; if (f==0) continue;
            for (int j=i;j<p;++j) A[r][j]-=f*A[i][j];
            b[r]-=f*b[i];
        }
    }
    return b;
}
static std::vector<double> wls(const std::vector<std::vector<double>>& X,
                               const std::vector<double>& w,
                               const std::vector<double>& y){
    std::vector<std::vector<double>> XtWX; std::vector<double> XtWy;
    xtwx_xtwy(X,w,y,XtWX,XtWy);
    return solve_spd(XtWX, XtWy);
}

// -------------- Weights (subset-aware) --------------
static void weights_univ_indices(const Merge1& m, const std::vector<double>& x_slope,
                                 double a, double h2,
                                 const std::vector<int>& idx, std::vector<double>& w){
    int k=(int)idx.size(); w.resize(k);
    for (int t=0;t<k;++t){
        int i=idx[t];
        double mu = a + h2 * x_slope[i];
        if (mu < 1.0) mu = 1.0;
        double wl = m.wld[i]; if (!(wl>0)) wl = 1.0;
        double var = 2.0 * mu * mu;
        double inv = 1.0 / (var * wl);
        if (!std::isfinite(inv) || inv<=0) inv = 1.0;
        w[t] = inv;
    }
}
static void weights_univ_fixed_indices(const Merge1& m, const std::vector<double>& x_slope,
                                       double a_fixed, double h2,
                                       const std::vector<int>& idx, std::vector<double>& w){
    weights_univ_indices(m, x_slope, a_fixed, h2, idx, w);
}

static void weights_bivar_indices(const Merge2& m,
                          const std::vector<double>& x1, const std::vector<double>& x2,
                          double a1, double h21,
                          double a2, double h22,
                          const std::vector<double>& x12,
                          double a12, double rho,
                          const std::vector<int>& idx,
                          std::vector<double>& w){
    int k=(int)idx.size(); w.resize(k);
    for (int t=0;t<k;++t){
        int i=idx[t];
        double mu1 = a1 + h21 * x1[i]; if (mu1 < 1.0) mu1 = 1.0;
        double mu2 = a2 + h22 * x2[i]; if (mu2 < 1.0) mu2 = 1.0;
        double mu12= a12 + rho * x12[i];
        double wl = m.wld[i]; if (!(wl>0)) wl=1.0;
        double var = mu1*mu2 + mu12*mu12;
        double inv = 1.0 / (var * wl);
        if (!std::isfinite(inv) || inv<=0) inv = 1.0;
        w[t] = inv;
    }
}

// -------------- Univariate LDSC --------------
struct H2Fit { double a, a_se, h2, h2_se, Nbar, Mtot; int step1_kept = -1; };

static void jackknife_fit_univ_free(const Merge1& m, double Mtot, int B, int iters,
                                    H2Fit& out){
    int n=(int)m.snp.size();
    auto cuts = block_boundaries(n, B);

    // Precompute slope regressor
    std::vector<double> xs(n);
    for (int i=0;i<n;++i) xs[i] = (m.N[i]/Mtot) * m.ld[i];

    // Moments init
    double Nbar = std::accumulate(m.N.begin(), m.N.end(), 0.0) / n;
    double mean_chi=std::accumulate(m.y_chi2.begin(), m.y_chi2.end(), 0.0) / n;
    double mean_x = std::accumulate(xs.begin(), xs.end(), 0.0) / n;
    double a=1.0, h2=std::max(0.0, (mean_chi - a)/std::max(1e-12, mean_x));

    // Full-data IRWLS (free intercept)
    std::vector<int> idx_all = all_idx(n);
    std::vector<double> w;
    for (int t=0;t<iters;++t){
        weights_univ_indices(m, xs, a, h2, idx_all, w);
        std::vector<std::vector<double>> X(2, std::vector<double>(n,1.0));
        X[1] = xs;
        auto beta = wls(X, w, m.y_chi2);
        a = beta[0]; h2 = beta[1];
    }

    // Delete-one-block jackknife
    std::vector<double> a_del, h2_del; a_del.reserve(B); h2_del.reserve(B);
    for (int b=0;b<B;++b){
        int lo=cuts[b], hi=cuts[b+1];
        std::vector<int> idx; idx.reserve(n-(hi-lo));
        for (int i=0;i<n;++i) if (i<lo || i>=hi) idx.push_back(i);

        double ab=a, h2b=h2;
        for (int t=0;t<iters;++t){
            weights_univ_indices(m, xs, ab, h2b, idx, w);
            std::vector<std::vector<double>> X(2);
            X[0].resize(idx.size(), 1.0);
            X[1].resize(idx.size());
            std::vector<double> y(idx.size());
            for (size_t k=0;k<idx.size();++k){ int i=idx[k]; X[1][k]=xs[i]; y[k]=m.y_chi2[i]; }
            auto bb = wls(X, w, y);
            ab = bb[0]; h2b = bb[1];
        }
        a_del.push_back(ab); h2_del.push_back(h2b);
    }

    auto se_from_delete = [&](const std::vector<double>& del){
        long double mean=0.0L; for (double v: del) mean+=v; mean/= (long double)del.size();
        long double var=0.0L; for (double v: del){ long double d=v-mean; var+=d*d; }
        var *= (del.size()-1.0L) / del.size();
        return std::sqrt((double)var);
    };

    out.a = a; out.h2 = h2; out.a_se = se_from_delete(a_del); out.h2_se = se_from_delete(h2_del);
    out.Nbar = Nbar; out.Mtot = Mtot; out.step1_kept = -1;
}

static void jackknife_fit_univ_twostep(const Merge1& m, double Mtot, int B, int iters,
                                       bool fix_intercept, double step1_thresh,
                                       H2Fit& out){
    int n=(int)m.snp.size();
    auto cuts = block_boundaries(n, B);

    // Precompute slope regressor
    std::vector<double> xs(n);
    for (int i=0;i<n;++i) xs[i] = (m.N[i]/Mtot) * m.ld[i];

    // Helper: SE from delete values
    auto se_from_delete = [&](const std::vector<double>& del){
        long double mean=0.0L; for (double v: del) mean+=v; mean/= (long double)del.size();
        long double var=0.0L; for (double v: del){ long double d=v-mean; var+=d*d; }
        var *= (del.size()-1.0L) / del.size();
        return std::sqrt((double)var);
    };

    // Step 1 (optional): Estimate intercept on filtered SNPs with free intercept
    double a_hat = 1.0;
    std::vector<double> a_del; a_del.reserve(B);
    int step1_kept = -1;

    if (!fix_intercept){
        // Filter: chi^2 <= threshold
        std::vector<int> idx1; idx1.reserve(n);
        for (int i=0;i<n;++i) if (m.y_chi2[i] <= step1_thresh && std::isfinite(xs[i])) idx1.push_back(i);

        // Fallback if too few SNPs — use all SNPs for Step1
        if ((int)idx1.size() < std::max(100, B+1)) { idx1 = all_idx(n); }

        step1_kept = (int)idx1.size();

        // Moments init on Step1 set
        double mean_y=0.0, mean_x=0.0;
        for (int i: idx1){ mean_y += m.y_chi2[i]; mean_x += xs[i]; }
        mean_y /= idx1.size(); mean_x /= idx1.size();
        double a=1.0, h2=std::max(0.0, (mean_y - a)/std::max(1e-12, mean_x));

        std::vector<double> w;
        for (int t=0;t<iters;++t){
            weights_univ_indices(m, xs, a, h2, idx1, w);
            std::vector<std::vector<double>> X(2);
            X[0].resize(idx1.size(), 1.0);
            X[1].resize(idx1.size());
            std::vector<double> y(idx1.size());
            for (size_t k=0;k<idx1.size();++k){ int i=idx1[k]; X[1][k]=xs[i]; y[k]=m.y_chi2[i]; }
            auto beta = wls(X, w, y);
            a = beta[0]; h2 = beta[1];
        }
        a_hat = a;

        // Delete-one-block intercepts on Step1 set (drop the block's indices from idx1)
        for (int b=0;b<B;++b){
            int lo=cuts[b], hi=cuts[b+1];
            std::vector<int> idxb; idxb.reserve(idx1.size());
            for (int i: idx1) if (i<lo || i>=hi) idxb.push_back(i);
            if ((int)idxb.size() < 2) { a_del.push_back(a_hat); continue; } // fallback
            double ab=a, h2b=h2;
            std::vector<double> wdel;
            for (int t=0;t<iters;++t){
                weights_univ_indices(m, xs, ab, h2b, idxb, wdel);
                std::vector<std::vector<double>> Xb(2);
                Xb[0].resize(idxb.size(), 1.0);
                Xb[1].resize(idxb.size());
                std::vector<double> yb(idxb.size());
                for (size_t k=0;k<idxb.size();++k){ int i=idxb[k]; Xb[1][k]=xs[i]; yb[k]=m.y_chi2[i]; }
                auto bb = wls(Xb, wdel, yb);
                ab = bb[0]; h2b = bb[1];
            }
            a_del.push_back(ab);
        }
    } else {
        // Fixed intercept path (a=1 for --no-intercept); we still produce delete-values = 1 for consistency
        a_hat = 1.0;
        a_del.assign(B, 1.0);
        step1_kept = -1;
    }

    // Step 2: slope only on ALL SNPs, intercept fixed to a_hat
    double Nbar = std::accumulate(m.N.begin(), m.N.end(), 0.0) / n;
    std::vector<int> idx_all = all_idx(n);
    std::vector<double> y_adj(n);
    for (int i=0;i<n;++i) y_adj[i] = m.y_chi2[i] - a_hat;

    // Init slope from moments with fixed intercept
    double mean_yadj = std::accumulate(y_adj.begin(), y_adj.end(), 0.0) / n;
    double mean_x    = std::accumulate(xs.begin(),   xs.end(),   0.0) / n;
    double h2 = std::max(0.0, mean_yadj / std::max(1e-12, mean_x));

    std::vector<double> w;
    for (int t=0;t<iters;++t){
        weights_univ_fixed_indices(m, xs, a_hat, h2, idx_all, w);
        std::vector<std::vector<double>> X(1); X[0] = xs;
        auto beta = wls(X, w, y_adj);
        h2 = beta[0];
    }

    // Delete-one-block for slope using a_del[b] per block
    std::vector<double> h2_del; h2_del.reserve(B);
    for (int b=0;b<B;++b){
        int lo=cuts[b], hi=cuts[b+1];
        std::vector<int> idxb; idxb.reserve(n-(hi-lo));
        for (int i=0;i<n;++i) if (i<lo || i>=hi) idxb.push_back(i);

        double ab = a_del[b];
        std::vector<double> yb; yb.resize(idxb.size());
        std::vector<double> xb; xb.resize(idxb.size());
        for (size_t k=0;k<idxb.size();++k){ int i=idxb[k]; yb[k]=m.y_chi2[i]-ab; xb[k]=xs[i]; }

        double h2b = h2;
        for (int t=0;t<iters;++t){
            weights_univ_fixed_indices(m, xs, ab, h2b, idxb, w);
            std::vector<std::vector<double>> Xb(1);
            Xb[0] = xb;
            auto bb = wls(Xb, w, yb);
            h2b = bb[0];
        }
        h2_del.push_back(h2b);
    }

    out.a = a_hat;
    out.a_se = se_from_delete(a_del);
    out.h2 = h2;
    out.h2_se = se_from_delete(h2_del);
    out.Nbar = Nbar; out.Mtot = Mtot; out.step1_kept = step1_kept;
}

// -------------- Bivariate LDSC --------------
struct RGFit { double a12, a12_se, rho, rho_se, rg, rg_se; double Nbar1, Nbar2; double h21, h22; };

static void jackknife_fit_bivar(const Merge2& m, double Mtot,
                                const H2Fit& f1, const H2Fit& f2,
                                int B, int iters, bool fix_int12,
                                RGFit& out){
    int n=(int)m.snp.size();
    auto cuts = block_boundaries(n, B);

    // Precompute regressors
    std::vector<double> x1(n), x2(n), x12(n);
    for (int i=0;i<n;++i){
        x1[i]  = (m.N1[i]/Mtot) * m.ld[i];
        x2[i]  = (m.N2[i]/Mtot) * m.ld[i];
        x12[i] = (std::sqrt(m.N1[i]*m.N2[i]) / Mtot) * m.ld[i];
    }
    // Response: z1*z2
    std::vector<double> y(n); for (int i=0;i<n;++i) y[i] = m.Z1[i]*m.Z2[i];

    // Design
    std::vector<std::vector<double>> X;
    if (!fix_int12){ X.assign(2, std::vector<double>(n,1.0)); X[1]=x12; }
    else            { X.assign(1, x12); }

    // Initial guesses
    double a12 = fix_int12 ? 0.0 : 0.0;
    double mean_y=std::accumulate(y.begin(), y.end(), 0.0)/n;
    double mean_x12=std::accumulate(x12.begin(), x12.end(), 0.0)/n;
    double rho = std::max(0.0, (mean_y - a12) / std::max(1e-12, mean_x12));

    // IRWLS with subset-aware weights
    std::vector<int> idx_all = all_idx(n);
    std::vector<double> w;
    for (int t=0;t<iters;++t){
        weights_bivar_indices(m, x1,x2, f1.a, f1.h2, f2.a, f2.h2, x12, a12, rho, idx_all, w);
        auto beta = wls(X, w, y);
        if (!fix_int12){ a12 = beta[0]; rho = beta[1]; }
        else            { rho = beta[0]; }
    }

    // Jackknife
    std::vector<double> a12_del, rho_del;
    for (int b=0;b<B;++b){
        int lo=cuts[b], hi=cuts[b+1];
        std::vector<int> idxb; idxb.reserve(n-(hi-lo));
        for (int i=0;i<n;++i) if (i<lo || i>=hi) idxb.push_back(i);

        std::vector<std::vector<double>> Xb(X.size());
        for (size_t c=0;c<X.size();++c){
            Xb[c].reserve(idxb.size());
            for (int j: idxb) Xb[c].push_back(X[c][j]);
        }
        std::vector<double> yb; yb.reserve(idxb.size());
        for (int j: idxb) yb.push_back(y[j]);

        double a12b=a12, rhob=rho;
        for (int t=0;t<iters;++t){
            weights_bivar_indices(m, x1,x2, f1.a, f1.h2, f2.a, f2.h2, x12, a12b, rhob, idxb, w);
            auto bb = wls(Xb, w, yb);
            if (!fix_int12){ a12b = bb[0]; rhob = bb[1]; }
            else            { rhob = bb[0]; }
        }
        if (!fix_int12) a12_del.push_back(a12b);
        rho_del.push_back(rhob);
    }

    auto se_from_delete = [&](const std::vector<double>& del){
        long double mean=0.0L; for (double v: del) mean+=v; mean/= (long double)del.size();
        long double var=0.0L; for (double v: del){ long double d=v-mean; var+=d*d; }
        var *= (del.size()-1.0L) / del.size();
        return std::sqrt((double)var);
    };

    out.a12 = fix_int12 ? 0.0 : a12;
    out.rho = rho;
    out.a12_se = fix_int12 ? 0.0 : se_from_delete(a12_del);
    out.rho_se = se_from_delete(rho_del);
    out.h21 = f1.h2; out.h22 = f2.h2;
    out.Nbar1 = f1.Nbar; out.Nbar2 = f2.Nbar;

    // rg
    double rg = rho / std::sqrt(std::max(1e-30, f1.h2 * f2.h2));
    std::vector<double> rg_del; rg_del.reserve(rho_del.size());
    for (double rdel : rho_del)
        rg_del.push_back(rdel / std::sqrt(std::max(1e-30, f1.h2 * f2.h2)));
    out.rg = rg;
    out.rg_se = se_from_delete(rg_del);
}

// ---------------------------- Output ----------------------------
static void write_summary_h2(const std::string& out, const H2Fit& f,
                             size_t snps_used, double mean_chi2, double lambda_gc, double max_chi2,
                             bool two_step, double two_step_thresh){
    std::ofstream s(out + ".ldsc.summary.txt");
    if (!s) throw std::runtime_error("Could not open summary file");
    s << "Field\tValue\n";
    s << "Mode\tH2\n";
    s << "SNPs_used\t" << snps_used << "\n";
    s << "N_bar\t" << f.Nbar << "\n";
    s << "M_5_50_total\t" << f.Mtot << "\n";
    s << "Intercept\t" << f.a << "\n";
    s << "Intercept_SE\t" << f.a_se << "\n";
    s << "h2_observed\t" << f.h2 << "\n";
    s << "h2_observed_SE\t" << f.h2_se << "\n";
    s << "Mean_Chi2\t" << mean_chi2 << "\n";
    s << "Lambda_GC\t" << lambda_gc << "\n";
    s << "Max_Chi2\t" << max_chi2 << "\n";
    if (two_step) {
        s << "TwoStep\t1\n";
        s << "TwoStep_Thresh\t" << two_step_thresh << "\n";
        if (f.step1_kept >= 0) s << "TwoStep_Step1_SNPs\t" << f.step1_kept << "\n";
    } else {
        s << "TwoStep\t0\n";
    }
}

static void write_summary_rg(const std::string& out,
                             const H2Fit& f1, const H2Fit& f2, const RGFit& g,
                             size_t snps_used, double mean_z1z2, bool two_step, double two_step_thresh){
    std::ofstream s(out + ".ldsc.summary.txt");
    if (!s) throw std::runtime_error("Could not open summary file");
    s << "Field\tValue\n";
    s << "Mode\tRG\n";
    s << "SNPs_used\t" << snps_used << "\n";
    s << "N_bar_trait1\t" << f1.Nbar << "\n";
    s << "N_bar_trait2\t" << f2.Nbar << "\n";
    s << "M_5_50_total\t" << f1.Mtot << "\n";
    s << "Intercept1\t" << f1.a << "\n";
    s << "Intercept1_SE\t" << f1.a_se << "\n";
    s << "Intercept2\t" << f2.a << "\n";
    s << "Intercept2_SE\t" << f2.a_se << "\n";
    s << "Intercept12\t" << g.a12 << "\n";
    s << "Intercept12_SE\t" << g.a12_se << "\n";
    s << "h2_trait1\t" << f1.h2 << "\n";
    s << "h2_trait1_SE\t" << f1.h2_se << "\n";
    s << "h2_trait2\t" << f2.h2 << "\n";
    s << "h2_trait2_SE\t" << f2.h2_se << "\n";
    s << "Genetic_covariance\t" << g.rho << "\n";
    s << "Genetic_covariance_SE\t" << g.rho_se << "\n";
    s << "Genetic_correlation\t" << g.rg << "\n";
    s << "Genetic_correlation_SE\t" << g.rg_se << "\n";
    s << "Mean_z1z2\t" << mean_z1z2 << "\n";
    if (two_step) {
        s << "TwoStep\t1\n";
        s << "TwoStep_Thresh\t" << two_step_thresh << "\n";
        if (f1.step1_kept >= 0) s << "TwoStep_Step1_SNPs_trait1\t" << f1.step1_kept << "\n";
        if (f2.step1_kept >= 0) s << "TwoStep_Step1_SNPs_trait2\t" << f2.step1_kept << "\n";
    } else {
        s << "TwoStep\t0\n";
    }
}

// ------------------------------- Main -------------------------------
int run_ldsc(int argc, char** argv){
    try{
        LdscArgs args = parse_ldsc(argc, argv);
        Logger log(args.out + ".log", /*append=*/true);
        echo_call(log, args);

        // Load LD (single column) and weights
        log.log("Reading reference LD ...");
        Table refld = args.ref_ld.empty()
            ? read_ld_fileset_1col(args.ref_ld_chr, ".l2.ldscore", /*for_weights=*/false)
            : read_ld_fileset_1col(args.ref_ld, "", /*for_weights=*/false);
        log.log("Ref-LD: ", refld.snp.size(), " rows, col=", refld.cols[0]);

        log.log("Reading regression weights (w-LD) ...");
        Table wld = args.w_ld.empty()
            ? read_ld_fileset_1col(args.w_ld_chr, ".l2.ldscore", /*for_weights=*/true)
            : read_ld_fileset_1col(args.w_ld, "", /*for_weights=*/true);
        if (wld.cols.size()!=1) throw std::runtime_error("w-LD must be exactly one column.");
        log.log("w-LD: ", wld.snp.size(), " rows, col=", wld.cols[0]);

        // Determine M_total
        long long Mtot = -1;
        if (!args.M.empty()){
            std::ifstream mf(args.M);
            if (mf){ double v; mf>>v; if (mf) Mtot = (long long)std::llround(v); }
            if (Mtot < 0) Mtot = std::llround(std::strtod(args.M.c_str(), nullptr));
        }
        if (Mtot <= 0) Mtot = detect_M_total(args.ref_ld_chr, refld.snp.size());
        if (Mtot <= 0) Mtot = (long long)refld.snp.size();
        log.log("Using M_5_50 total = ", (long long)Mtot);

        // --- Mode: h2 ---
        if (!args.h2.empty()){
            Sumstats ss = read_sumstats(args.h2, log);
            Merge1 m = merge_univ(ss, refld, wld);

            // QC summaries
            std::vector<double> chisq; chisq.reserve(m.y_chi2.size());
            for (double c: m.y_chi2) if (std::isfinite(c)) chisq.push_back(c);
            if (chisq.empty()) throw std::runtime_error("No finite chi^2.");
            auto tmp = chisq; std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end());
            double mean_chi=0.0; for (double v:chisq) mean_chi+=v; mean_chi/=chisq.size();
            double med = tmp[tmp.size()/2]; double lambda_gc = med / 0.4549;
            double mx = *std::max_element(chisq.begin(), chisq.end());

            H2Fit fit{};
            if (args.two_step || args.fix_int){
                if (args.fix_int) log.log("Univariate: --no-intercept set; intercept fixed at 1.0 (Step1 skipped).");
                else              log.log("Univariate: two-step enabled; Step1 chi^2 threshold = ", args.two_step_thresh);
                jackknife_fit_univ_twostep(m, (double)Mtot, std::min(args.n_blocks,(int)m.snp.size()),
                                           args.irwls_iter, args.fix_int, args.two_step_thresh, fit);
            } else {
                jackknife_fit_univ_free(m, (double)Mtot, std::min(args.n_blocks,(int)m.snp.size()),
                                        args.irwls_iter, fit);
            }

            write_summary_h2(args.out, fit, m.snp.size(), mean_chi, lambda_gc, mx, args.two_step||args.fix_int, args.two_step_thresh);
            log.log("Total Observed scale h2: ", fit.h2, " (", fit.h2_se, ")");
            log.log("Intercept: ", fit.a, " (", fit.a_se, ")");
            if (args.two_step && fit.step1_kept>=0) log.log("Two-step: Step1 kept ", fit.step1_kept, " SNPs (chi^2 ≤ ", args.two_step_thresh, ").");
            log.log("Lambda GC: ", lambda_gc, " ; Mean Chi^2: ", mean_chi, " ; Max Chi^2: ", mx);
            log.log("ldsc h2: completed.");
            return 0;
        }

        // --- Mode: rg ---
        Sumstats s1 = read_sumstats(args.rg1, log);
        Sumstats s2 = read_sumstats(args.rg2, log);
        Merge1 m1 = merge_univ(s1, refld, wld);
        Merge1 m2 = merge_univ(s2, refld, wld);
        Merge2 mb = merge_bivar(s1, s2, refld, wld);

        auto qc_mean_lambda = [&](const std::vector<double>& z)->std::pair<double,double>{
            std::vector<double> chi; chi.reserve(z.size());
            for (double v: z){ double c=v*v; if (std::isfinite(c)) chi.push_back(c); }
            auto tmp = chi; std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end());
            double mean=0.0; for (double v:chi) mean+=v; mean/=chi.size();
            double med=tmp[tmp.size()/2]; return {mean, med/0.4549};
        };
        auto [mean1, lam1] = qc_mean_lambda(s1.Z);
        auto [mean2, lam2] = qc_mean_lambda(s2.Z);
        log.log("Trait1: mean chi^2=", mean1, " ; lambdaGC=", lam1);
        log.log("Trait2: mean chi^2=", mean2, " ; lambdaGC=", lam2);

        // Univariate fits first (benefit from two-step if enabled)
        H2Fit f1{}, f2{};
        if (args.two_step || args.fix_int) {
            jackknife_fit_univ_twostep(m1, (double)Mtot, std::min(args.n_blocks,(int)m1.snp.size()),
                                       args.irwls_iter, args.fix_int, args.two_step_thresh, f1);
            jackknife_fit_univ_twostep(m2, (double)Mtot, std::min(args.n_blocks,(int)m2.snp.size()),
                                       args.irwls_iter, args.fix_int, args.two_step_thresh, f2);
        } else {
            jackknife_fit_univ_free(m1, (double)Mtot, std::min(args.n_blocks,(int)m1.snp.size()),
                                    args.irwls_iter, f1);
            jackknife_fit_univ_free(m2, (double)Mtot, std::min(args.n_blocks,(int)m2.snp.size()),
                                    args.irwls_iter, f2);
        }

        // Bivariate fit
        RGFit g{};
        jackknife_fit_bivar(mb, (double)Mtot, f1, f2,
                            std::min(args.n_blocks,(int)mb.snp.size()),
                            args.irwls_iter, args.fix_int12, g);

        // Report mean z1z2
        double mean_z1z2 = 0.0; for (size_t i=0;i<mb.snp.size();++i) mean_z1z2 += mb.Z1[i]*mb.Z2[i];
        mean_z1z2 /= (double)mb.snp.size();

        write_summary_rg(args.out, f1, f2, g, mb.snp.size(), mean_z1z2, args.two_step||args.fix_int, args.two_step_thresh);
        log.log("rg: rho_g=", g.rho, " (", g.rho_se, "), r_g=", g.rg, " (", g.rg_se, ")");
        if (!args.fix_int12) log.log("Intercept_12: ", g.a12, " (", g.a12_se, ")");
        if (args.two_step) {
            if (f1.step1_kept>=0) log.log("Two-step: Trait1 Step1 kept ", f1.step1_kept, " SNPs (chi^2 ≤ ", args.two_step_thresh, ").");
            if (f2.step1_kept>=0) log.log("Two-step: Trait2 Step1 kept ", f2.step1_kept, " SNPs (chi^2 ≤ ", args.two_step_thresh, ").");
        }
        log.log("ldsc rg: completed.");
        return 0;

    } catch (const std::exception& e){
        std::cerr << "ldsc ERROR: " << e.what() << "\n";
        return 1;
    }
}
