// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
// ==========================================================================

/**********************************************************************
 * Computes the growth factors fro different norms
 * Usage:   L.sms R.sms P.sms
 * References:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Towards automated generation of fast and accurate algorithms
 *     for recursive matrix multiplication.
 *     Journal of Symbolic Computation 134:102524, 2026.
 *     (https://hal.science/hal-04995684) ]
 **********************************************************************/

#include "plinopt_library.h"

#ifndef GIVABS
#define GIVABS(a) ((a)>0?(a):-(a))
#endif

// ============================================
// Element or second in pair
template<typename T> double access(const T& a) { return double(a); }
template<> double access(const std::pair<size_t, Givaro::Rational>& a) {
    return double(a.second);
}

// ============================================
// Different vector norms
template<typename Vect_t> size_t norm0(const Vect_t& v) {
    size_t s(0); for(const auto& it: v) ++s; return s;
}

template<typename Vect_t> double norm1(const Vect_t& v) {
    double s(0); for(const auto& it: v) s += GIVABS(access(it));
    return s;
}

template<typename Vect_t> double norm2(const Vect_t& v) {
    double s(0.); for(const auto& it: v) s += access(it)*access(it);
    return std::sqrt(s);
}

template<typename Vect_t> double norminfty(const Vect_t& v) {
    double s(0.); for(const auto& it: v) {
        double r = GIVABS(access(it));
        if (r>s) s=r;
    }
    return s;
}

// ============================================
// Gamma factors

// (n1 L * n1 R) * abs Pij
std::vector<double> GPinf(const PLinOpt::Matrix& L,
                          const PLinOpt::Matrix& R,
                          const PLinOpt::Matrix& P) {
    std::vector<double> r(P.rowdim(),0);
    for(size_t i(0); i<P.coldim(); ++i) {
        const double n1LRi(norm1(L[i])*norm1(R[i]));
        for(size_t j(0); j<P.rowdim(); ++j) {
            r[j] += n1LRi*GIVABS(double(P.getEntry(j,i)));
        }
    }
    return r;
}

// max GPinf
double Ginfinf(const PLinOpt::Matrix& L,
               const PLinOpt::Matrix& R,
               const PLinOpt::Matrix& P) {
    auto r( GPinf(L,R,P) );
    return *std::max_element(r.begin(),r.end());
}

// n2 GPinf
double G2inf(const PLinOpt::Matrix& L,
             const PLinOpt::Matrix& R,
             const PLinOpt::Matrix& P) {
    return norm2(GPinf(L,R,P) );
}

// (n2 L * n2 R) * abs Pij
std::vector<double> GP2(const PLinOpt::Matrix& L,
                        const PLinOpt::Matrix& R,
                        const PLinOpt::Matrix& P) {
    std::vector<double> r(P.rowdim(),0.);
    for(size_t i(0); i<P.coldim(); ++i) {
        const double n2LRi(norm2(L[i])*norm2(R[i]));
        for(size_t j(0); j<P.rowdim(); ++j) {
            r[j] += n2LRi*GIVABS(double(P.getEntry(j,i)));
        }
    }
    return r;
}

// max GP2
double Ginf2(const PLinOpt::Matrix& L,
             const PLinOpt::Matrix& R,
             const PLinOpt::Matrix& P) {
    auto r( GP2(L,R,P) );
    return *std::max_element(r.begin(),r.end());
}


// n2 GP2
double G22(const PLinOpt::Matrix& L,
           const PLinOpt::Matrix& R,
           const PLinOpt::Matrix& P) {
    return norm2(GP2(L,R,P) );
}

// n2 L * n2 R * n2 P
double G2(const PLinOpt::Matrix& L,
          const PLinOpt::Matrix& R,
          const PLinOpt::Matrix& P) {
    PLinOpt::Matrix Pt(P.field()); PLinOpt::Transpose(Pt,P);
    double s(0.); for(size_t i(0); i<P.coldim(); ++i) {
        s += norm2(L[i])*norm2(R[i])*norm2(Pt[i]);
    }
    return s;
}

// Q0 factor
double Q0(const PLinOpt::Matrix& L,
          const PLinOpt::Matrix& R,
          const PLinOpt::Matrix& P) {

    std::vector<double> n0LR(P.coldim(),0.);
    for(size_t i(0); i<P.coldim(); ++i) n0LR[i] = norm0(L[i])*norm0(R[i]);

    std::vector<double> r(P.rowdim(),0.);
    for(size_t j(0); j<P.rowdim(); ++j) {
        for(const auto& it: P[j]) {
            const auto& n0LRi(n0LR[it.first]);
            if (n0LRi > r[j]) r[j]= n0LRi;
        }
        r[j] += norm0(P[j]);
    }
    return *std::max_element(r.begin(),r.end());
}

#define Qk(q0,gamma,k) ((q0)*(gamma)/std::abs((gamma)-(k)))


// ===============================================================
// argv[1-3]: L.sms R.sms P.sms
int main(int argc, char ** argv) {

    if ((argc <=3) || (std::string(argv[1]) == "-h")) {
        std::clog << "Usage:" << argv[0] << " L.sms R.sms P.sms\n";
        exit(-1);
    }

        // =============================================
        // Reading matrices
	std::ifstream left (argv[1]), right (argv[2]), product(argv[3]);

    using PLinOpt::FileFormat;
    PLinOpt::QRat QQ;
    PLinOpt::QMstream ls(QQ, left), rs(QQ, right), ss(QQ, product);
    PLinOpt::Matrix L(ls), R(rs), P(ss);

    if ( (L.rowdim() != R.rowdim()) || (L.rowdim() != P.coldim()) ) {
        std::cerr << "# \033[1;31m****** ERROR, inner dimension mismatch: "
                  << L.rowdim() << "(.)" << R.rowdim() << '|' << P.coldim()
                  << " ******\033[0m"
                  << std::endl;
        return 2;
    }

#ifdef VERBATIM_PARSING
    L.write(std::clog << "L:=",FileFormat::Maple) << ';' << std::endl;
    R.write(std::clog << "R:=",FileFormat::Maple) << ';' << std::endl;
    P.write(std::clog << "P:=",FileFormat::Maple) << ';' << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

    PLinOpt::Tricounter mkn(PLinOpt::LRP2MM(L,R,P));
    const size_t& m(std::get<0>(mkn)), k(std::get<1>(mkn)), n(std::get<2>(mkn));

    std::clog <<"# Norms of "
              << m << 'x' << k << 'x' << n
              << " Matrix-Multiplication:" << std::endl;

        // =============================================
        // Different norms

    const double ginfinf(Ginfinf(L,R,P));
    const double ginf2(Ginf2(L,R,P));
    const double g2inf(G2inf(L,R,P));
    const double g22(G22(L,R,P));
    const double g2(G2(L,R,P));
    const double q0(Q0(L,R,P));

    const double qkinfinf(Qk(q0,ginfinf,k));
    const double q1inf2(Qk(q0,ginf2,1));
    const double q12inf(Qk(q0,g2inf,1));
    const double sqrtk(std::sqrt(double(k)));
    const double kth(sqrtk*sqrtk*sqrtk);
    const double qk2inf(Qk(q0,g2inf,kth));
    const double q122(Qk(q0,g22,1));


    std::clog << std::fixed << std::setw(8)
              << "#  \t\tGamma \t\tlog_" << k << std::endl;
    std::clog << "## Ginfinf:\t" << ginfinf << '\t'
              << std::log(ginfinf)/std::log(k) << std::endl;
    std::clog << "## Ginf2:\t" << ginf2 << '\t'
              << std::log(ginf2 )/std::log(k) << std::endl;
    std::clog << "## G2inf:\t" << g2inf << '\t'
              << std::log(g2inf)/std::log(k) << std::endl;
    std::clog << "## G22:\t\t" << g22 << '\t'
              << std::log(g22)/std::log(k) << std::endl;
    std::clog << "## G2:\t\t" << g2 << '\t'
              << std::log(g2)/std::log(k) << std::endl;
    std::clog << "## Q0:\t\t" << q0 << '\t'
              << std::log(q0)/std::log(k) << std::endl;
    std::clog << "## Qkinfinf:\t" << qkinfinf << '\t'
              << std::log(qkinfinf)/std::log(k) << std::endl;
    std::clog << "## Q1inf2:\t" << q1inf2 << '\t'
              << std::log(q1inf2)/std::log(k) << std::endl;
    std::clog << "## Qk12inf:\t" << q12inf << '\t'
              << std::log(q12inf)/std::log(k) << std::endl;
    std::clog << "## Qk2inf:\t" << qk2inf << '\t'
              << std::log(qk2inf)/std::log(k) << std::endl;
    std::clog << "## Q122:\t" << q122 << '\t'
              << std::log(q122)/std::log(k) << std::endl;

    return 0;
}
