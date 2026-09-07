// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library
 * References:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Towards automated generation of fast and accurate algorithms
 *     for recursive matrix multiplication.
 *     J. of Symb. Comput. Vol. 134, nUM. 102524, 2026.
 *     (https://doi.org/10.1016/j.jsc.2025.102524) ]
 ****************************************************************/

#include "plinopt_norms.h"

// ============================================

namespace PLinOpt {
// ============================================
// ============================================

#ifndef GIVABS
#define GIVABS(a) ((a)>0?(a):-(a))
#endif

// ============================================
// Different (sparse) vector norms


// ============================================
// Element or second in pair
template<typename T> inline double access(const T& a) { return double(a); }
template<> inline double access(const std::pair<size_t, Givaro::Rational>& a) {
    return double(a.second);
}

// ============================================
// Different vector norms
template<typename Vect_t> inline size_t norm0(const Vect_t& v) {
    return v.size();
}

template<typename Vect_t> inline double norm1(const Vect_t& v) {
    double s(0); for(const auto& it: v) s += GIVABS(access(it));
    return s;
}

template<typename Vect_t> inline double norm2(const Vect_t& v) {
    double s(0.); for(const auto& it: v) s += access(it)*access(it);
    return std::sqrt(s);
}

template<typename Vect_t> inline double norminfty(const Vect_t& v) {
    double s(0.); for(const auto& it: v) {
        double r = GIVABS(access(it));
        if (r>s) s=r;
    }
    return s;
}

// ============================================
// Different Gamma factors



// (n1 L * n1 R) * abs Pij
template<typename _Mat>
inline std::vector<double> GPinf(const _Mat& L, const _Mat& R, const _Mat& P) {
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
template<typename _Mat>
inline double Ginfinf(const _Mat& L, const _Mat& R, const _Mat& P) {
    auto r( GPinf(L,R,P) );
    return *std::max_element(r.begin(),r.end());
}

// n2 GPinf
template<typename _Mat>
inline double G2inf(const _Mat& L, const _Mat& R, const _Mat& P) {
    return norm2(GPinf(L,R,P) );
}

// (n2 L * n2 R) * abs Pij
template<typename _Mat>
inline std::vector<double> GP2(const _Mat& L, const _Mat& R, const _Mat& P) {
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
template<typename _Mat>
inline double Ginf2(const _Mat& L, const _Mat& R, const _Mat& P) {
    auto r( GP2(L,R,P) );
    return *std::max_element(r.begin(),r.end());
}


// n2 GP2
template<typename _Mat>
inline double G22(const _Mat& L, const _Mat& R, const _Mat& P) {
    return norm2(GP2(L,R,P) );
}

// n2 L * n2 R * n2 P
template<typename _Mat>
inline double G2(const _Mat& L, const _Mat& R, const _Mat& P) {
    _Mat Pt(P.field()); PLinOpt::Transpose(Pt,P);
    double s(0.); for(size_t i(0); i<P.coldim(); ++i) {
        s += norm2(L[i])*norm2(R[i])*norm2(Pt[i]);
    }
    return s;
}


// ============================================
// Q0 factor
template<typename _Mat>
inline double Q0(const _Mat& L, const _Mat& R, const _Mat& P) {

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

} // End of namespace PLinOpt
// ============================================
