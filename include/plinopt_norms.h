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

#ifndef _PLINOPT_NORMS_H_
#define _PLINOPT_NORMS_H_

#include "plinopt_library.h"

// ============================================

namespace PLinOpt {
// ============================================
// ============================================

// ============================================
// Different (sparse) vector norms
template<typename Vect_t> size_t norm0(const Vect_t& v);
template<typename Vect_t> double norm1(const Vect_t& v);
template<typename Vect_t> double norm2(const Vect_t& v);
template<typename Vect_t> double norminfty(const Vect_t& v);

// ============================================
// Different Gamma factors

// (n1 L * n1 R) * abs Pij
template<typename _Mat>
std::vector<double> GPinf(const _Mat& L, const _Mat& R, const _Mat& P);

// max GPinf
template<typename _Mat>
double Ginfinf(const _Mat& L, const _Mat& R, const _Mat& P);

// n2 GPinf
template<typename _Mat>
double G2inf(const _Mat& L, const _Mat& R, const _Mat& P);

// (n2 L * n2 R) * abs Pij
template<typename _Mat>
std::vector<double> GP2(const _Mat& L, const _Mat& R, const _Mat& P);

// max GP2
template<typename _Mat>
double Ginf2(const _Mat& L, const _Mat& R, const _Mat& P);

// n2 GP2
template<typename _Mat>
double G22(const _Mat& L, const _Mat& R, const _Mat& P);

// n2 L * n2 R * n2 P
template<typename _Mat>
double G2(const _Mat& L, const _Mat& R, const _Mat& P);

// ============================================
// Different Q factors

// Qk generic factor
#define _PLO_Qk_(q0,gamma,k) ((q0)*(gamma)/std::abs((gamma)-(k)))

// Q0 factor
template<typename _Mat>
double Q0(const _Mat& L, const _Mat& R, const _Mat& P);

} // End of namespace PLinOpt
// ============================================


#include "plinopt_norms.inl"
#endif
