// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Sparsification definitions
 * References:
 *   [ G. Beniamini, N. Cheng, O. Holtz, E. Karstadt, O. Schwartz;
 *     Sparsifying the Operators of Fast Matrix Multiplication
 *     Algorithms (https://arxiv.org/abs/2008.03759) ]
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *     Strassen's algorithm is not optimally accurate
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/

// ============================================================
// Define to print comments during parsing
//#define VERBATIM_PARSING
// ============================================================

// ============================================================
// Defalt maximal depth of search:
//   number of coefficients tried within linear combinations,
//   at each step
//   (search within COEFFICIENT_SEARCH^4 possibilities)
// ============================================================
#ifndef COEFFICIENT_SEARCH
#define COEFFICIENT_SEARCH 11u
#endif
// ============================================================

#include "plinopt_library.h"
#include <linbox/algorithms/gauss.h>


#ifndef _PLINOPT_LIBRARY_SPARSIFY_H_
#define _PLINOPT_LIBRARY_SPARSIFY_H_

// ============================================================
// Sparsifying a Matrix M, in-place
//   uses a limited number of coefficients for the linear comb.
//   (that is at most: COEFFICIENT_SEARCH)
//   Returns the (inverse) change of basis and the sparsified M
Matrix& Sparsifier(Matrix& TC, Matrix& TM);
// ============================================================


// ============================================================
// Factoring out some coefficients:
//    from a row of TM, towards a column of TC
Matrix& FactorDiagonals(Matrix& TC, Matrix& TM);
// ============================================================

// ============================================================
// Alternating sparsification and column factoring
//   starting with only 3 coefficents (thus -1,0,1) ...
//   ... increasing this number by increment ...
//   ... until threshold.
size_t SparseFactor(Matrix& TICoB, Matrix& TM,
                    const size_t start=3u, const size_t increment=4u,
                    const size_t threshold=COEFFICIENT_SEARCH);
// ============================================================


// ============================================================
// Utilities

    // Testing afor a potential better sparsity
std::pair<int,int>& testLinComb(std::pair<int,int>& weight,
                                Matrix& LCoB, Matrix& Cand,
                                const size_t num,
                                const QArray& w, const Matrix& TM);


    // Consistency check of M == R.C
template<typename AMatrix>
std::ostream& consistency(std::ostream& out, const AMatrix& M,
                          const Matrix& R, const DenseMatrix& C);

    // Prints out density profile of M
    // returns total density
std::ostream& densityProfile(std::ostream& out, size_t& s, const Matrix& M);

	// Computes the transposed inverse of A
Matrix& inverseTranspose(Matrix& I, const Matrix& A);

	// Computes the rank of A
size_t& rank(size_t& r, const Matrix& A);

    // Find most frequent Rational between begin and end
    // element is suppoded to be the second in a pair
template <typename Fwd>
Givaro::Rational most_frequent_element(const Fwd& begin, const Fwd& end);

// ============================================================

#include "plinopt_sparsify.inl"
#endif
