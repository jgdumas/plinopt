// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Sparsification definitions
 * References:
 *   [ G. Beniamini, N. Cheng, O. Holtz, E. Karstadt, O. Schwartz;
 *     Sparsifying the Operators of Fast Matrix Multiplication
 *     Algorithms (https://arxiv.org/abs/2008.03759) ]
 *   [ Gal Beniamini, Oded Schwartz;
 *     Faster Matrix Multiplication via Sparse Decomposition.
 *     SPAA 2019: 11-22]
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/


#ifndef _PLINOPT_LIBRARY_SPARSIFY_H_
#define _PLINOPT_LIBRARY_SPARSIFY_H_

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

namespace PLinOpt {

// ============================================================
// Factorizing into an Alternate basis time a CoB
template<typename _Mat>
int Factorizer(_Mat& Alt, _Mat& CoB, const _Mat& M,
               const size_t randomloops, const size_t selectinnerdim=0,
               const bool progressreport=true);
// ============================================================


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
//   starting with only start coefficents (at least -1,0,1) ...
//   ... increasing this number by increment ...
//   ... until threshold.
size_t SparseFactor(Matrix& TICoB, Matrix& TM,
                    const size_t start=3u, const size_t increment=4u,
                    const size_t threshold=COEFFICIENT_SEARCH);
// ============================================================

// ============================================================
// First:  FactorDiagonal
// Second: SparseFactor with default parameters
// Third:  SparseFactor with maxnumcoeff, 1, maxnumcoeff
Givaro::Timer& sparseAlternate(Givaro::Timer& chrono, Matrix& CoB, Matrix& Res,
                               const Matrix& M, const size_t maxnumcoeff);
// ============================================================


// ============================================================
// Sparsifying and reducing coefficient diversity of a matrix
// by sparse QLUP elimination, followed by block sparsification
int blockSparsifier(Givaro::Timer& elapsed, Matrix& CoB, Matrix& Res,
                    const Matrix& M, const size_t blocksize,
                    const FileFormat& matformat, const size_t maxnumcoeff);


// ============================================================
// Decomposing the matrix into:
//   (Res=[identity,lower part])*(CoB=[upperpart])
//   with prescribed inner dimension
template<typename _Mat>
Pair<size_t> backSolver(_Mat& Res, _Mat& CoB, const _Mat& iM);


// ============================================================
// Utilities

    // Testing afor a potential better sparsity
bool testLinComb(Pair<int>& weight, Matrix& LCoB, Matrix& Cand,
                 const size_t num, const QArray& w, const Matrix& TM);


    // Consistency check of M == R.C
template<typename _Mat1, typename _Mat2, typename _DMat>
std::ostream& consistency(std::ostream& out, const _Mat1& M,
                          const _Mat2& R, const _DMat& C);

    // Prints out density profile of M
    // returns total density
template<typename _Mat>
std::ostream& densityProfile(std::ostream& out, size_t& s, const _Mat& M);

	// Computes the transposed inverse of A
template<typename _Mat1, typename _Mat2>
_Mat1& inverseTranspose(_Mat1& TI, const _Mat2& A);

    // ============================================================
    // Compute R, s.t. A == T . R
template<typename _Mat1, typename _Mat2>
DenseMatrix& applyInverse(DenseMatrix& R, const _Mat1& T, const _Mat2& A,
                          const LinBox::MatrixDomain<QRat>& BMD);

	// Computes the rank of A
template<typename _Mat>
size_t& rank(size_t& r, const _Mat& A);


// ============================================================
// QLUP Gaussian elimination of A
//    A  <-- (QL)^{-1} . A
//    QL <-- Q.L
// Updates QL and A only if resulting A = UP is sparser
inline bool sparseLU(Matrix& QL, Matrix& A, const size_t sparsity);


// ============================================================
// QLUP Gaussian elimination of A
//    A  <--  UP = (QL)^{-1} . A
//    TC <--  (QL)^{-1} . TC
// Updates TC and A only if resulting A is sparser
inline bool sparseILU(Matrix& TC, Matrix& A, const size_t sparsity);



    // Cut matrix by blocks of columns
std::vector<Matrix>& separateColumnBlocks(std::vector<Matrix>&,
                                          const Matrix&, const size_t);

	// Build block diagonal matrix from vector of blocks
Matrix& diagonalMatrix(Matrix& M, const std::vector<Matrix>& V);

    // Append columns by blocks of columns
template<typename _Mat>
Matrix& augmentedMatrix(Matrix&, const std::vector<_Mat>&);

    // Find most frequent Rational between begin and end
    // element is suppoded to be the second in a pair
template <typename Fwd>
Givaro::Rational most_frequent_element(const Fwd& begin, const Fwd& end);


} // End of namespace PLinOpt
// ============================================

#include "plinopt_sparsify.inl"
#endif
