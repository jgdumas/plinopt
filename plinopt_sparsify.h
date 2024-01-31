// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Sparsification definitions
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
#define COEFFICIENT_SEARCH 11
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
Matrix& Sparsifier(Matrix& CoB, Matrix& M);
// ============================================================


// ============================================================
// Factoring out some coefficients:
//    from a column of M, towards a row of CoB
Matrix& FactorDiagonals(Matrix& CoB, Matrix& M);
// ============================================================



// ============================================================
// Utilities

    // Consistency check of M == R.C
std::ostream& consistency(std::ostream& out, const Matrix& M,
                          const DenseMatrix& R, const Matrix& C);

    // Prints out density profile of M
    // returns total density
std::ostream& densityProfile(std::ostream& out, size_t& s, const Matrix& M);


	// Find most frequent Rational between begin and end
    // element is suppoded to be the second in a pair
template <typename Fwd>
Givaro::Rational most_frequent_element(const Fwd& begin, const Fwd& end);

	// Replace row i of A, by v
template<typename Vector>
Matrix& setRow(Matrix& A, size_t i, const Vector& v, const QRat& QQ);

	// copy dense matrix M into sparse matrix A
Matrix& dense2sparse(Matrix& A, const DenseMatrix& M, const QRat& QQ);

	// copy sparse matrix B into sparse matrix A
Matrix& sparse2sparse(Matrix& A, const Matrix& B);

	// v is augmented by i, -i, 1/i and -1/i
template<typename Vector>
void augment(Vector& v, const Givaro::Integer& i, const QRat& QQ);

	// Computes the rank of A
size_t& rank(size_t& r, const Matrix& A);

// ============================================================


#include "plinopt_sparsify.inl"
#endif
