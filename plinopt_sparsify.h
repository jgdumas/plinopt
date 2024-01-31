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


	// Find most frequent Rational;

// ============================================================


#include "plinopt_sparsify.inl"
#endif
