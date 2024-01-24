// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Optimization definitions
 ****************************************************************/

// ============================================================
// Define to print comments during parsing
//#define VERBATIM_PARSING
// ============================================================

// ============================================================
// Common Subexpressions selection:
//   will determine the maximal size of CSE
//   then among CSE with the maximal size, select:
//     by default via a score
//     otherwise define the following for random choice
//#define RANDOM_TIES
// ============================================================


#include "plinopt_library.h"
#include <linbox/algorithms/gauss.h>

#ifndef _PLINOPT_LIBRARY_OPTIMIZE_H_
#define _PLINOPT_LIBRARY_OPTIMIZE_H_


typedef std::tuple<size_t, size_t, Rational> triple;
typedef std::vector<triple> VTriple;
typedef std::map<triple,size_t> STriple;

bool operator==(const triple& u, const triple&v);


std::ostream& operator<<(std::ostream& out, const std::pair<long unsigned int, Givaro::Rational>& p);

std::ostream& operator<<(std::ostream& out, const triple& t);


std::ostream& operator<<(std::ostream& out, const VTriple& v);


template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::map<T1,T2>& v);

// Build pairs of indices with normalized coefficient (ratio of the two values)
template<typename Iter>
VTriple listpairs(const Iter& start, const Iter& end);

// If cse is present in row, add square of row density to score
size_t score(const std::vector<VTriple>& AllPairs,
             const std::vector<size_t>& Density,
             const triple& cse);

// Removing one pair
bool OneSub(Matrix& M, VTriple& multiples, size_t& nbmul,
            const char tev, const char rav);


// Factors out same coefficient in a column
template<typename Iter>
void FactorOutColumns(Matrix& T, VTriple& multiples, size_t& nbmul,
                      const char tev, const char rav,
                      const size_t j, const Iter& start, const Iter& end) ;

// Factors out same coefficient in a row
template<typename Iter>
void FactorOutRows(Matrix& M, size_t& nbadd, const char tev,
                   const size_t i, const Iter& start, const Iter& end);

// Sets new temporaries with the input values
void input2Temps(const size_t N, const char inv, const char tev);


// Global optimization function (pairs and factors)
std::pair<size_t,size_t> Optimizer(Matrix& M,
                                   const char inv, const char ouv,
                                   const char tev, const char rav);


	// Precondition _Matrix A is upper triangularized
std::pair<size_t,size_t> nullspacedecomp(Matrix& x, Matrix& A) ;


#include "plinopt_optimize.inl"
#endif
