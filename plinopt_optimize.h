// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Optimization definitions
 * Reference: [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *              Strassen's algorithm is not optimally accurate
 *              (https://hal.science/hal-04441653) ]
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
//     with DEFAULT_RANDOM_LOOPS for the default random loops
// ============================================================


#include "plinopt_library.h"
#include <givaro/modular.h>
#include <linbox/algorithms/gauss.h>

#ifndef _PLINOPT_LIBRARY_OPTIMIZE_H_
#define _PLINOPT_LIBRARY_OPTIMIZE_H_

template<typename Ring>
using Etriple = std::tuple<size_t, size_t, typename Ring::Element>;

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::map<T1,T2>& v);

// Build pairs of indices with normalized coefficient (ratio of the two values)
template<typename Container, typename Ring>
std::vector<Etriple<Ring>> listpairs(const Container& c, const Ring& F);

// If cse is present in row, add square of row density to score
template<typename triple>
size_t score(const std::vector<std::vector<triple>>& AllPairs,
             const std::vector<size_t>& Density,
             const triple& cse);

// Direct program generateur from a matrix
template<typename _Mat, typename triple>
std::ostream& ProgramGen(std::ostream& sout, _Mat& M,
                         std::vector<triple>& multiples,
                         size_t& addcount, size_t& nbmul,
                         const char inv, const char ouv,
                         const char tev, const char rav);

// Removing one CSE
template<typename triple, typename _Mat>
bool OneSub(std::ostream& sout, _Mat& M, std::vector<triple>& multiples,
            size_t& nbmul, const char tev, const char rav);

// Recusive search for the best CSE
template<typename triple,typename _Mat>
bool RecSub(std::vector<std::string>& out, _Mat& Mat,
            std::vector<triple>& multiples,
            size_t& nbadd, size_t& nbmul, const size_t lvl,
            const char tev, const char rav);

// Factors out same coefficient in a column
template<typename Iter, typename triple, typename _Mat>
void FactorOutColumns(std::ostream& sout, _Mat& T,
                      std::vector<triple>& multiples, size_t& nbmul,
                      const char tev, const char rav,
                      const size_t j, const Iter& start, const Iter& end) ;

// Factors out same coefficient in a row
template<typename Iter, typename _Mat>
void FactorOutRows(std::ostream& sout, _Mat& M, size_t& nbadd, const char tev,
                   const size_t i, const Iter& start, const Iter& end);

// Sets new temporaries with the input values
void input2Temps(std::ostream& sout, const size_t N,
                 const char inv, const char tev);
template<typename _Mat>
void input2Temps(std::ostream& sout, const size_t N,
                 const char inv, const char tev, const _Mat& trsp);


// Global random optimization function (pairs and factors)
template<typename _Mat>
Pair<size_t> Optimizer(std::ostream& sout, _Mat& M,
                       const char inv, const char ouv,
                       const char tev, const char rav);

// Global exhaustive optimization function (pairs and factors)
template<typename _Mat>
Pair<size_t> RecOptimizer(std::ostream& sout, _Mat& M,
                          const char inv, const char ouv,
                          const char tev, const char rav);

// Precondition _Matrix A is upper triangular
template<typename _Mat>
Pair<size_t> nullspacedecomp(std::ostream& sout, _Mat& x, _Mat& A,
                             const bool mostCSE=false);

template<typename _Mat>
Pair<size_t> nullspacedecomp(std::ostream& sout, _Mat& x, _Mat& A,
                             std::vector<size_t>& l,
                             const bool mostCSE=false);



// prints c[i] * e, or c[i] / b for rational e=1/b
// updates nbmul if e not in {-1,1}
template<typename Ring>
std::ostream& printmulorjustdiv(std::ostream& out,
                                const char c, const size_t i,
                                const typename Ring::Element& e,
                                size_t& nbmul, const Ring& F);

#include "plinopt_optimize.inl"
#endif
