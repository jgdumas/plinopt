// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Optimization definitions
 * Reference:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
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


#ifndef _PLINOPT_LIBRARY_OPTIMIZE_H_
#define _PLINOPT_LIBRARY_OPTIMIZE_H_

#include "plinopt_library.h"
#include <givaro/modular.h>
#include <linbox/algorithms/gauss.h>

namespace PLinOpt {
// ===============================================================
template<typename Ring>
using Etriple = std::tuple<size_t, size_t, typename Ring::Element>;

	// Printing tuples
template<typename Ring>
std::ostream& printEtriple(std::ostream& out,
                           const Ring& R, const Etriple<Ring>& t) {
    return R.write(out << '{' << std::get<0>(t)<< ','
                   << std::get<1>(t)<< ',', std::get<2>(t)) << '}';
}

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

// Factors out triangles:
// < ab | b > replaced by < 0 | b | b > then < 0 | 0 | 0 | 1>
// < a  | . >             < 0 | . | 1 >      < 0 | . | 1 | 0>
// With only 2 multiplications, by a, then by b, instead of 3
template<typename triple, typename _Mat>
bool Triangle(std::ostream& sout, _Mat& M, _Mat& T,
              std::vector<triple>& multiples, size_t& nbadd, size_t& nbmul,
              const char tev, const char rav, const size_t j);


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


} // End of namespace PLinOpt
// ============================================

#include "plinopt_optimize.inl"
#endif
