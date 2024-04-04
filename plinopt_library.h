// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library
 * References:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *     Strassen's algorithm is not optimally accurate
 *     (https://hal.science/hal-04441653) ]
 *   [ J-G. Dumas, B. Grenet; Jul. 2023
 *     In-place accumulation of fast multiplication formulae
 *     (https://hal.science/hal-04167499) ]
 ****************************************************************/

#ifndef _PLINOPT_LIBRARY_H_
#define _PLINOPT_LIBRARY_H_

#include <iostream>

// ============================================

#ifndef DEBUG
// otherwise DenseMatrix.resize issues a warning when resizing a matrix
#  define NDEBUG
#endif

#ifdef RANDOM_TIES
#  define DORANDOMSEARCH true
#else
#  define DORANDOMSEARCH false
#endif

#ifndef DEFAULT_RANDOM_LOOPS
#  define DEFAULT_RANDOM_LOOPS 100u
#endif

#include <givaro/givrational.h>
#include <givaro/givprint.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/blackbox/permutation.h>

// ============================================

using LinBox::Tag::FileFormat;

typedef Givaro::QField<Givaro::Rational> QRat;

typedef LinBox::MatrixStream<QRat> QMstream;

typedef LinBox::SparseMatrix<QRat,
                             LinBox::SparseMatrixFormat::SparseSeq > Matrix;
typedef LinBox::DenseMatrix<QRat> DenseMatrix;
typedef LinBox::Permutation<QRat> Permutation;

typedef std::vector<Givaro::Rational> QArray;
typedef LinBox::DenseVector<QRat> QVector;

template<typename _T> using Pair = std::pair<_T,_T>;

// ============================================

	// Copy the transposed  matrix
template<typename _Mat1, typename _Mat2>
inline _Mat1& Transpose(_Mat1& T, const _Mat2& A);

	// Copy the negation of the transposed  matrix
template<typename _Mat1, typename _Mat2>
inline _Mat1& NegTranspose(_Mat1& T, const _Mat2& A);

	// Copy (and convert) a matrix
template<typename _Mat1, typename _Mat2, typename _Field>
_Mat1& matrixCopy(_Mat1&, const _Mat2&, const _Field&);

	// Replace row i of A, by row j of B
template<typename _Mat1, typename _Mat2, typename _Field>
_Mat1& setRow(_Mat1& A, size_t i, const _Mat2& B, size_t j, const _Field&);

	// Replace row i of A, by v
template<typename _Mat, typename _Vector, typename _Field>
_Mat& setRow(_Mat& A, size_t i, const _Vector& v, const _Field& QQ);

	// Negate row i of A
template<typename _Mat, typename _Field>
_Mat& negRow(_Mat& A, size_t i, const _Field& F);

    // M[i] <- M[i] + c * s
template<typename _Mat, typename Field>
inline void opRow(_Mat& M, const size_t i, const typename _Mat::Row& s,
                  const typename Field::Element& c, const Field& F);

	// permute rows
template<typename _Mat1, typename _Mat2>
_Mat1& permuteRows(_Mat1& R, const Permutation& P, const _Mat2& A,
                   const QRat& QQ);

// ============================================

	// Printing tuples
typedef std::tuple<size_t, size_t, size_t> Tricounter;
std::ostream& operator<<(std::ostream& out, const Tricounter& t) {
    return out << '{' << std::get<0>(t)<< ','
               << std::get<1>(t)<< ',' << std::get<2>(t) << '}';
}

// ============================================

    // From HM representation(L,R,P) to represented matrix-multiplication
template<typename _Mat>
Tricounter LRP2MM(const _Mat& L, const _Mat& R, const _Mat& P) {
    const size_t n(std::sqrt(R.coldim()*P.rowdim()/L.coldim()));
    return Tricounter { n, R.coldim()/n, P.rowdim()/n };
}


// ============================================

    // Tools testing units
template<typename Ring>
bool isAbsOne(const Ring& F, const typename Ring::Element& e) {
    return (F.isOne(e) || F.isMOne(e));
}

template<typename Ring>
bool notAbsOne(const Ring& F, const typename Ring::Element& e) {
    return ( (!F.isOne(e)) && (!F.isMOne(e)) );
}


Givaro::Rational Fabs(const QRat& QQ, const Givaro::Rational& r) {
    return abs(r);
}

template<typename Element>
Element Fabs(const Givaro::Modular<Element>& F, const Element& e) {
    Element a; F.init(a);
    F.neg(a,e);
    return (a<e ? a : F.assign(a,e));
}

int Fsign(const QRat& QQ, const Givaro::Rational& r) {
    return sign(r);
}

template<typename Element>
int Fsign(const Givaro::Modular<Element>& F, const Element& e) {
    if (F.isZero(e)) return 0;
    Element a; F.init(a);
    F.neg(a,e);
    return (a<e ? -1 : 1);
}

// ============================================


// Vector of factorial up to n
std::vector<long> factorial(const size_t n) {
    std::vector<long> a(n); a[0]=1;
    for(size_t j=1; j<n; ++j) { a[j] = a[j-1]*(j+1); } return a;
}

// Returns the kth-permutation of 0..(n-1) via the factoradic
// Fn is the factorial vector up to n-1
std::vector<size_t> kthpermutation(const size_t k, const size_t n,
                                   const std::vector<long>& Fn );

// ============================================

#include "plinopt_library.inl"
#endif
