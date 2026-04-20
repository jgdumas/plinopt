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
 *   [ J-G. Dumas, B. Grenet;
 *     In-place accumulation of fast multiplication formulae
 *     ISSAC 2024, Raleigh, NC USA, pp. 16-25.
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


namespace PLinOpt {
// ============================================
// === PLinOpt basis types

using LinBox::Tag::FileFormat;

typedef Givaro::QField<Givaro::Rational> QRat;

typedef LinBox::MatrixStream<QRat> QMstream;

typedef LinBox::SparseMatrix<QRat,
                             LinBox::SparseMatrixFormat::SparseSeq > Matrix;
typedef LinBox::DenseMatrix<QRat> DenseMatrix;
typedef LinBox::Permutation<QRat> Permutation;



template<typename _Field>
using SMatrix=LinBox::SparseMatrix<_Field,
                                   LinBox::SparseMatrixFormat::SparseSeq >;
template<typename _Field>
using DMatrix=LinBox::DenseMatrix<_Field>;


typedef std::vector<Givaro::Rational> QArray;
typedef LinBox::DenseVector<QRat> QVector;

template<typename _T> using Pair = std::pair<_T,_T>;
inline Pair<size_t>& operator +=(Pair<size_t>& x, const Pair<size_t>& y) {
    x.first += y.first; x.second += y.second;
    return x;
}
inline Pair<size_t>& operator -=(Pair<size_t>& x, const Pair<size_t>& y) {
    x.first -= y.first; x.second -= y.second;
    return x;
}

// ============================================
// === PLinOpt matrix operations

	// Copy the transposed  matrix
template<typename _Mat1, typename _Mat2>
inline _Mat1& Transpose(_Mat1& T, const _Mat2& A);

	// Copy the negation of the transposed  matrix
template<typename _Mat1, typename _Mat2>
inline _Mat1& NegTranspose(_Mat1& T, const _Mat2& A);

	// Copy (and convert) a matrix
template<typename _Mat1, typename _Mat2>
_Mat1& matrixCopy(_Mat1&, const _Mat2&);

    // copy sparse matrix M into dense matrix A
template<typename _DMat, typename _SMat>
inline _DMat& any2dense(_DMat& A, const _SMat& M);

template<typename _SMat>
inline std::ostream& print(std::ostream& o, const _SMat& M,
                           const FileFormat _ff=FileFormat::Maple) {
    DenseMatrix T(M.field());
    any2dense(T, M);
    return T.write(o, _ff);
}

	// Replace row i of A, by row j of B
template<typename _Mat1, typename _Mat2>
_Mat1& setRow(_Mat1& A, size_t i, const _Mat2& B, size_t j);

	// Replace row i of A, by v
template<typename _Mat, typename _Vector>
_Mat& setRow(_Mat& A, size_t i, const _Vector& v);

	// Update col j of A, by v
template<typename _Mat, typename _Vector>
_Mat& updCol(_Mat& A, size_t j, const _Vector& v);

	// Negate row i of A
template<typename _Mat>
_Mat& negRow(_Mat& A, size_t i);

	// Multiply row i of (sparse) A
template<typename _Mat>
inline _Mat& mulRow(_Mat& A, size_t i, const typename _Mat::Element& e);

    // M[i] <- M[i] + c * s
template<typename _Mat>
inline void opRow(_Mat& M, const size_t i, const typename _Mat::Row& s,
                  const typename _Mat::Element& c);

	// permute rows
template<typename _Mat1, typename _Mat2>
_Mat1& permuteRows(_Mat1& R, const Permutation& P, const _Mat2& A,
                   const QRat& QQ);


// ===============================================================
// === PLinOpt Tensor product of matrices
template<typename _Matrix>
_Matrix& Tensor(_Matrix& T, const _Matrix& A, const _Matrix& B);



// ============================================
// === PLinOpt tuples

typedef std::tuple<size_t, size_t, size_t> Tricounter;

	// Printing tuples
std::ostream& operator<<(std::ostream& out, const Tricounter& t) {
    return out << '{' << std::get<0>(t)<< ','
               << std::get<1>(t)<< ',' << std::get<2>(t) << '}';
}

// ============================================

    // From HM representation(L,R,P) to represented matrix-multiplication
template<typename _Mat>
Tricounter LRP2MM(const _Mat& L, const _Mat& R, const _Mat& P) {
    const size_t n(std::sqrt(R.coldim()*P.rowdim()/L.coldim()));
    return Tricounter { P.rowdim()/n, R.coldim()/n, n };
}


// ============================================
// === PLinOpt helpers

    // Tools testing units
template<typename Ring>
bool isAbsOne(const Ring& F, const typename Ring::Element& e) {
    return (F.isOne(e) || F.isMOne(e));
}

template<typename Ring>
bool notAbsOne(const Ring& F, const typename Ring::Element& e) {
    return ( (!F.isOne(e)) && (!F.isMOne(e)) );
}

template<typename Ring>
bool notZeroNotAbsOne(const Ring& F, const typename Ring::Element& e) {
    return ( (!F.isZero(e)) && (!F.isOne(e)) && (!F.isMOne(e)) );
}

Givaro::Rational Fabs(const QRat& QQ, const Givaro::Rational& r) {
    return abs(r);
}

    // Modular genrralizations
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
// === PLinOpt Sets new temporaries with the input values
std::ostream& input2Temps(std::ostream& sout, const size_t N,
                          const char inv, const char tev);

template<typename _Mat>
std::ostream& input2Temps(std::ostream& sout, const size_t N,
                          const char inv, const char tev, const _Mat& trsp);

std::ostream& input2Temps(std::ostream& sout, const size_t N,
                          const char inv, const char tev,
                          const LinBox::LightContainer<long>& perm);

// ============================================
// === PLinOpt programs
//		prints c[i] * e, or c[i] / b for rational e=1/b
//		updates nbmul if e not in {-1,1}
template<typename Ring>
std::ostream& printmulorjustdiv(std::ostream& out,
                                const char c, const size_t i,
                                const typename Ring::Element& e,
                                size_t& nbmul, const Ring& F);



// ============================================
// === PLinOpt Vector of factorial up to n
std::vector<size_t> factorial(const size_t n) {
    std::vector<size_t> a(n); a[0]=1u;
    for(size_t j=1; j<n; ++j) { a[j] = a[j-1]*(j+1); }
    return a;
}

// ============================================
// === PLinOpt permutation
//		Returns the kth-permutation of 0..(n-1) via the factoradic
//		Fn is the factorial vector up to n-1
std::vector<size_t> kthpermutation(const size_t k, const size_t n,
                                   const std::vector<long>& Fn );

// ============================================
// === PLinOpt Generic verification, within Field, of matrices over Base
template<typename Field, typename _Matrix>
int MMchecker(const Field& FF, const size_t bitsize,
              const _Matrix& L, const _Matrix& R, const _Matrix& P);

// ============================================

} // End of namespace PLinOpt
// ============================================


#include "plinopt_library.inl"
#endif
