// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
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


using LinBox::Tag::FileFormat;

typedef Givaro::QField<Givaro::Rational> QRat;

typedef LinBox::MatrixStream<QRat> QMstream;

typedef LinBox::SparseMatrix<QRat,
                             LinBox::SparseMatrixFormat::SparseSeq > Matrix;
typedef LinBox::DenseMatrix<QRat> DenseMatrix;
typedef LinBox::Permutation<QRat> Permutation;

typedef std::vector<Givaro::Rational> QArray;
typedef LinBox::DenseVector<QRat> QVector;

	// Copy the transposed  matrix
template<typename _Mat1, typename _Mat2>
inline _Mat1& Transpose(_Mat1& T, const _Mat2& A);

	// Copy the negation of the transposed  matrix
Matrix& NegTranspose(Matrix& T, const Matrix& A);

	// Copy (and convert) a matrix
template<typename _Mat1, typename _Mat2>
_Mat1& matrixCopy(_Mat1&, const _Mat2&, const QRat&);

	// Replace row i of A, by row j of B
Matrix& setRow(Matrix& A, size_t i, const Matrix& B, size_t j, const QRat&);

	// Negate row i of A
Matrix& negRow(Matrix& A, size_t i, const QRat& QQ);

	// Replace row i of A, by v
template<typename _Mat, typename Vector>
_Mat& setRow(_Mat& A, size_t i, const Vector& v, const QRat& QQ);

	// permute rows
template<typename _Mat1, typename _Mat2>
_Mat1& permuteRows(_Mat1& R, const Permutation& P, const _Mat2& A,
                   const QRat& QQ);

	// Printing tuples
typedef std::tuple<size_t, size_t, size_t> Tricounter;
std::ostream& operator<<(std::ostream& out, const Tricounter& t) {
    return out << '{' << std::get<0>(t)<< ','
               << std::get<1>(t)<< ',' << std::get<2>(t) << '}';
}

    // From HM representation(L,R,P) to represented matrix-multiplication
template<typename _Mat>
Tricounter LRP2MM(const _Mat& L, const _Mat& R, const _Mat& P) {
    const size_t n(std::sqrt(R.coldim()*P.rowdim()/L.coldim()));
    return Tricounter { n, R.coldim()/n, P.rowdim()/n };
}

#include "plinopt_library.inl"
#endif
