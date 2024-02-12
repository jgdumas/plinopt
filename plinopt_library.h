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

#include <givaro/givrational.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/util/matrix-stream.h>

using Givaro::Rational;
using LinBox::Tag::FileFormat;

typedef Givaro::QField<Givaro::Rational> QRat;

typedef LinBox::MatrixStream<QRat> QMstream;

typedef LinBox::SparseMatrix<QRat,
                             LinBox::SparseMatrixFormat::SparseSeq > Matrix;
typedef LinBox::DenseMatrix<QRat> DenseMatrix;

typedef std::vector<Givaro::Rational> QArray;
typedef LinBox::DenseVector<QRat> QVector;


Matrix& Transpose(Matrix& T, const Matrix& A);

Matrix& NegTranspose(Matrix& T, const Matrix& A);


#include "plinopt_library.inl"
#endif
