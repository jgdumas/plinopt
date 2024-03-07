// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, inline implementations
 ****************************************************************/

#include "plinopt_library.h"

template<typename _Mat1, typename _Mat2>
inline _Mat1& Transpose(_Mat1& T, const _Mat2& A) {
    T.resize(0,0);
    T.resize(A.coldim(), A.rowdim());
    for(auto it = A.IndexedBegin(); it != A.IndexedEnd(); ++it)
        T.setEntry(it.colIndex(),it.rowIndex(), it.value());
    return T;
}

inline Matrix& NegTranspose(Matrix& T, const Matrix& A) {
    T.resize(0,0);
    T.resize(A.coldim(), A.rowdim());
    typename Matrix::Element tmp; T.field().init(tmp);
    for(auto it = A.IndexedBegin(); it != A.IndexedEnd(); ++it)
        T.setEntry(it.colIndex(),it.rowIndex(), T.field().neg(tmp,it.value()));
    return T;
}


	// Replace row i of A, by row j of B
inline Matrix& setRow(Matrix& A, size_t i, const Matrix& B, size_t j,
                      const QRat& QQ) {
    A[i].resize(0);
    for(const auto& iter: B[j]) {
        if (! QQ.isZero(iter.second)) A.setEntry(i,iter.first,iter.second);
    }
    return A;
}

	// Replace row i of A, by v
template<typename Vector>
inline Matrix& setRow(Matrix& A, size_t i, const Vector& v, const QRat& QQ) {
    A[i].resize(0);
    for(size_t j=0; j<v.size(); ++j) {
        if (! QQ.isZero(v[j])) A.setEntry(i,j,v[j]);
    }
    return A;
}

template<typename Vector>
inline DenseMatrix& setRow(DenseMatrix& A, size_t i,
                           const Vector& v, const QRat& QQ) {
    for(size_t j=0; j<v.size(); ++j) {
        A.setEntry(i,j,v[j]);
    }
    return A;
}


	// Negate row i of A
inline Matrix& negRow(Matrix& A, size_t i, const QRat& QQ) {
    for(auto& iter: A[i]) {
        QQ.negin(iter.second);
    }
    return A;
}

	// copy dense matrix M into sparse matrix A
inline Matrix& dense2sparse(Matrix& A, const DenseMatrix& M, const QRat& QQ) {
    A.resize(M.rowdim(), M.coldim());
    for(size_t i=0; i<A.rowdim(); ++i) {
        setRow(A,i,M[i],QQ);
    }
    return A;
}

	// copy sparse matrix M into dense matrix A
template<typename _Mat>
inline DenseMatrix& sparse2dense(DenseMatrix& A, const _Mat& M) {
    A.resize(M.rowdim(), M.coldim());
    for(auto indices = M.IndexedBegin(); indices != M.IndexedEnd(); ++indices) {
        A.setEntry(indices.rowIndex(), indices.colIndex(), indices.value());
    }
    return A;
}

	// copy sparse matrix B into sparse matrix A
inline Matrix& sparse2sparse(Matrix& A, const Matrix& B) {
    A.resize(B.rowdim(), B.coldim());
    std::copy(B.rowBegin(), B.rowEnd(), A.rowBegin());
    return A;
}


template<>
inline Matrix& matrixCopy(Matrix& C, const Matrix& A, const QRat& QQ) {
    return sparse2sparse(C,A);
}

template<>
inline Matrix& matrixCopy(Matrix& C, const DenseMatrix& A, const QRat& QQ) {
    return dense2sparse(C,A,QQ);
}

template<>
inline DenseMatrix& matrixCopy(DenseMatrix& C, const Matrix& A, const QRat& QQ) {
    return sparse2dense(C,A);
}


template<>
inline DenseMatrix& matrixCopy(DenseMatrix& C, const DenseMatrix& A, const QRat& QQ) {
    return sparse2dense(C,A);
}


	// permute rows
template<>
inline DenseMatrix& permuteRows(DenseMatrix& R, const Permutation& P,
                         const DenseMatrix& A, const QRat& QQ) {
    return P.applyRight(R,A);
}

template<>
inline DenseMatrix& permuteRows(DenseMatrix& R, const Permutation& P,
                         const Matrix& A, const QRat& QQ) {
    DenseMatrix dA(QQ,A.rowdim(),A.coldim()); matrixCopy(dA, A, QQ);
    return P.applyRight(R,dA);
}

template<>
inline Matrix& permuteRows(Matrix& R, const Permutation& P,
                    const DenseMatrix& A, const QRat& QQ) {
    DenseMatrix dR(QQ,A.rowdim(),A.coldim());
    P.applyRight(dR, A);
    return dense2sparse(R, dR, QQ);
}

template<>
inline Matrix& permuteRows(Matrix& R, const Permutation& P,
                    const Matrix& A, const QRat& QQ) {
    DenseMatrix dA(QQ,A.rowdim(),A.coldim()); matrixCopy(dA, A, QQ);
    DenseMatrix dR(QQ,A.rowdim(),A.coldim());
    P.applyRight(dR, dA);
    return dense2sparse(R, dR, QQ);
}
