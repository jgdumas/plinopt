// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
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

template<typename _Mat1, typename _Mat2>
inline _Mat1& NegTranspose(_Mat1& T, const _Mat2& A) {
    T.resize(0,0);
    T.resize(A.coldim(), A.rowdim());
    typename _Mat1::Element tmp; T.field().init(tmp);
    for(auto it = A.IndexedBegin(); it != A.IndexedEnd(); ++it)
        T.setEntry(it.colIndex(),it.rowIndex(), T.field().neg(tmp,it.value()));
    return T;
}


	// Replace row i of A, by row j of B
template<typename _Mat1, typename _Mat2, typename _Field>
inline _Mat1& setRow(_Mat1& A, size_t i, const _Mat2& B, size_t j,
                     const _Field& F) {
    A[i].resize(0);
    for(const auto& iter: B[j]) {
        if (! F.isZero(iter.second)) A.setEntry(i,iter.first,iter.second);
    }
    return A;
}

	// Replace row i of A, by v
template<typename _Mat, typename _Vector, typename _Field>
inline _Mat& setRow(_Mat& A, size_t i, const _Vector& v, const _Field& F) {
    A[i].resize(0);
    for(size_t j=0; j<v.size(); ++j) if (! F.isZero(v[j])) A.setEntry(i,j,v[j]);
    return A;
}

template<typename _Vector, typename _Field>
inline DenseMatrix& setRow(DenseMatrix& A, size_t i,
                           const _Vector& v, const _Field& F) {
    for(size_t j=0; j<v.size(); ++j) A.setEntry(i,j,v[j]);
    return A;
}


	// Negate row i of (sparse) A
template<typename _Mat, typename _Field>
inline _Mat& negRow(_Mat& A, size_t i, const _Field& F) {
    for(auto& iter: A[i]) F.negin(iter.second);
    return A;
}


// M[i] <- M[i] + c * s
template<typename _Mat, typename Field>
inline void opRow(_Mat& M, const size_t i, const typename _Mat::Row& s,
                  const typename Field::Element& c, const Field& F) {
    for(auto e: s) {
        const size_t j(e.first);
        typename Field::Element t; F.init(t);
        F.assign(t,M.getEntry(i,j)); // might be zero (can't use refEntry)
        F.axpyin(t,c,e.second);
        M.setEntry(i,j,t);
    }
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
template<typename _Mat>
inline _Mat& sparse2sparse(_Mat& A, const _Mat& B) {
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



// Returns the kth-permutation of 0..(n-1) via the factoradic
// Fn is the factorial vector up to n-1 -- from 'factorial(m)'
std::vector<size_t> kthpermutation(const size_t k, const size_t n,
                                   const std::vector<long>& Fn ) {
    std::vector<size_t> l;
    std::vector<size_t> Q(n);
    std::iota(Q.begin(), Q.end(), 0); // 0,1,2,3,4
    ldiv_t result;
    long r(k);
    for(size_t i=n-1; i>0; --i) {
        result = div(r, Fn[i-1]);
        l.push_back(Q[result.quot]);
        Q.erase(Q.begin()+result.quot);
        r = result.rem;
    }
    l.push_back(Q.front());
    return l;
}
