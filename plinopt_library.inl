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
    T.resize(0,0); T.resize(A.coldim(), A.rowdim());
    if (A.rowdim() != 0)
      for(auto it = A.IndexedBegin(); it != A.IndexedEnd(); ++it)
        T.setEntry(it.colIndex(),it.rowIndex(), it.value());
    return T;
}

template<typename _Mat1, typename _Mat2>
inline _Mat1& NegTranspose(_Mat1& T, const _Mat2& A) {
    T.resize(0,0); T.resize(A.coldim(), A.rowdim());
    typename _Mat1::Element tmp; T.field().init(tmp);
    if (A.rowdim() != 0)
      for(auto it = A.IndexedBegin(); it != A.IndexedEnd(); ++it)
        T.setEntry(it.colIndex(),it.rowIndex(), T.field().neg(tmp,it.value()));
    return T;
}


	// Replace row i of A, by row j of B
template<typename _Mat1, typename _Mat2>
inline _Mat1& setRow(_Mat1& A, size_t i, const _Mat2& B, size_t j) {
    using Field = typename _Mat1::Field;
    const Field& F(A.field());
    A[i].resize(0);
    for(const auto& iter: B[j]) {
        if (! F.isZero(iter.second)) A.setEntry(i,iter.first,iter.second);
    }
    return A;
}

	// Replace row i of A, by v
template<typename _Mat, typename _Vector>
inline _Mat& setRow(_Mat& A, size_t i, const _Vector& v) {
    using Field = typename _Mat::Field;
    const Field& F(A.field());
    A[i].resize(0);
    for(size_t j=0; j<v.size(); ++j) if (! F.isZero(v[j])) A.setEntry(i,j,v[j]);
    return A;
}

template<typename _Vector>
inline DenseMatrix& setRow(DenseMatrix& A, size_t i, const _Vector& v) {
    for(size_t j=0; j<v.size(); ++j) A.setEntry(i,j,v[j]);
    return A;
}


	// Negate row i of (sparse) A
template<typename _Mat>
inline _Mat& negRow(_Mat& A, size_t i) {
    using Field = typename _Mat::Field;
    const Field& F(A.field());
    for(auto& iter: A[i]) F.negin(iter.second);
    return A;
}

	// Multiply row i of (sparse) A
template<typename _Mat>
inline _Mat& mulRow(_Mat& A, size_t i, const typename _Mat::Element& e) {
    using Field = typename _Mat::Field;
    const Field& F(A.field());
    for(auto& iter: A[i]) F.mulin(iter.second, e);
    return A;
}


// M[i] <- M[i] + c * s
template<typename _Mat>
inline void opRow(_Mat& M, const size_t i, const typename _Mat::Row& s,
                  const typename _Mat::Element& c) {
    using Field = typename _Mat::Field;
    const Field& F(M.field());
    for(auto e: s) {
        const size_t j(e.first);
        typename Field::Element t; F.init(t);
        F.assign(t,M.getEntry(i,j)); // might be zero (can't use refEntry)
        F.axpyin(t,c,e.second);
        M.setEntry(i,j,t);
    }
}


	// copy dense matrix M into sparse matrix A
template<typename _Mat, typename _DMat>
inline _Mat& dense2sparse(_Mat& A, const _DMat& M) {
    A.resize(0,0); A.resize(M.rowdim(), M.coldim());
    for(size_t i=0; i<A.rowdim(); ++i) {
        setRow(A,i,M[i]);
    }
    return A;
}

	// copy sparse matrix B into sparse matrix A
template<typename _Mat>
inline _Mat& sparse2sparse(_Mat& A, const _Mat& B) {
    A.resize(0,0); A.resize(B.rowdim(), B.coldim());
    std::copy(B.rowBegin(), B.rowEnd(), A.rowBegin());
    return A;
}


	// copy matrix M into dense matrix A
template<typename _DMat, typename _SMat>
inline _DMat& any2dense(_DMat& A, const _SMat& M) {
    A.init(M.rowdim(), M.coldim());
    if (M.rowdim() != 0)
      for(auto indices = M.IndexedBegin(); indices != M.IndexedEnd(); ++indices)
        A.setEntry(indices.rowIndex(), indices.colIndex(), indices.value());
    return A;
}

	// Copy (and convert) a matrix
template<typename _Field>
inline SMatrix<_Field>& matrixCopy(SMatrix<_Field>& C,
                                   const SMatrix<_Field>& A){
    return sparse2sparse(C,A);
}

template<typename _Field>
inline SMatrix<_Field>& matrixCopy(SMatrix<_Field>& C,
                                   const DMatrix<_Field>& A){
    return dense2sparse(C,A);
}


template<typename _Field>
inline DMatrix<_Field>& matrixCopy(DMatrix<_Field>& C,
                                   const SMatrix<_Field>& A){
    return any2dense(C,A);
}

template<typename _Field>
inline DMatrix<_Field>& matrixCopy(DMatrix<_Field>& C,
                                   const DMatrix<_Field>& A){
    return any2dense(C,A);
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
    DenseMatrix dA(QQ,A.rowdim(),A.coldim()); any2dense(dA, A);
    return P.applyRight(R,dA);
}

template<>
inline Matrix& permuteRows(Matrix& R, const Permutation& P,
                    const DenseMatrix& A, const QRat& QQ) {
    DenseMatrix dR(QQ,A.rowdim(),A.coldim());
    P.applyRight(dR, A);
    return dense2sparse(R, dR);
}

template<>
inline Matrix& permuteRows(Matrix& R, const Permutation& P,
                    const Matrix& A, const QRat& QQ) {
    DenseMatrix dA(QQ,A.rowdim(),A.coldim()); any2dense(dA, A);
    DenseMatrix dR(QQ,A.rowdim(),A.coldim());
    P.applyRight(dR, dA);
    return dense2sparse(R, dR);
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





// Sets new temporaries with the input values
void input2Temps(std::ostream& sout, const size_t N,
                 const char inv, const char tev) {
    // Inputs to temporaries
    for(size_t i=0; i<N; ++i) {
        sout << tev << i << ":="
                  << inv << i << ';' << std::endl;
    }
}
// Sets new temporaries with the input values
template<typename _Mat>
void input2Temps(std::ostream& sout, const size_t N,
                 const char inv, const char tev,
                 const _Mat& trsp) {
    // Inputs to temporaries
    for(size_t i=0; i<N; ++i) {
        if (trsp[i].size()) {
            sout << tev << i << ":="
                      << inv << i << ';' << std::endl;
        } // otherwise variable is not used
    }
}




template<typename Ring>
std::ostream& printmulorjustdiv(std::ostream& out,
                                const char c, const size_t i,
                                const typename Ring::Element& e,
                                size_t& nbmul, const Ring& F) {
    out << c << i;
    if (notAbsOne(F,e)) {
        ++nbmul;
        out << '*' << e;
    }
    return out;
}

template<>
std::ostream& printmulorjustdiv(std::ostream& out,
                                const char c, const size_t i,
                                const Givaro::Rational& r,
                                size_t& nbmul, const QRat& QQ) {
    out << c << i;
    if (!QQ.isOne(r)) {
        ++nbmul;
        if (Givaro::isOne(r.nume()))
            out << '/' << r.deno();
        else
            out << '*' << r;
    }
    return out;
}




std::ostream& printSCA(std::ostream& out,
                       const char c, const size_t i, const char p,
                       const Givaro::Rational& r,
                       size_t& nbmul, const QRat& QQ) {
    if (p == '*')
        return printmulorjustdiv(out, c, i, r, nbmul, QQ);
    else {
        Givaro::Rational u; QQ.inv(u,r);
        return printmulorjustdiv(out, c, i, u, nbmul, QQ);
    }
}
