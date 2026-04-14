// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, inline implementations
 ****************************************************************/

#include "plinopt_library.h"

namespace PLinOpt {
// ===============================================================
// ==================
// Matrix operations

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
    auto& Ai(A[i]); const auto& Bj(B[j]);
    Ai.resize(0); Ai.reserve(Bj.size());
    for(const auto& biter: Bj) {
        if (! F.isZero(biter.second))
            Ai.emplace_back(biter.first,biter.second);
    }
    return A;
}

	// Replace row i of A, by v
template<typename _Mat, typename _Vector>
inline _Mat& setRow(_Mat& A, size_t i, const _Vector& v) {
    using Field = typename _Mat::Field;
    const Field& F(A.field());
    auto& Ai(A[i]); Ai.resize(0);
    for(size_t j=0; j<v.size(); ++j)
        if (! F.isZero(v[j])) Ai.emplace_back(j,v[j]);
    return A;
}

template<typename _Vector>
inline DenseMatrix& setRow(DenseMatrix& A, size_t i, const _Vector& v) {
    for(size_t j=0; j<v.size(); ++j) A.setEntry(i,j,v[j]);
    return A;
}


	// Update col j of A, by v
template<typename _Mat, typename _Vector>
inline _Mat& updCol(_Mat& A, size_t j, const _Vector& v) {

    using Field = typename _Mat::Field;
    const Field& F(A.field());
    for(size_t i=0; i<v.size(); ++i)
        if (! F.isZero(v[i])) A.setEntry(i,j,v[i]);
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
    _Mat B(F,1,M.coldim()); setRow(B,0,M,i); // copy row i of M
    for(auto e: s) {
        const size_t j(e.first);
        typename Field::Element t; F.init(t);
        F.assign(t,M.getEntry(i,j));    // might be zero (can't use refEntry)
        F.axpyin(t,c,e.second);
        B.setEntry(0,j,t);              // might have become zero ...
    }
    setRow(M,i,B,0);                    // ... gets only the non-zeroes in B
}


	// copy dense matrix M into sparse matrix A
template<typename _Mat, typename _DMat>
inline _Mat& dense2sparse(_Mat& A, const _DMat& M) {
    A.resize(M.rowdim(), M.coldim());
    for(size_t i=0; i<A.rowdim(); ++i) {
        setRow(A,i,M[i]);
    }
    return A;
}

	// copy sparse matrix B into sparse matrix A
template<typename _Mat>
inline _Mat& sparse2sparse(_Mat& A, const _Mat& B) {
    A.resize(B.rowdim(), B.coldim());
    for(size_t i=0; i<A.rowdim(); ++i) {
        setRow(A,i,B,i);
    }
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

template<> inline
Matrix& permuteRows(Matrix& R, const Permutation& P,
                    const Matrix& A, const QRat& QQ) {
    DenseMatrix dA(QQ,A.rowdim(),A.coldim()); any2dense(dA, A);
    DenseMatrix dR(QQ,A.rowdim(),A.coldim());
    P.applyRight(dR, dA);
    return dense2sparse(R, dR);
}


// ===============================================================
// Tensor product of matrices
template<typename _Matrix>
_Matrix& Tensor(_Matrix& T, const _Matrix& A, const _Matrix& B) {
    T.resize(A.rowdim()*B.rowdim(), A.coldim()*B.coldim());
    using Element=typename _Matrix::Element;
    Element tmp; A.field().init(tmp);
    for(auto itA = A.IndexedBegin(); itA != A.IndexedEnd(); ++itA) {
    for(auto itB = B.IndexedBegin(); itB != B.IndexedEnd(); ++itB) {
        T.setEntry(itA.rowIndex()*B.rowdim()+itB.rowIndex(),
                   itA.colIndex()*B.coldim()+itB.colIndex(),
                   A.field().mul(tmp,itA.value(),itB.value()));
    }
    }
    return T;
}

// Matrix operations
template<typename _Mat> inline
Pair<size_t> naiveOps(const _Mat& M) {
    const auto& FF(M.field());
    size_t nbadd(0), nbmul(0);
    for(auto iter=M.rowBegin(); iter != M.rowEnd(); ++iter)
        nbadd += std::max((long)iter->size()-1l,0l);
    for(auto it = M.IndexedBegin(); it != M.IndexedEnd(); ++it)
        if (notAbsOne(FF,it.value())) ++nbmul;
    return Pair<size_t>{nbadd,nbmul};
}

	// Number of non-zero elements (sparse)
template<typename _Matrix> inline
size_t density(const _Matrix& A) {
    size_t nnz(0);
    for(auto row=A.rowBegin(); row != A.rowEnd(); ++row) {
        nnz += row->size(); // # of non-zero elements
    }
    return nnz;
}

	// Number of non-zero elements (dense)
template<typename Field> inline
size_t density(const LinBox::DenseMatrix<Field>& A) {
    size_t nnz(0);
    for(auto itA = A.IndexedBegin(); itA != A.IndexedEnd(); ++itA) {
        if (! A.field().isZero(itA.value())) ++nnz;
    }
    return nnz;

}

template<typename _Mat> inline
Pair<size_t> nonzeroes(const _Mat& M) {
    const auto& F(M.field());
    size_t nnz(0), nno(0);
    for(auto it=M.IndexedBegin();it!=M.IndexedEnd();++it) {
        if (! F.isZero(it.value()) ) ++nnz;
        if (! isAbsOne(F, it.value())) ++nno;
    }
    return Pair<size_t>{nnz,nno};
}

// ==================
// Returns the kth-permutation of 0..(n-1) via the factoradic
// Fn is the factorial vector up to n-1 -- from 'factorial(m)'
inline std::vector<size_t> kthpermutation(const size_t k, const size_t n,
                                          const std::vector<size_t>& Fn ) {
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


// ==================
// Sets new temporaries with the input values
inline std::ostream& input2Temps(std::ostream& sout, const size_t N,
                                 const char inv, const char tev) {
        // Inputs to temporaries
    for(size_t i=0; i<N; ++i) {
        sout << tev << i << ":="
             << inv << i << ';' << std::endl;
    }
    return sout;
}
// Sets new temporaries with the input values, checking usage
template<typename _Mat> inline
std::ostream& input2Temps(std::ostream& sout, const size_t N,
                          const char inv, const char tev,
                          const _Mat& trsp) {
        // Inputs to temporaries
    for(size_t i=0; i<N; ++i) {
        if (trsp[i].size()) {
            sout << tev << i << ":="
                 << inv << i << ';' << std::endl;
        } // otherwise variable is not used
    }
    return sout;
}
// ==================
// Sets new temporaries with the permuted input values
inline std::ostream& input2Temps(std::ostream& sout, const size_t N,
                                 const char inv, const char tev,
                                 const LinBox::LightContainer<long>& P) {
    // Inputs to temporaries
    for(size_t i=0; i<N; ++i) {
        sout << tev << i << ":="
                  << inv << P[i] << ';' << std::endl;
    }
    return sout;
}

// ==================
// Program outputs
template<typename Ring> inline
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

template<> inline
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

inline std::ostream& printSCA(std::ostream& out,
                              const char c, const size_t i, const char p,
                              const Givaro::Rational& r,
                              size_t& nbmul, const QRat& QQ) {
    if (p == '*')
        return printmulorjustdiv(out, c, i, r, nbmul, QQ);
    else { // operation is division --> multiply by the inverse
        Givaro::Rational u; QQ.inv(u,r);
        return printmulorjustdiv(out, c, i, u, nbmul, QQ);
    }
}


} // End of namespace PLinOpt
// ============================================


#include <givaro/modular.h>
#include <givaro/givquotientdomain.h>

// ===============================================================
// Generic random element with bitsize
template<typename Domain> inline
typename Domain::Element& RandomElt(typename Domain::Element& e,
                                    const Domain& D,
                                    Givaro::GivRandom& generator,
                                    const size_t bitsize) {
    return D.random(generator, e, bitsize);
}

// Specialization for Modular, as bitsize is meanigless for its random
template<> inline
typename Givaro::Integer& RandomElt(Givaro::Integer& e,
                                    const Givaro::Modular<Givaro::Integer>& D,
                                    Givaro::GivRandom& generator,
                                    const size_t bitsize) {
    return D.random(generator, e);
}

// Generic {-1,0,1} element
template<typename Domain>
typename Domain::Element& zoRandomElt(typename Domain::Element& e,
                                      const Domain& D,
                                      Givaro::GivRandom& generator) {
    D.init(e, (generator() % 3));
    return D.subin(e, D.one);
}
// ===============================================================


// ===============================================================
// Missing LinBox specialization of quotient homomorphism
template <class Field>
class LinBox::Hom<Field, Givaro::QuotientDom<Field> > {
public:
    typedef Field Source;
    typedef typename Givaro::QuotientDom<Field> Target;
    typedef typename Source::Element SrcElt;
    typedef typename Target::Element Elt;

    Hom(const Source& S, const Target& T) : _source(S), _target(T) {}

    Elt& image(Elt& t, const SrcElt& s) {
        _target.assign (t, s); // this will compute the modular image
        return t;
    }

    SrcElt& preimage(SrcElt& s, const Elt& t) {
        _source.assign (s, t); // this will just copy
        return s;
    }

    const Source& source() { return _source;}
    const Target& target() { return _target;}

private:
    const Source& _source;
    const Target& _target;
}; // end Hom
// ===============================================================

namespace PLinOpt {
// ===============================================================
// Generic verification, within Field, of matrices over Base
template<typename Field, typename _Matrix>
int MMchecker(const Field& FF, const size_t bitsize,
              const _Matrix& L, const _Matrix& R, const _Matrix& P) {
    using PLinOpt::FileFormat;
    using FMatrix=_Matrix;
    typedef LinBox::DenseVector<Field> FVector;

    PLinOpt::Tricounter mkn(PLinOpt::LRP2MM(L,R,P));
    const size_t& m(std::get<0>(mkn)), k(std::get<1>(mkn)), n(std::get<2>(mkn));

        // =============================================
        // Random inputs
    FVector va(FF,L.rowdim()), ua(FF,L.coldim());
    FVector vb(FF,R.rowdim()), ub(FF,R.coldim());
    FVector wc(FF,P.rowdim()), vc(FF,P.coldim());

    if ( (ua.size()!=(m*k)) || (ub.size()!=(k*n)) || (wc.size()!=(m*n)) ) {
        std::cerr << "# \033[1;31m****** ERROR, outer dimension mismatch: "
                  << L.coldim() << ':' << m << 'x' << k << ' '
                  << R.coldim() << ':' << k << 'x' << n << ' '
                  << P.rowdim() << ':' << m << 'x' << n
                  << " ******\033[0m"
                  << std::endl;
        return 3;
    }

    Givaro::GivRandom generator;
    Givaro::Integer::seeding(generator.seed());
    for(auto &iter:ua) RandomElt(iter, FF, generator, bitsize);
    for(auto &iter:ub) RandomElt(iter, FF, generator, bitsize);

        // =============================================
        // Compute matrix product via the HM algorithm
    L.apply(va,ua);
    R.apply(vb,ub);

    for(size_t i=0; i<vc.size(); ++i) FF.mul(vc[i],va[i],vb[i]);

    P.apply(wc,vc);

        // =============================================
        // Compute the matrix product directly
    LinBox::DenseMatrix<Field> Ma(FF,m,k), Mb(FF,k,n), Delta(FF,m,n);

        // row-major vectorization
    for(size_t i=0; i<ua.size(); ++i) Ma.setEntry(i/k,i%k,ua[i]);
    for(size_t i=0; i<ub.size(); ++i) Mb.setEntry(i/n,i%n,ub[i]);
    for(size_t i=0; i<wc.size(); ++i) Delta.setEntry(i/n,i%n,wc[i]);

    LinBox::MatrixDomain<Field> BMD(FF);
    LinBox::DenseMatrix<Field> Mc(FF,m,n);
    BMD.mul(Mc,Ma,Mb); // Direct matrix multiplication

    BMD.subin(Delta, Mc);

        // =============================================
        // Both computations should agree
    if (BMD.isZero (Delta)) {
        FF.write(std::clog <<"# \033[1;32mSUCCESS: correct "
                 << m << 'x' << k << 'x' << n
                 << " Matrix-Multiplication \033[0m") << std::endl;
    } else {
        FF.write(std::cerr << "# \033[1;31m****** ERROR, not a "
                 << m << 'x' << k << 'x' << n
                 << " MM algorithm******\033[0m") << std::endl;

        ua.write(std::clog << "Ua:=", FileFormat::Maple ) << ';' << std::endl;
        va.write(std::clog << "Va:=", FileFormat::Maple ) << ';' << std::endl;
        ub.write(std::clog << "Ub:=", FileFormat::Maple ) << ';' << std::endl;
        vb.write(std::clog << "Vb:=", FileFormat::Maple ) << ';' << std::endl;
        vc.write(std::clog << "Vc:=", FileFormat::Maple ) << ';' << std::endl;
        wc.write(std::clog << "wc:=", FileFormat::Maple ) << ';' << std::endl;

        Ma.write(std::clog << "Ma:=", FileFormat::Maple) << ';' << std::endl;
        Mb.write(std::clog << "Mb:=", FileFormat::Maple) << ';' << std::endl;
            // Correct value
        Mc.write(std::clog << "Mc:=", FileFormat::Maple) << ';' << std::endl;
            // Difference with computed value
        Delta.write(std::clog << "Df:=", FileFormat::Maple) << ';' << std::endl;
        BMD.addin(Delta, Mc);
            // Computed value
        Delta.write(std::clog << "Rc:=", FileFormat::Maple) << ';' << std::endl;

        return 1;
    }
    return 0;
}

// ===============================================================


} // End of namespace PLinOpt
// ============================================
