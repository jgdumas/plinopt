// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Sparsification  inline implementations
 ****************************************************************/

#include "plinopt_sparsify.h"


    // Find most frequent non-zero Rational between begin and end
    // element is suppoded to be the second in a pair
template <typename Fwd>
Givaro::Rational most_frequent_nonzero(const Fwd& begin, const Fwd& end)
{
    std::map<Givaro::Rational, int> count;

    for (auto it = begin; it != end; ++it) ++count[it->second];

    return std::max_element(
        count.begin(), count.end(),
        [](const auto& a, const auto& b) { return a.second < b.second;})->first;
}

	// Replace row i of A, by v
template<typename Vector>
Matrix& setRow(Matrix& A, size_t i, const Vector& v, const QRat& QQ) {
    A[i].resize(0);
    for(size_t j=0; j<v.size(); ++j) {
        if (! QQ.isZero(v[j])) A.setEntry(i,j,v[j]);
    }
    return A;
}

	// copy dense matrix M into sparse matrix A
Matrix& dense2sparse(Matrix& A, const DenseMatrix& M, const QRat& QQ) {
    A.resize(M.rowdim(), M.coldim());
    for(size_t i=0; i<A.rowdim(); ++i) {
        setRow(A,i,M[i],QQ);
    }
    return A;
}

	// copy sparse matrix B into sparse matrix A
Matrix& sparse2sparse(Matrix& A, const Matrix& B) {
    A.resize(B.rowdim(), B.coldim());
    std::copy(B.rowBegin(), B.rowEnd(), A.rowBegin());
    return A;
}

	// v is augmented by i, -i, 1/i and -1/i
template<typename Vector>
void augment(Vector& v, const Givaro::Integer& i, const QRat& QQ) {
    const Givaro::Rational r(i);
    if (std::find(v.begin(), v.end(), r) == v.end()) {
        v.push_back(r);
        v.push_back(-r);
        Givaro::Rational tmp;
        QQ.inv(tmp, r);
        v.push_back(tmp);
        QQ.negin(tmp);
        v.push_back(tmp);
    }
}

	// Computes the rank of A
size_t& rank(size_t& r, const Matrix& A) {
    static QRat QQ;
    static LinBox::GaussDomain<QRat> GD(QQ);
    Matrix copyA(QQ); sparse2sparse(copyA, A);
    return GD.rankInPlace(r, copyA);
}



	// Computes the transposed inverse of A
Matrix& inverseTranspose(Matrix& TI, const Matrix& A) {
    assert(A.rowdim() == A.coldim());
    const size_t n(A.rowdim());

    static QRat QQ;
    static LinBox::GaussDomain<QRat> GD(QQ);

    TI.resize(n,n);
    Matrix U(QQ); sparse2sparse(U, A);

    Givaro::Rational Det;
    size_t Rank;
    Matrix L(QQ, n, n);
    LinBox::Permutation<QRat> Q(QQ,n);
    LinBox::Permutation<QRat> P(QQ,n);

    GD.QLUPin(Rank, Det, Q, L, U, P, n, n );

    QVector x(QQ,n), w(QQ,n), Iv(QQ, n);
    for(size_t i=0; i<Iv.size(); ++i) Iv[i] = QQ.zero;

    for(int i=0; i<n; ++i) {
        Iv[i] = QQ.one;
        GD.solve(x, w, Rank, Q, L, U, P, Iv);
        Iv[i] = QQ.zero;
        setRow(TI, i, x, QQ);
    }

    return TI;
}

    // Prints out density profile of M
    // returns total density
std::ostream& densityProfile(std::ostream& out, size_t& ss, const Matrix& M) {
    ss = 0;
    for(auto it=M.rowBegin();it!=M.rowEnd();++it) {
        ss += it->size();
        out << it->size() << ' ';
    }
    return out << '=' << ss;
}



// ============================================================
// Testing whether w produces a sparsest row of TM
//    previous sparsity is in weight
//    w has to be independent of Cand
//    if sparsity is the same then test w's sparsity
//    If w is better, then replaces the line #num of LCoB
std::pair<int,int>& testLinComb(std::pair<int,int>& weight,
                                Matrix& LCoB, Matrix& Cand,
                                const size_t num,
                                const QArray& w, const Matrix& TM) {
    static QRat QQ;
    static QArray v(TM.coldim());
    static auto zeroTest = [](const auto& e) { return isZero(e);};

        // Check independency
    setRow(Cand,num,w,QQ);
    size_t r; rank(r, Cand);

    if (r>num) {
        TM.applyTranspose(v,w); // transpose of M.apply(v,w)

            // Count number of produced zeroes
        int hwl = std::count_if(v.begin(), v.end(), zeroTest );
        int vwl = std::count_if(w.begin(), w.end(), zeroTest );

            // Find the best one so far
        if ((hwl > weight.first) ||
            ( (hwl==weight.first) && (vwl>weight.second) ) ) {
#ifdef VERBATIM_PARSING
            std::clog << "# Found(" << num << "): " << w << ' ' << hwl
                      << " >= " << weight.first << std::endl;
#endif
            weight.first = hwl;
            weight.second = vwl;
                // Register best linear combination so far
            setRow(LCoB,num,w,QQ);
        }
    }
    return weight;
}


// ============================================================
// Sparsifying a Matrix TM
//   uses a limited number of coefficients for the linear comb.
//   (that is at most: maxnumcoeff)
Matrix& Sparsifier(Matrix& TCoB, Matrix& TM, const size_t maxnumcoeff) {
    static QRat QQ;

    assert(TCoB.rowdim() == TCoB.coldim());
    assert(TCoB.coldim() == TM.rowdim());
    const size_t n(TCoB.rowdim());

        // ========================================
        // Read Matrix of Linear Transformation
    Matrix LCoB(QQ,n,n);

        // ========================================
        // Try 0, 1, -1 and some coefficients in M
    std::vector<Givaro::Rational> Coeffs{0, 1, -1};
    for(auto row=TM.rowBegin(); row != TM.rowEnd(); ++row) {
        for(auto it=row->begin(); it != row->end(); ++it) {
            augment(Coeffs, it->second.nume(), QQ);
            augment(Coeffs, it->second.deno(), QQ);
         }
    }

        // ========================================
        // reduce to at most maxnumcoeff
    if (Coeffs.size()>maxnumcoeff) Coeffs.resize(maxnumcoeff);

        // ========================================
        // Try first 4 rows of TCoB, one at a time
    const size_t firstcolumns(std::min(size_t(4u),LCoB.rowdim()));

    for(size_t num=0; num<firstcolumns; ++num) {
        QArray v(TM.coldim());
        Matrix A(QQ); sparse2sparse(A, LCoB);
        std::pair<int,int> weight{-1,-1}; // Best Hamming weight so far

            // ========================================
            // Each row has 4 coefficients
        for(size_t i=0; i<Coeffs.size(); ++i) {
        for(size_t j=0; j<Coeffs.size(); ++j) {
        for(size_t k=0; k<Coeffs.size(); ++k) {
        for(size_t l=0; l<Coeffs.size(); ++l) {
                // Try linear combination
            QArray w{Coeffs[i],Coeffs[j],Coeffs[k],Coeffs[l]};
            w.resize(TM.rowdim());

            testLinComb(weight, LCoB, A, num, w, TM);
        }}}}
    }

        // Set bottom-right to identity
    for(size_t num=firstcolumns; num<LCoB.rowdim(); ++num) {
        LCoB.setEntry(num,num,QQ.one);
    }

        // ========================================
        // Supplementary change of basis found
        // Apply it:
        //   - to the sparsified matrix
        //   - to previous change of Basis
    static LinBox::MatrixDomain<QRat> BMD(QQ);
    DenseMatrix TR(QQ,TM.rowdim(), TM.coldim()),
        TS(QQ, TCoB.rowdim(),TCoB.coldim());
    BMD.mul(TR,LCoB,TM);   // Direct matrix multiplication
    BMD.mul(TS,LCoB,TCoB); // Direct matrix multiplication

        // Need sparse representation
    dense2sparse(TM,TR,QQ);
    dense2sparse(TCoB,TS,QQ);

    return TCoB;
}


// ============================================================
// Factoring out some coefficients:
//    from most frequent in a row of TM, towards a column of TCoB
Matrix& FactorDiagonals(Matrix& TCoB, Matrix& TM) {
    static QRat QQ;
    for(size_t i=0; i<TM.rowdim(); ++i) {
        if (TM[i].size()>0) {
            Givaro::Rational r(most_frequent_nonzero(TM[i].begin(),
                                                     TM[i].end() ) );
            if (! QQ.isOne(r)) {
                for(auto& iter : TM[i]) iter.second /= r;   // scale TM
                for(auto& iter : TCoB[i]) iter.second /= r; // scale TCoB
            }
        }
    }

    return TCoB;
}




// ============================================================
// Alternating sparsification and column factoring
//   starting with only 3 coefficents (thus -1,0,1) ...
//   ... increasing this number by increment ...
//   ... until threshold.
size_t SparseFactor(Matrix& TICoB, Matrix& TM,
                    const size_t start, const size_t increment,
                    const size_t threshold) {

        // ============================================================
        // Prints and computes density profile of TM
    size_t s2;
    densityProfile(std::clog << "# Columns profile: ", s2, TM) << std::endl;

#ifdef DEBUG
    static QRat QQ;
    static LinBox::MatrixDomain<QRat> BMD(QQ);
    static DenseMatrix TR(QQ,TM.rowdim(),TM.coldim());
    static Matrix CoB(QQ,TICoB.rowdim(), TICoB.coldim());
    static Matrix TCoB(QQ,CoB.coldim(), CoB.rowdim());
    inverseTranspose(CoB, TICoB); // CoB  == TICoB^{-T}
    Transpose(TCoB, CoB);         // TCoB == TICoB^{-1}
    BMD.mul(TR,TCoB,TM);          // TM == TICoB . TR
#endif

    size_t numcoeffs(start);
    size_t ss(s2);

        // ============================================================
        // Main loop, alternating sparsification and column factoring
    do {
        ss = s2;
        Sparsifier(TICoB, TM, numcoeffs);
        FactorDiagonals(TICoB, TM);
        densityProfile(std::clog << "# Density profile: ", s2, TM) << std::endl;
        if (numcoeffs<threshold) numcoeffs += increment;

#ifdef DEBUG
            // Check that a factorization R=M.CoB is preserved
            //              via TM = TICoB.TR
        static LinBox::MatrixDomain<QRat> BMD(QQ);
        consistency(std::clog, TM, TICoB, TR) << std::endl;
#endif
    } while ( s2 < ss );

    return s2;
}


// ============================================================
// Consistency check of M == R.C
template<typename AMatrix>
std::ostream& consistency(std::ostream& out, const AMatrix& M,
                          const Matrix& R, const DenseMatrix& C) {
    static QRat QQ;
    static LinBox::MatrixDomain<QRat> BMD(QQ);

    assert( M.rowdim() == R.rowdim() );
    assert( M.coldim() == C.coldim() );
    assert( R.coldim() == C.rowdim() );

    DenseMatrix A(QQ, R.rowdim(), C.coldim() );

    BMD.mul(A,R,C);
    BMD.subin(A,M);

    if (BMD.isZero (A))
        out <<"# \033[1;32mOK: consistent factorization!\033[0m";
    else{
        std::cerr << "# \033[1;31m****** ERROR inconsistency ******\033[0m"
                  << std::endl;
        M.write(out,FileFormat::Maple) << std::endl;
        out << " != " << std::endl;
        R.write(out,FileFormat::Maple) << std::endl;
        out << " * " << std::endl;
        C.write(out,FileFormat::Maple) << std::endl;
        out << std::string(30,'#') << std::endl;
    }

    return out;
}
