// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Sparsification  inline implementations
 ****************************************************************/

#include "plinopt_sparsify.h"


auto zeroTest = [](const auto& e) { return isZero(e);};
auto sizeSup = [](const auto& a, const auto& b) { return a.size() > b.size(); };
auto secondInf = [](const auto& a, const auto& b) { return a.second < b.second;};

    // Find most frequent non-zero Rational between begin and end
    // element is suppoded to be the second in a pair
template <typename Fwd>
inline Givaro::Rational most_frequent_nonzero(const Fwd& begin, const Fwd& end)
{
    std::map<Givaro::Rational, int> count;
    for (auto it = begin; it != end; ++it) ++count[it->second];
    return std::max_element(count.begin(), count.end(), secondInf)->first;
}

	// v is augmented by i, -i, 1/i and -1/i
template<typename Vector>
inline Vector& augment(Vector& v, const Givaro::Integer& i, const QRat& QQ) {
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
    return v;
}

	// Computes the rank of A
inline size_t& rank(size_t& r, const Matrix& A) {
    static QRat QQ;
    static LinBox::GaussDomain<QRat> GD(QQ);
    Matrix copyA(QQ); sparse2sparse(copyA, A);
    return GD.rankInPlace(r, copyA);
}



	// Computes the transposed inverse of A
inline Matrix& inverseTranspose(Matrix& TI, const Matrix& A) {
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

    for(size_t i=0; i<n; ++i) {
        Iv[i] = QQ.one;
        GD.solve(x, w, Rank, Q, L, U, P, Iv);
        Iv[i] = QQ.zero;
        setRow(TI, i, x, QQ);
    }

    return TI;
}


	// Build block diagonal matrix from vector of blocks
inline Matrix& diagonalMatrix(Matrix& M, const std::vector<Matrix>& V) {
    M.resize(0,0);
    for(const auto& mat: V) {
        const size_t m(M.rowdim()), n(M.coldim());
        M.resize(m+mat.rowdim(),n+mat.coldim());
        for( auto indices = mat.IndexedBegin();
             (indices != mat.IndexedEnd()) ; ++indices ) {
            M.setEntry(m+indices.rowIndex(),
                       n+indices.colIndex(),
                       indices.value());
        }
    }
    return M;
}

    // Append columns by blocks of columns
inline Matrix& augmentedMatrix(Matrix& M, const std::vector<Matrix>& V) {
    const size_t m(M.rowdim());
    M.resize(m,0);
    for(const auto& mat: V) {
        const size_t n(M.coldim());
        M.resize(m,n+mat.coldim());
        for( auto indices = mat.IndexedBegin();
             (indices != mat.IndexedEnd()) ; ++indices ) {
            assert(m == mat.rowdim());
            M.setEntry(indices.rowIndex(),
                       n+indices.colIndex(),
                       indices.value());
        }
    }
    return M;
}

    // Cut matrix by blocks of columns
inline std::vector<Matrix>& separateColumnBlocks(
    std::vector<Matrix>& V, const Matrix& A, const size_t blocksize) {
    static QRat QQ;
    const size_t numlargeblocks(A.coldim()/blocksize);
    const size_t lastblock(A.coldim()-numlargeblocks*blocksize);
    const bool hassmallblock(lastblock>0);

    for(size_t i=0; i<numlargeblocks; ++i)
        V.emplace_back(QQ,A.rowdim(),blocksize);
    if (hassmallblock)
        V.emplace_back(QQ,A.rowdim(),lastblock);

    for( auto indices = A.IndexedBegin(); (indices != A.IndexedEnd()) ;
         ++indices ) {
        V[indices.colIndex()/blocksize].setEntry(
            indices.rowIndex(),
            indices.colIndex() % blocksize,
            indices.value());
    }

    return V;
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
inline std::pair<int,int>& testLinComb(
    std::pair<int,int>& weight, Matrix& LCoB, Matrix& Cand,
    const size_t num,const QArray& w, const Matrix& TM) {

    static QRat QQ;
    static QArray v(TM.coldim()); v.resize(TM.coldim());

        // Check independency
    setRow(Cand,num,w,QQ);
    size_t r; rank(r, Cand);

    if (r>num) {
        TM.applyTranspose(v,w); // transpose of M.apply(v,w)

            // Count number of produced zeroes
        int rlHw = std::count_if(v.begin(), v.end(), zeroTest );
        int clHw = std::count_if(w.begin(), w.end(), zeroTest );

            // Find the best one so far
        if ((rlHw > weight.first) ||
            ( (rlHw==weight.first) && (clHw>weight.second) ) ) {
#ifdef VERBATIM_PARSING
            std::clog << "# Found(" << num << "): " << w << ' ' << rlHw
                      << " >= " << weight.first << std::endl;
#endif
            weight.first = rlHw;
            weight.second = clHw;
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
    static LinBox::GaussDomain<QRat> GD(QQ);

    const size_t n(TCoB.rowdim());
    assert(n == TCoB.rowdim());
    assert(n == TM.rowdim());
        // ========================================
        // Read Matrix of Linear Transformation
    Matrix LCoB(QQ,n,n);
    int cnHw(-1),rnHw(-1);

        // ========================================
        // Start by computing a nullspace vector:
        //   of the first denser rows
    if (TM.rowdim()>1) {
        Matrix N(QQ,TM.coldim(),TM.rowdim()); Transpose(N, TM);
        std::sort(N.rowBegin(), N.rowEnd(), sizeSup);	// get denser rows
        size_t r; while( (N.rowdim()>0) && (rank(r,N) == N.coldim())) {
            N.resize(N.rowdim()-1,N.coldim());          // select rank-1
        }
        if (N.rowdim()>0) {
            Matrix x(QQ, n, n), Tx(QQ,n,n);
            GD.nullspacebasisin(x, N);						// nullspace vector
            for(size_t i=0; i<n; ++i) {
                const auto& value( x.refEntry(i,0) );
                if (! QQ.isZero(value)) LCoB.setEntry(0,i,value);
            }
            QArray v(TM.coldim());
            TM.applyTranspose(v, LCoB[0]);					// result Hamming weight
            cnHw = LCoB[0].size();
            rnHw = std::count_if(v.begin(), v.end(), zeroTest );

#ifdef VERBATIM_PARSING
            std::clog << "# nullspace vector (" << cnHw << "): " << LCoB[0] << std::endl;
            std::clog << "# reduces residuum (" << rnHw << "): " << v << std::endl;
#endif
        }
    }

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
    const size_t numlargeblocks(TM.rowdim()>>2);
    const size_t lastblock(TM.rowdim()-(numlargeblocks<<2));
    const size_t numblocks(lastblock?numlargeblocks+1:numlargeblocks);
    const size_t multiple(numblocks<<2);
    QArray w(multiple);

    Matrix A(QQ);

    for(size_t block=0; block < numblocks; ++block) {
        w.resize(0); w.resize(multiple);
        const size_t firstcolumns(std::min(size_t(4u),LCoB.rowdim()-(block<<2)));

        for(size_t num=0; num<firstcolumns; ++num) {

            sparse2sparse(A, LCoB);
            std::pair<int,int> weight{-1,-1}; // Best Hamming weight so far
            if (num ==0) {
                weight.first=rnHw;
                weight.second=cnHw;
            }

                // ========================================
                // Each row has 4 coefficients
            for(size_t i=0; i<Coeffs.size(); ++i) {
            for(size_t j=0; j<Coeffs.size(); ++j) {
            for(size_t k=0; k<Coeffs.size(); ++k) {
            for(size_t l=0; l<Coeffs.size(); ++l) {
                    // Try linear combination

                w.resize(multiple);
                w[0+(block<<2)] = Coeffs[i];
                w[1+(block<<2)] = Coeffs[j];
                w[2+(block<<2)] = Coeffs[k];
                w[3+(block<<2)] = Coeffs[l];

                w.resize(TM.rowdim());
                testLinComb(weight, LCoB, A, num+(block<<2), w, TM);
            }}}}
        }
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
inline Matrix& FactorDiagonals(Matrix& TCoB, Matrix& TM) {
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
    DenseMatrix TR(QQ,TM.rowdim(),TM.coldim());
    Matrix CoB(QQ,TICoB.rowdim(), TICoB.coldim());
    Matrix TCoB(QQ,CoB.coldim(), CoB.rowdim());
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
        consistency(std::clog, TM, TICoB, TR) << std::endl;
#endif
    } while ( s2 < ss );

    return s2;
}



// ============================================================
// First:  FactorDiagonal
// Second: SparseFactor with default parameters
// Third:  SparseFactor with maxnumcoeff, 1, maxnumcoeff
Givaro::Timer& sparseAlternate(
    Givaro::Timer& chrono, Matrix& CoB, Matrix& Res, const Matrix& M,
    const FileFormat& matformat, const size_t maxnumcoeff) {
    chrono.start();

    static QRat QQ;
    Matrix TM(QQ,M.coldim(),M.rowdim()); Transpose(TM, M);

        // ============================================================
        // Initialize TICoB to identity
    Matrix TICoB(QQ,TM.rowdim(), TM.rowdim());
    for(size_t i=0; i<TICoB.coldim(); ++i) TICoB.setEntry(i,i,QQ.one);

        // ============================================================
        // Alternating sparsification and column factoring
        //    start by diagonals
    FactorDiagonals(TICoB, TM);
        //    default alternate to sparsify/factor simple things first
    SparseFactor(TICoB, TM);
        //    now try harder (with more potential combination coeffs)
    SparseFactor(TICoB, TM, maxnumcoeff, 1u, maxnumcoeff);

        // ============================================================
        // CoB = TICoB^{-T}, transposed inverse
        // Res = TM^T
    inverseTranspose(CoB, TICoB);
    size_t sc;
    densityProfile(std::clog << "# CoBasis profile: ", sc, CoB) << std::endl;

    Transpose(Res, TM);
    chrono.stop();

#ifdef VERBATIM_PARSING
        // Transposed Inverse change of basis to stdlog
    TICoB.write(std::clog, matformat)<< std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif
    return chrono;
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
        out <<"# \033[1;32mSUCCESS: consistent factorization!\033[0m";
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
