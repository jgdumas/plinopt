// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Sparsification  inline implementations
 ****************************************************************/

#include "plinopt_sparsify.h"


namespace PLinOpt {
// ===============================================================
auto zeroTest { [](const auto& e) { return isZero(e);} };
auto sizeSup { [](const auto& a, const auto& b) { return a.size() > b.size();} };
auto secondInf {[](const auto& a, const auto& b) { return a.second < b.second;}};

	// v is augmented by i, -i, 1/i and -1/i
template<typename Vector, typename _Field>
inline Vector& augment(Vector& v, const typename _Field::Element& i,
                       const _Field& FF) {
    using Element = typename _Field::Element;
    const Element& r(i);
    if (std::find(v.begin(), v.end(), r) == v.end()) {
        v.push_back(r);
        v.push_back(-r);
        Element tmp;
        FF.inv(tmp, r);
        v.push_back(tmp);
        FF.negin(tmp);
        v.push_back(tmp);
    }
    return v;
}

	// Computes the rank of A
template<typename _Mat>
inline size_t& rank(size_t& r, const _Mat& A) {
    using Field = typename _Mat::Field;
    const Field& FF(A.field());
    LinBox::GaussDomain<Field> GD(FF);
    _Mat copyA(FF); matrixCopy(copyA, A);
    return GD.rankInPlace(r, copyA);
}


	// Build block diagonal matrix from vector of blocks
template<typename _Mat>
inline _Mat& diagonalMatrix(_Mat& M, const std::vector<_Mat>& V) {
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
template<typename _FMat, typename _Mat2>
inline _FMat& augmentedMatrix(_FMat& M, const std::vector<_Mat2>& V) {
    using Field = typename _FMat::Field;
    const Field& F(M.field());
    const size_t m(M.rowdim());
    M.resize(m,0);
    for(const auto& mat: V) {
        const size_t n(M.coldim());
        M.resize(m,n+mat.coldim());
        for( auto indices = mat.IndexedBegin();
             (indices != mat.IndexedEnd()) ; ++indices ) {
            assert(m == mat.rowdim());
            if (! F.isZero(indices.value()))
                M.setEntry(indices.rowIndex(),
                           n+indices.colIndex(),
                           indices.value());
        }
    }
    return M;
}

    // Cut matrix by blocks of columns
template<typename _Mat>
inline std::vector<_Mat>& separateColumnBlocks(std::vector<_Mat>& V,
                                               const _Mat& A,
                                               const size_t blocksize) {
    using FMatrix = _Mat;
    using Field = typename FMatrix::Field;
    const Field& FF(A.field());

    const size_t numlargeblocks(A.coldim()/blocksize);
    const size_t lastblock(A.coldim()-numlargeblocks*blocksize);
    const bool hassmallblock(lastblock>0);

    for(size_t i=0; i<numlargeblocks; ++i)
        V.emplace_back(FF,A.rowdim(),blocksize);
    if (hassmallblock)
        V.emplace_back(FF,A.rowdim(),lastblock);

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
template<typename _Mat>
std::ostream& densityProfile(std::ostream& out, size_t& ss, const _Mat& M) {
    ss = 0;
    for(auto it=M.rowBegin();it!=M.rowEnd();++it) {
        ss += it->size();
        out << it->size() << ' ';
    }
    return out << '=' << ss;
}


// Show profile & test for consistency
// show 1 --> std::cout, show 0 --> noshow, show # --> std::clog
template<typename _Mat>
size_t profileConsistency(const FileFormat& matformat, Givaro::Timer& elapsed,
                          const _Mat& C, size_t sc, const std::string& fcode,
                          const _Mat& A, int showA,
                          const _Mat& B, int showB) {
    size_t sa,sb;
    densityProfile(std::clog << "# " << fcode
                   << " chgobase profile: \033[1;36m",
                   sb, B) << "\033[0m" << std::endl;
    if (showB != 0) {
        B.write( showB==1? std::cout : std::clog, matformat) << std::endl;
    }
    densityProfile(std::clog << "# " << fcode
                   << " residuum profile: \033[1;36m",
                   sa, A) << "\033[0m" << std::endl;
    if (showA != 0) {
        A.write( showA==1? std::cout : std::clog, matformat) << std::endl;
    }
    consistency(std::clog, C, A, B)
        << " \033[1;36m" << A.rowdim() << 'x' << A.coldim()
        << " by " << B.rowdim() << 'x' << B.coldim() << " with "
        << sa << " non-zeroes (" << sb << " alt.) instead of " << sc
        << "\033[0m:" << ' ' << elapsed << std::endl;

    return sa;
}



// ============================================================
// Testing whether w produces a sparsest row of TM
//    previous sparsity is in weight
//    w has to be independent of Cand
//    if sparsity is the same then test w's sparsity
//    If w is better, then replaces the line #num of LCoB
template<typename _Mat, typename _Array>
inline bool testLinComb(Pair<int>& weight, _Mat& LCoB, _Mat& Cand,
                        const size_t num, const _Array& w, const _Mat& TM) {
    _Array v(TM.coldim()); v.resize(TM.coldim());

        // Check independency
    setRow(Cand,num,w);
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
            setRow(LCoB,num,w);
            return true;
        }
    }
    return false;
}



// ============================================================
// Sparsifying a Matrix TM
//   uses a limited number of coefficients for the linear comb.
//   (that is at most: maxnumcoeff)
template<typename _Mat>
_Mat& Sparsifier(_Mat& TCoB, _Mat& TM, const size_t maxnumcoeff) {
    using FMatrix = _Mat;
    using Field = typename FMatrix::Field;
    using Element = typename Field::Element;
    using DenseFMatrix = LinBox::DenseMatrix<Field>;
    using FArray = std::vector<Element>;
    const Field& FF(TM.field());

    LinBox::GaussDomain<Field> GD(FF);

    const size_t n(TCoB.rowdim());
    assert(n == TCoB.rowdim());
    assert(n == TM.rowdim());
        // ========================================
        // Read Matrix of Linear Transformation
    FMatrix LCoB(FF,n,n);
    int cnHw(-1),rnHw(-1);

        // ========================================
        // Start by computing a nullspace vector:
        //   of the first denser rows
    if (TM.rowdim()>1) {
        FMatrix N(FF,TM.coldim(),TM.rowdim()); Transpose(N, TM);
        std::sort(N.rowBegin(), N.rowEnd(), sizeSup);	// get denser rows
        size_t r; while( (N.rowdim()>0) && (rank(r,N) == N.coldim())) {
            N.resize(N.rowdim()-1,N.coldim());          // select rank-1
        }
        if (N.rowdim()>0) {
            FMatrix x(FF, n, n), Tx(FF,n,n);
            GD.nullspacebasisin(x, N);						// nullspace vector
            for(size_t i=0; i<n; ++i) {
                const auto& value( x.refEntry(i,0) );
                if (! FF.isZero(value)) LCoB.setEntry(0,i,value);
            }
            FArray v(TM.coldim());
            TM.applyTranspose(v, LCoB[0]);					// result Hamming weight
            cnHw = LCoB[0].size();
            rnHw = std::count_if(v.begin(), v.end(), zeroTest );

#ifdef VERBATIM_PARSING
            std::clog << "# [SPRF] nullspace vector (" << cnHw << "): " << LCoB[0] << std::endl;
            std::clog << "# [SPRF] reduces residuum (" << rnHw << "): " << v << std::endl;
#endif
        }
    }

        // ========================================
        // Try 0, 1, -1 and some coefficients in M
    std::vector<Element> Coeffs{0, 1, -1};
    for(auto row=TM.rowBegin(); row != TM.rowEnd(); ++row) {
        for(auto it=row->begin(); it != row->end(); ++it) {
            augment(Coeffs, it->second, FF);
         }
    }

    for(size_t i=2; Coeffs.size() < maxnumcoeff; ++i) {
        augment(Coeffs, Element(i), FF);
    }
        // ========================================
        // reduce to at most maxnumcoeff
    if (Coeffs.size()>maxnumcoeff) Coeffs.resize(maxnumcoeff);
    std::clog << "# [SPRF] linear combination coefficients: " << Coeffs << std::endl;

        // ========================================
        // Try first 4 rows of TCoB, one at a time
    const size_t numlargeblocks(TM.rowdim()>>2);
    const size_t lastblock(TM.rowdim()-(numlargeblocks<<2));
    const size_t numblocks(lastblock?numlargeblocks+1:numlargeblocks);
    const size_t multiple(numblocks<<2);
    FArray w(multiple);

    FMatrix A(FF);

    for(size_t block=0; block < numblocks; ++block) {
        w.resize(0); w.resize(multiple);
        const size_t offsetblock(block<<2);
        const size_t firstcolumns(std::min(size_t(4u),
                                           LCoB.rowdim()-(offsetblock)));

        for(size_t num=0; num<firstcolumns; ++num) {
            matrixCopy(A, LCoB);
            Pair<int> weight{-1,-1}; // Best Hamming weight so far
            bool found((block == 0) && (num == 0));
            if (found) {
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
                w[0+(offsetblock)] = Coeffs[i];
                w[1+(offsetblock)] = Coeffs[j];
                w[2+(offsetblock)] = Coeffs[k];
                w[3+(offsetblock)] = Coeffs[l];

                w.resize(TM.rowdim());
                found |= testLinComb(weight, LCoB, A, num+(offsetblock), w, TM);
            }}}}

                // If not enough lin. comb. just add an indep. canonical one
            for(size_t p=0; ! found; ++p) {
                weight = {-1,-1};
                w.resize(0); w.resize(TM.rowdim());
                w[p]=1;
                found |= testLinComb(weight, LCoB, A, num+(offsetblock), w, TM);
#ifdef VERBATIM_PARSING
                if (found) std::clog << "# [SPRF] Using canonical " << p << std::endl;
#endif
            }

        }
    }

        // ========================================
        // Supplementary change of basis found
        // Apply it:
        //   - to the sparsified matrix
        //   - to previous change of Basis
    LinBox::MatrixDomain<Field> BMD(FF);
    DenseFMatrix TR(FF,TM.rowdim(), TM.coldim()),
        TS(FF, TCoB.rowdim(),TCoB.coldim());
    BMD.mul(TR,LCoB,TM);   // Direct matrix multiplication
    BMD.mul(TS,LCoB,TCoB); // Direct matrix multiplication

        // Need sparse representation
    dense2sparse(TM,TR);
    dense2sparse(TCoB,TS);

    return TCoB;
}



// ============================================================
// Factoring out some coefficients:
//    from most frequent in a row of TM, towards a column of TCoB
template<typename _Mat>
inline _Mat& FactorDiagonals(_Mat& TCoB, _Mat& TM) {
    using Field = typename _Mat::Field;
    using Element = typename Field::Element;
    const Field& FF(TM.field());
    for(size_t i=0; i<TM.rowdim(); ++i) {
        if (TM[i].size()>0) {
            std::map<Element, int> count;
            auto begin(TM[i].begin());
            auto end(TM[i].end());
            for (auto it = begin; it != end; ++it) ++count[it->second];
            Element r(std::max_element(count.begin(), count.end(),
                                       secondInf)->first);
            if (! FF.isOne(r)) {
                for(auto& iter : TM[i]) FF.divin(iter.second,r);   // scale TM
                for(auto& iter : TCoB[i]) FF.divin(iter.second,r); // scale TCoB
            }
        }
    }

    return TCoB;
}



	// Computes the inverse of A
template<typename _Mat1, typename _Mat2>
inline _Mat1& inverse(_Mat1& T, const _Mat2& A) {
    assert(A.rowdim() == A.coldim());
    const size_t n(A.rowdim());

    using Field = typename _Mat1::Field;
    using Element = typename Field::Element;
    using FVector = LinBox::DenseVector<Field>;
    const Field& FF(A.field());

    LinBox::GaussDomain<Field> GD(FF);

    T.resize(n,n);
    _Mat1 U(FF); matrixCopy(U, A);

    Element Det;
    size_t Rank;
    _Mat1 L(FF, n, n);
    LinBox::Permutation<Field> Q(FF,n);
    LinBox::Permutation<Field> P(FF,n);

    GD.QLUPin(Rank, Det, Q, L, U, P, n, n );

    FVector x(FF,n), w(FF,n), Iv(FF, n);
    for(size_t j=0; j<Iv.size(); ++j) FF.assign(Iv[j],FF.zero);

    for(size_t j=0; j<n; ++j) {
        FF.assign(Iv[j],FF.one);
        GD.solve(x, w, Rank, Q, L, U, P, Iv);
        FF.assign(Iv[j],FF.zero);
        updCol(T, j, x);
    }

    return T;
}
    // ============================================================
    // Compute R, s.t. A == TICoB . R
template<typename _Mat1, typename _Mat2, typename _DMat>
_DMat& applyInverse(_DMat& R, const _Mat1& TICoB, const _Mat2& A,
                    const LinBox::MatrixDomain<typename _Mat1::Field>& BMD) {
    using FMatrix = _Mat1;
    using Field = typename FMatrix::Field;
    const Field& FF(TICoB.field());
    FMatrix CoB(FF,TICoB.rowdim(), TICoB.coldim());
    FMatrix TCoB(FF,CoB.coldim(), CoB.rowdim());
    inverseTranspose(CoB, TICoB); // CoB  == TICoB^{-T}
    Transpose(TCoB, CoB);         // TCoB == TICoB^{-1}
    return BMD.mul(R,TCoB,A);     // A == TICoB . R
}

	// Computes the transposed inverse of A
template<typename _Mat1, typename _Mat2>
inline _Mat1& inverseTranspose(_Mat1& TI, const _Mat2& A) {
    assert(A.rowdim() == A.coldim());
    const size_t n(A.rowdim());

    using Field = typename _Mat1::Field;
    using Element = typename Field::Element;
    using FVector = LinBox::DenseVector<Field>;
    const Field& FF(A.field());

    LinBox::GaussDomain<Field> GD(FF);

    TI.resize(n,n);
    _Mat1 U(FF); matrixCopy(U, A);

    Element Det;
    size_t Rank;
    _Mat1 L(FF, n, n);
    LinBox::Permutation<Field> Q(FF,n);
    LinBox::Permutation<Field> P(FF,n);

    GD.QLUPin(Rank, Det, Q, L, U, P, n, n );

    FVector x(FF,n), w(FF,n), Iv(FF, n);
    for(size_t i=0; i<Iv.size(); ++i) Iv[i] = FF.zero;

    for(size_t i=0; i<n; ++i) {
        Iv[i] = FF.one;
        GD.solve(x, w, Rank, Q, L, U, P, Iv);
        Iv[i] = FF.zero;
        setRow(TI, i, x);
    }

    return TI;
}


// ============================================================
// Alternating sparsification and column factoring
//   starting with only 3 coefficents (thus -1,0,1) ...
//   ... increasing this number by increment ...
//   ... until threshold.
template<typename _Mat>
size_t SparseFactor(_Mat& TICoB, _Mat& TM,
                    const size_t start, const size_t increment,
                    const size_t threshold) {
        // ============================================================
        // Prints and computes density profile of TM
    size_t s2;
    densityProfile(std::clog << "# [SpFc] Columns profile: ", s2, TM) << std::endl;

#ifdef DEBUG
    using Field = typename _Mat::Field;
    using DenseFMatrix = LinBox::DenseMatrix<Field>;
    LinBox::MatrixDomain<Field> BMD(FF);
    DenseFMatrix TR(FF,TM.rowdim(),TM.coldim());
    applyInverse(TR, TICoB, TM, BMD);	// TM == TICoB . TR
#endif

    size_t numcoeffs(start);
    size_t ss(s2);

        // ============================================================
        // Main loop, alternating sparsification and column factoring
    do {
        ss = s2;
        Sparsifier(TICoB, TM, numcoeffs);
        FactorDiagonals(TICoB, TM);
        densityProfile(std::clog << "# [SpFc] Density profile: ", s2, TM) << std::endl;
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
// QLUP Gaussian elimination of A
//    A  <-- (QL)^{-1} . A
//    QL <-- Q.L
// Updates QL and A only if resulting A = UP is sparser
template<typename _Mat>
inline bool sparseLU(_Mat& QL, _Mat& A, const size_t sparsity) {
    using FMatrix = _Mat;
    using Field = typename FMatrix::Field;
    using Element = typename Field::Element;
    using DenseFMatrix = LinBox::DenseMatrix<Field>;
    const Field& FF(A.field());

    LinBox::GaussDomain<Field> GD(FF);
    LinBox::MatrixDomain<Field> BMD(FF);

    const size_t m(A.rowdim()), n(A.coldim());
#ifdef DEBUG
    DenseFMatrix TR(FF,m,n); any2dense(TR, A);
#endif

    Element Det;
    size_t Rank;
    FMatrix U(FF,m,n); matrixCopy(U, A);
    FMatrix L(FF,m,m);
    LinBox::Permutation<Field> Q(FF,m);
    LinBox::Permutation<Field> P(FF,n);
    GD.QLUPin(Rank, Det, Q, L, U, P, m, n );

    bool sparser(density(U)<sparsity);

    if (sparser) {
            // use it
        DenseFMatrix R(FF,U.rowdim(), U.coldim()), B(FF,U.rowdim(), U.coldim()),
            S(FF, L.rowdim(),L.coldim()), C(FF, L.rowdim(),L.coldim());
        matrixCopy(R, U);
        matrixCopy(S, L);
        P.applyLeft(B, R);
        Q.applyRight(C, S);

        dense2sparse(A, B);
        dense2sparse(QL, C);
    }

#ifdef DEBUG
            // Check that the factorization is preserved
            //              via TR = QL.A
        consistency(std::clog, TR, QL, A) << std::endl;
#endif
    return sparser;
}


// ============================================================
// QLUP Gaussian elimination of A with pre-application of TC
//    A  <--  UP = (QL)^{-1} . A
//    TC <--  (QL)^{-1} . TC
// Updates TC and A only if resulting A is sparser
template<typename _Mat>
inline bool sparseILU(_Mat& TC, _Mat& A, const size_t sparsity) {
    using FMatrix = _Mat;
    using Field = typename FMatrix::Field;
    using DenseFMatrix = LinBox::DenseMatrix<Field>;
    const Field& FF(A.field());
    LinBox::MatrixDomain<Field> BMD(FF);

    const size_t m(A.rowdim());
#ifdef DEBUG
    DenseFMatrix TR(FF,m,A.coldim());
    applyInverse(TR, TC, A, BMD);	// TR s.t., A == TC . TR
#endif

    FMatrix QL(FF,m,m); for(size_t i=0; i<m; ++i) QL.setEntry(i,i,FF.one);
    bool sparser( sparseLU(QL, A, sparsity) );
    if (sparser) {
        DenseFMatrix K(FF,m,m);
        applyInverse(K, QL, TC, BMD);
        dense2sparse(TC, K);
    }

#ifdef DEBUG
            // Check that the factorization is preserved
            //              via A = TC.TR
        consistency(std::clog, A, TC, TR) << std::endl;
#endif
    return sparser;
}
// ============================================================
// First:  FactorDiagonal
// Second: SparseFactor with default parameters
// Third:  SparseFactor with maxnumcoeff, 1, maxnumcoeff
template<typename _Mat>
Givaro::Timer& sparseAlternate(Givaro::Timer& chrono, _Mat& CoB, _Mat& Res,
                               const _Mat& M, const FileFormat& matformat,
                               const size_t maxnumcoeff) {
    using Field = typename _Mat::Field;
    const Field& F(M.field());

    chrono.start();

    const size_t m(M.rowdim()), n(M.coldim());
    _Mat TM(F,n,m); Transpose(TM, M);

        // ============================================================
        // Initialize TICoB to identity
    _Mat TICoB(F,n,n);
    for(size_t i=0; i<n; ++i) TICoB.setEntry(i,i,F.one);

        // ============================================================
        // Alternating sparsification and column factoring
        //    start by diagonals
    FactorDiagonals(TICoB, TM);
        //    Then use QLUP factorization, if result is sparser
    bool reduced = sparseILU(TICoB, TM, density(TM));
    if (reduced) {
        size_t sl,su;
        densityProfile(std::clog << "# [sALT] GaussLo profile: ", sl, TICoB) << std::endl;
        densityProfile(std::clog << "# [sALT] GaussUp profile: ", su, TM) << std::endl;
    }
        //    default alternate to sparsify/factor simple things first
    SparseFactor(TICoB, TM);
        //    now try harder (with more potential combination coeffs)
    SparseFactor(TICoB, TM, maxnumcoeff, 1u, maxnumcoeff);

        // ============================================================
        // CoB = TICoB^{-T}, transposed inverse
        // Res = TM^T
    inverseTranspose(CoB, TICoB);
    size_t sc;
    densityProfile(std::clog << "# [sALT] CoBasis profile: ", sc, CoB) << std::endl;

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
// Sparsifying and reducing coefficient diversity of a matrix
// by sparse QLUP elimination, followed by block sparsification
template<typename _Mat>
int blockSparsifier(Givaro::Timer& elapsed, _Mat& CoB, _Mat& Res,
                    const _Mat& M, const size_t blocksize,
                    const FileFormat& matformat, const size_t maxnumcoeff,
                    const bool initialElimination) {
    using FMatrix = _Mat;
    using Field = typename FMatrix::Field;
    using DenseFMatrix = LinBox::DenseMatrix<Field>;
    const Field& F(M.field());

        // ============================================================
        // Sparsify matrix as a whole
    if (blocksize <= 1) {
        sparseAlternate(elapsed, CoB, Res, M, matformat, maxnumcoeff);
    } else {
        Givaro::Timer chrono; chrono.start();
        const size_t m(M.rowdim()), n(M.coldim());
        FMatrix U(F,n,m), L(F,n,n);
        bool reduced(initialElimination);

        if (initialElimination) {
            Transpose(U, M);
                // ============================================================
                // Initialize L to identity
            for(size_t i=0; i<n; ++i) L.setEntry(i,i,F.one);
            reduced = sparseLU(L, U, density(U));

            size_t sl,su;
            densityProfile(std::clog << "# [bSpr] IGaussL profile: ", sl, L) << std::endl;
            densityProfile(std::clog << "# [bSpr] IGaussU profile: ", su, U) << std::endl;
        }

        FMatrix TU(F,m,n);
        const FMatrix& A( reduced? Transpose(TU,U) : M);

        chrono.stop();
        elapsed += chrono;

            // ============================================================
            // Deal with blocks of columns
        std::vector<FMatrix> vC, vR, vA;
        separateColumnBlocks(vA, A, blocksize);
        for(const auto& mat: vA) {
            vC.emplace_back(F,mat.coldim(), mat.coldim());
            vR.emplace_back(F,mat.rowdim(), mat.coldim());

            sparseAlternate(chrono, vC.back(), vR.back(),
                            mat, matformat, maxnumcoeff);

            elapsed += chrono;
#ifdef DEBUG
            std::clog << std::string(30,'#') << std::endl;
            consistency(std::clog, mat, vR.back(), vC.back()) << ' ' << chrono << std::endl;
#endif
        }

            // Build resulting matrices
        augmentedMatrix(Res, vR);

        if (reduced) {
            LinBox::MatrixDomain<Field> BMD(F);
            std::vector<FMatrix> vL;
            std::vector<DenseFMatrix> vB;
            separateColumnBlocks(vL, L, blocksize);
            for(size_t i=0; i<vL.size(); ++i) {
                vB.emplace_back(F,vL[i].rowdim(), vL[i].coldim());
                FMatrix TvC(F, vC[i].coldim(), vC[i].rowdim());
                Transpose(TvC, vC[i]);
                BMD.mul(vB.back(),vL[i],TvC);
            }
            FMatrix TCoB(F,CoB.coldim(), CoB.rowdim());
            augmentedMatrix(TCoB, vB);
            Transpose(CoB, TCoB);
        } else {
            diagonalMatrix(CoB, vC);
        }
    }

    return 0;
}

// ============================================================
// Decomposing the matrix into:
//   (Res=[identity,lower part])*(CoB=[upperpart])
//   with prescribed inner dimension
// precondition: upper part is full-rank
template<typename _Mat>
Tricounter backSolver(_Mat& CoB, _Mat& Res, const _Mat& iM) {
    using FMatrix = _Mat;
    using Field = typename FMatrix::Field;
    using Element = typename Field::Element;
    using DenseFMatrix = LinBox::DenseMatrix<Field>;
    using FVector = LinBox::DenseVector<Field>;


    const Field& FF(iM.field());

    const size_t r(iM.rowdim()), n(iM.coldim()), k(CoB.rowdim());
    const size_t s(r-k);
    assert( Res.rowdim() == r );  // Res is (r x k)
    assert( Res.coldim() == k );
    assert( CoB.coldim() == n );  // CoB is (k x n)

    FMatrix M(FF, iM.rowdim(), iM.coldim());
#ifdef RANDOM_TIES
    DenseFMatrix dM(FF, iM.rowdim(), iM.coldim());
    LinBox::Permutation<Field> S(FF,r);
    S.random(Givaro::BaseTimer::seed());
    S.applyRight(dM, iM);
    matrixCopy(M,dM);
#else
    matrixCopy(M,iM);
#endif

    LinBox::Permutation<Field> T(FF,r);

        // Select n independent rows
    size_t rk(0);
    for(size_t i=0; i<n; ++i) {
        for(size_t j=i;j<r;++j) {
            setRow(CoB, i, M, j);
            if (rank(rk,CoB) == (i+1)) {
                if (i != j) {
                    T.permute(i,j);
                    std::swap(M[i],M[j]);
                }
                break;
            }
        }
        if (rk != (i+1)) {
            CoB[i].resize(0);
            break;
        }
    }

        // Add up to k rows
    for(size_t j=rk; j<k; ++j) {
        setRow(CoB,j, M, j);
    }
        // Other rows to be solved for
    FMatrix A2(FF,s,n);
    for(size_t j=k; j<r; ++j)
        setRow(A2,j-k, M, j);

#ifdef VERBATIM_PARSING
    T.write(std::clog   << "# [bSol] Initial perm.: ") << std::endl;
    CoB.write(std::clog << "# [bSol] Full row rank: ",FileFormat::Pretty) << std::endl;
    A2.write(std::clog  << "# [bSol] Free profile : ",FileFormat::Pretty) << std::endl;
#endif

    FMatrix U(FF,n,k); Transpose(U, CoB); // U is (n x k)
    FMatrix B(FF,n,s); Transpose(B, A2);  // B is (n x s)

    Element Det;
    size_t Rank;
    FMatrix L(FF, n, n);
    LinBox::Permutation<Field> P(FF,n);
    LinBox::Permutation<Field> Q(FF,k);

        // Gaussian elimination of the CoB upper part
    LinBox::GaussDomain<Field> GD(FF);
    GD.QLUPin(Rank, Det, P, L, U, Q, n, k );


#ifdef VERBATIM_PARSING
    P.write(std::clog << "CoB P: ") << std::endl;
    L.write(std::clog << "CoB L: ", FileFormat::Pretty) << std::endl;
    U.write(std::clog << "CoB U: ", FileFormat::Pretty) << std::endl;
    Q.write(std::clog << "CoB Q: ") << std::endl;
#endif

        // Upper part is identity
    for(size_t i=0; i<k; ++i) Res.setEntry(i,i,FF.one);

        // Solving for each vector of the lower part
    DenseFMatrix TY(FF,s,n); Transpose(TY,B);
    FVector x(FF,k), w(FF,k), bi(FF,n);
    for(size_t i=0; i<s; ++i) {
        for(size_t j=0; j<n; ++j)
            bi[j] = TY[i][j];
        GD.solve(x, w, r, P, L, U, Q, bi);
        setRow(Res, k+i, x);
    }

    DenseFMatrix R(FF, r, k);
    T.solveRight(R, Res);
#ifdef RANDOM_TIES
    DenseFMatrix R2(FF, r, k);
    S.solveRight(R2, R);
    dense2sparse(Res, R2);
#else
    dense2sparse(Res, R);
#endif

    const auto rops {nonzeroes(Res) };
    return Tricounter { rops.first, rops.second, density(CoB) };
}

// ============================================================
// Consistency check of M == R.C
template<typename _Mat1, typename _Mat2, typename _DMat>
std::ostream& consistency(std::ostream& out, const _Mat1& M,
                          const _Mat2& R, const _DMat& C) {

    using Field = typename _Mat1::Field;
    using DenseFMatrix = LinBox::DenseMatrix<Field>;

    const Field& FF(M.field());

    LinBox::MatrixDomain<Field> BMD(FF);

    assert( M.rowdim() == R.rowdim() );
    assert( M.coldim() == C.coldim() );
    assert( R.coldim() == C.rowdim() );

    DenseFMatrix A(FF, R.rowdim(), C.coldim() );

    BMD.mul(A,R,C);
    BMD.subin(A,M);

    if (BMD.isZero (A))
        out <<"# \033[1;32mSUCCESS: consistent factorization!\033[0m";
    else{
        std::cerr << "# \033[1;31m****** ERROR inconsistency ******\033[0m"
                  << std::endl;
        DenseFMatrix dM(FF, M.rowdim(), M.coldim()); any2dense(dM,M);
        dM.write(out,FileFormat::Maple) << std::endl;
        out << " != " << std::endl;
        DenseFMatrix dR(FF, R.rowdim(), R.coldim()); any2dense(dR,R);
        dR.write(out,FileFormat::Maple) << std::endl;
        out << " * " << std::endl;
        C.write(out,FileFormat::Maple) << std::endl;
        out << std::string(30,'#') << std::endl;
    }

    return out;
}


// ============================================================
// Factorizing into an Alternate basis time a CoB

// Optimizing for non-zeroes of Alt, then non-ones, then nnz of CoB
auto tricOpCount {[](const auto& a, const auto& b) {
    return (std::get<0>(a) < std::get<0>(b))
        || ( (std::get<0>(a) == std::get<0>(b))
             && (std::get<1>(a) < std::get<1>(b)) )
        || ( (std::get<0>(a) == std::get<0>(b))
             && (std::get<1>(a) == std::get<1>(b))
             && (std::get<2>(a) < std::get<2>(b)) )
; } };


template<typename _Mat>
int Factorizer(_Mat& Alt, _Mat& CoB, const _Mat& M,
               const size_t randomloops, const size_t selectinnerdim,
               const bool progressreport) {
    using FMatrix = _Mat;
    using Field = typename FMatrix::Field;
    const Field& F(M.field());


#ifdef VERBATIM_PARSING
    M.write(std::clog,FileFormat::Pretty) << std::endl;
#endif
    const size_t innerdim(selectinnerdim == 0 ? M.coldim() : selectinnerdim);
    if ( (innerdim > M.rowdim()) ||
         (innerdim < M.coldim()) ) {
        std::cerr << "# \033[1;36mFail: inner dimension has to be between " << M.coldim() << " and " << M.rowdim() << ".\033[0m\n";
        return -1;
    }

    Alt.resize(M.rowdim(), innerdim);
    if (M.rowdim() == M.coldim()) {
        matrixCopy(CoB, M);
        for(size_t i(0); i<Alt.rowdim(); ++i) Alt.setEntry(i, i, F.one);
        std::clog << std::string(30,'#') << std::endl;
        std::clog <<"# \033[1;36mWARNING: identity factorization\033[0m\n";
        return 0;
    }

    sparse2sparse(Alt, M);
    CoB.resize(innerdim, M.coldim());
    for(size_t i(0); i<CoB.coldim(); ++i) CoB.setEntry(i, i, F.one);

    const auto sc(nonzeroes(M));

	// Start with M and Identity
    Tricounter nbops{ sc.first, sc.second, M.coldim()};

#pragma omp parallel for shared(Alt,CoB,M,F,nbops,innerdim,progressreport)
    for(size_t i=0; i<randomloops; ++i) {
        FMatrix lCoB(F, innerdim, M.coldim());
        FMatrix lAlt(F, M.rowdim(), innerdim);
        auto bSops{backSolver(lCoB, lAlt, M)};

#pragma omp critical
        {
        if (tricOpCount(bSops, nbops)) {
#ifdef VERBATIM_PARSING
            std::clog << "# Alt/CoB profile[" << i << "]: "
                      << bSops << std::endl;
#endif
            nbops = bSops;
            sparse2sparse(CoB, lCoB);
            sparse2sparse(Alt, lAlt);

            if (progressreport)
                std::clog << "# Found [" << i << "], R/CB profile: "
                          << bSops << std::endl;
        }
        }
    }


    return 0;


}

// ============================================================







} // End of namespace PLinOpt
// ============================================
