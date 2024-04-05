// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet
// ==========================================================================

/****************************************************************
 * PLinOpt Library, In-place inline implementations
 ****************************************************************/

#include "plinopt_inplace.h"


// ===============================================================
// Directed selection of accumulation i/o variables
Matrix::Row::const_iterator orientindex(const size_t preci,
                                        const Matrix::Row& L,
                                        const bool oriented) {
    if (! oriented) return L.begin();

            // Try to reuse previous variable
    auto nexti(std::find_if(L.begin(),L.end(),
                            [preci](const auto&a) { return a.first == preci; } )
               );

//         // But always prefer a variable with coefficient One
//     if ( (nexti == L.end()) || (! isOne(nexti->second)) ) {
//         // Otherwise, prefer a variable with coefficient One
    if (nexti == L.end()) {
        std::vector<Matrix::Row::const_iterator> vnext;
        for(auto iter=L.begin(); iter!=L.end(); ++iter) {
            if (isOne(iter->second)) vnext.push_back(iter);
        }
        if (vnext.size()>0) {
#ifdef RANDOM_TIES
            std::shuffle (vnext.begin(), vnext.end(),
                          std::default_random_engine(Givaro::BaseTimer::seed()));
#endif
            nexti = vnext.front();
        }
    }

#ifdef VERBATIM_PARSING
    if ( (nexti != L.end()) && (nexti != L.begin()) ) {
        std::clog << "# in " << L << ", prefered " << nexti->first
                  << " to " << L.begin()->first << std::endl;
    }
#endif

    return (nexti != L.end()) ? nexti : L.begin() ;
}
// ===============================================================


// ===============================================================
// Random selection of accumulation i/o variables
Matrix::Row::const_iterator nextindex(const size_t preci, const Matrix::Row& L,
                                      const bool oriented) {
#ifdef RANDOM_TIES
    if (oriented) return orientindex(preci, L, oriented);
    std::vector<Matrix::Row::const_iterator> vnext;
    for(auto iter=L.begin(); iter!=L.end(); ++iter) vnext.push_back(iter);
    std::shuffle (vnext.begin(), vnext.end(),
                  std::default_random_engine(Givaro::BaseTimer::seed()));
    return vnext.front();
#else
    return orientindex(preci, L, oriented);
#endif
}
// ===============================================================


// ===============================================================
// In-place program realizing a bilinear function
Tricounter BiLinearAlgorithm(std::ostream& out,
                             const Matrix& A, const Matrix& B,
                             const Matrix& T, const bool oriented) {

    const QRat& QQ = T.field();

    Tricounter opcount; // 0:ADD, 1:SCA, 2:MUL
    size_t preci(A.coldim()), precj(B.coldim()), preck(T.coldim());

        // Buffering clauses, to enable optimizing them out
    std::ostringstream saout, sbout, scout;
    std::map<std::pair<size_t,Givaro::Rational>,std::ostringstream>
        maaout, mabout, macout;

    const size_t m(A.rowdim());
    for(size_t l=0; l<m; ++l) {
//         std::clog << "# BEGIN parsing of: AXPY[" << l
//                   << "], astr: " << rmnl(saout.str())
//                   << ", bstr: " <<  rmnl(sbout.str()) << std::endl;

            /******************
             * Left hand side *
             ******************/
        const auto aiter { nextindex(preci, A[l], oriented) }; // choose var
        const size_t i(aiter->first);
        if ( (i == preci) && (l>0) && (aiter->second == A.getEntry(l-1,i)) ) {
            if (! QQ.isOne(aiter->second)) {
#ifdef VERBATIM_PARSING
                std::clog << "# Optimized out: scalar (div;mul), "
                          << 'a' << i << " /* " << VALPAR(aiter->second)
                          << std::endl;
#endif
                --std::get<1>(opcount);
                saout.clear(); saout.str(std::string());
            }
        } else {
            SCA(saout, 'a',i,'*',aiter->second, opcount); // scale choosen var
        }
        const bool aSCA(saout.tellp() != std::streampos(0));

        for(auto iter = A[l].begin(); iter != A[l].end(); ++iter) {
            if (iter != aiter) {
                const auto asearch(maaout.find(*iter));
                if ( (i == preci) && (l>0) && (! aSCA) &&
                     ( asearch != maaout.end() ) ) {
#ifdef VERBATIM_PARSING
                    std::clog << "# Optimized out: (add;sub), "
                              << 'a' << i << " +- "
                              << VALPAR(iter->second) << 'a' << iter->first
                              << std::endl;
#endif
                    maaout.erase(asearch);
                    --std::get<0>(opcount);
                } else {
                    ADD(saout, 'a',i,'+',iter->second, iter->first, opcount);
                }
            }
        }

            // Print remaining LHS clauses
        for(const auto& [key,value]: maaout) {
            out << value.str() << std::flush;
        }
        maaout.clear();
        out << saout.str() << std::flush;
        saout.clear(); saout.str(std::string());

            /*******************
             * Right hand side *
             *******************/
        const auto bjter( nextindex(precj, B[l], oriented) );
        const size_t j(bjter->first);
        if ( (j == precj) && (l>0) && (bjter->second == B.getEntry(l-1,j)) ) {
            if (! QQ.isOne(bjter->second)) {
#ifdef VERBATIM_PARSING
                std::clog << "# Optimized out, scalar (div;mul): "
                          << 'b' << j << " /* " << VALPAR(bjter->second)
                          << std::endl;
#endif
                --std::get<1>(opcount);
                sbout.clear(); sbout.str(std::string());
            }
        } else {
            SCA(sbout, 'b',j,'*',bjter->second, opcount);
        }
        const bool bSCA(sbout.tellp() != std::streampos(0));

        for(auto jter = B[l].begin(); jter != B[l].end(); ++jter) {
            if (jter != bjter) {
                const auto bsearch(mabout.find(*jter));
                if ( (j == precj) && (l>0) && (! bSCA) &&
                     ( bsearch != mabout.end() ) ) {
#ifdef VERBATIM_PARSING
                    std::clog << "# Optimized out: (add;sub), "
                              << 'b' << j << " +- "
                              << VALPAR(jter->second) << 'b' << jter->first
                              << std::endl;
#endif
                    mabout.erase(bsearch);
                    --std::get<0>(opcount);
                } else {
                    ADD(sbout, 'b',j,'+',jter->second, jter->first, opcount);
                }
            }
        }

            // Print remaining RHS clauses
        for(const auto& [key,value]: mabout) {
            out << value.str() << std::flush;
        }
        mabout.clear();
        out << sbout.str() << std::flush;
        sbout.clear(); sbout.str(std::string());


            /***********
             * Product *
             ***********/
        const auto ckter( nextindex(preck, T[l], oriented) );
        const size_t k(ckter->first);

        if (notAbsOne(QQ,ckter->second)) {
            if ( (k == preck) && (l>0)
                 && (ckter->second == T.getEntry(l-1,k)) ) {
                if (! QQ.isOne(ckter->second)) {
#ifdef VERBATIM_PARSING
                    std::clog << "# Optimized out, scalar (div;mul): "
                              << 'c' << k << " /* " << VALPAR(ckter->second)
                              << std::endl;
#endif
                    --std::get<1>(opcount);
                    scout.clear(); scout.str(std::string());
               }
            } else {
                SCA(scout, 'c', k, '/', ckter->second, opcount);
            }
        }
        const bool cSCA(scout.tellp() != std::streampos(0));

        for(auto kter = T[l].begin(); kter != T[l].end(); ++kter) {
            if (kter != ckter) {
                auto cadd(*kter); cadd.second *= ckter->second;
                const auto csearch(macout.find(cadd));
                if ( (k == preck) && (l>0) && (! cSCA) &&
                     ( csearch != macout.end()) ) {
#ifdef VERBATIM_PARSING
                    std::clog << "# Optimized out: (add;sub), "
                              << 'c' << kter->first << " +- "
                              << VALPAR(kter->second) << 'c' << k
                              << std::endl;
#endif
                    macout.erase(csearch);
                    --std::get<0>(opcount);
                } else {
                    ADD(scout, 'c', kter->first, MONEOP('-',ckter->second),
                        kter->second, k, opcount);
                }
            }
        }

            // Print remaining Product clauses
        for(const auto& [key,value]: macout) {
            out << value.str() << std::flush;
        }
        macout.clear();
        out << scout.str() << std::flush;
        scout.clear(); scout.str(std::string());



            /****************************************
             * the recursive (multiplicative) calls *
             ****************************************/
        MUL(out, 'c',k,MONEOP('+',ckter->second),'a',i,'b',j, opcount);



            /********************
             * RESTORE: Product *
             ********************/
        for(auto kter = T[l].begin(); kter != T[l].end(); ++kter) {
            if (kter != ckter) {
                auto cadd(*kter); cadd.second *= ckter->second;
                ADD(macout[cadd], 'c', kter->first, MONEOP('+',ckter->second),
                    kter->second, k, opcount);
            }
        }
        if (notAbsOne(QQ,ckter->second)) {
            SCA(scout, 'c', ckter->first, '*', ckter->second, opcount);
        }

            /****************************
             * RESTORE: Right hand side *
             ****************************/
        for(auto jter = B[l].begin(); jter != B[l].end(); ++jter) {
            if (jter != bjter) {
                ADD(mabout[*jter], 'b',j,'-',jter->second, jter->first, opcount);
            }
        }
        SCA(sbout, 'b',j,'/',bjter->second, opcount);

            /***************************
             * RESTORE: Left hand side *
             ***************************/
        for(auto iter = A[l].begin(); iter != A[l].end(); ++iter) {
            if (iter != aiter) {
                ADD(maaout[*iter], 'a',i,'-',iter->second, iter->first, opcount);
            }
        }
        SCA(saout, 'a',i,'/',aiter->second, opcount);

        preci=i; precj=j; preck=k;
//         std::clog << "# END of parsing of: AXPY[" << l << ']'<<  std::endl;
    }

        // Print last Product clauses
    for(const auto& [key,value]: macout) { out << value.str() << std::flush; }
    out << scout.str() << std::flush;
        // Print last RHS clauses
    for(const auto& [key,value]: mabout) { out << value.str() << std::flush; }
    out << sbout.str() << std::flush;
        // Print last LHS clauses
    for(const auto& [key,value]: maaout) { out << value.str() << std::flush; }
    out << saout.str() << std::flush;

    return opcount;
}
// ===============================================================


// ===============================================================
// Duplicate intermediate products
void DoubleExpand(Matrix& AA, Matrix& BB, Matrix& TT,
                  const Matrix& A, const Matrix& B, const Matrix& T) {

    AA.resize(A.rowdim()<<1,A.coldim());
    BB.resize(B.rowdim()<<1,B.coldim());
    TT.resize(T.rowdim()<<1,T.coldim()+1);

        /* Inflate Left */
    for(size_t l=0; l<A.rowdim();++l) {
        for(auto iter(A[l].begin()); iter != A[l].end(); ++iter) {
            AA.setEntry(l<<1,iter->first,iter->second);
            AA.setEntry(1+(l<<1),iter->first,iter->second);
        }
    }

        /* Inflate Right */
    for(size_t l=0; l<B.rowdim();++l) {
        for(auto iter(B[l].begin()); iter != B[l].end(); ++iter) {
            BB.setEntry(l<<1,iter->first,iter->second);
            BB.setEntry(1+(l<<1),iter->first,iter->second);
        }
    }

        /* Inflate & Interleave Product */
    for(size_t l=0; l<T.rowdim();++l) {
        for(auto iter(T[l].begin()); iter != T[l].end(); ++iter) {
            TT.setEntry(l<<1,iter->first,iter->second);
            TT.setEntry(1+(l<<1),1+(iter->first),iter->second);
        }
    }

#ifdef VERBATIM_PARSING
    AA.write(std::clog << "A:=",FileFormat::Maple) << ';' << std::endl;
    BB.write(std::clog << "B:=",FileFormat::Maple) << ';' << std::endl;
    TT.write(std::clog << "T:=",FileFormat::Maple) << ';' << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

}
// ===============================================================


// ===============================================================
// Searching the space of in-place bilinear programs
Tricounter SearchBiLinearAlgorithm(std::ostream& out,
                                   const Matrix& A, const Matrix& B,
                                   const Matrix& T, size_t randomloops) {

    if ( (A.rowdim() != B.rowdim()) || (A.rowdim() != T.rowdim()) ) {
        std::cerr << "Incorrect dimensions :" << std::endl;
		std::cerr << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;
		std::cerr << "B is " << B.rowdim() << " by " << B.coldim() << std::endl;
		std::cerr << "C is " << T.coldim() << " by " << T.rowdim() << std::endl;
        return Tricounter{0,0,0};
    }

    const QRat& QQ = T.field();

    Givaro::Timer elapsed;
    std::ostringstream sout, matout;
    Tricounter nbops(BiLinearAlgorithm(sout, A, B, T, true));
    std::string res(sout.str());
    std::clog << "# Oriented number of operations: " << nbops << std::endl;
    omp_lock_t writelock; omp_init_lock(&writelock);

#pragma omp parallel for shared(A,B,T,res,nbops)
    for(size_t i=0; i<randomloops; ++i) {

            // =============================================
            // random permutation of the rows
        LinBox::Permutation<QRat> P(QQ,A.rowdim());
        Givaro::GivRandom generator;
        P.random(generator.seed());

        Matrix pA(QQ,A.rowdim(), A.coldim()),
            pB(QQ,B.rowdim(), B.coldim()),
            pT(QQ,T.rowdim(), T.coldim()),
            TpT(QQ,T.coldim(), T.rowdim());
        permuteRows(pA,P,A,QQ);
        permuteRows(pB,P,B,QQ);
        permuteRows(pT,P,T,QQ);


            // =============================================
            // random coherent negation of the rows
        for(size_t i=0; i<pA.rowdim(); ++i) {
            const bool negA( generator.brand() ), negB( generator.brand() );
            if (negA) negRow(pA, i, QQ);
            if (negB) negRow(pB, i, QQ);
            if (negA != negB) negRow(pT, i, QQ);
        }

        std::ostringstream lout, oout;

            // =============================================
            // trying first a random selection of variables
        Tricounter lops(BiLinearAlgorithm(lout, pA, pB, pT,false));

        omp_set_lock(&writelock);
        if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
             ( (std::get<0>(lops)==std::get<0>(nbops))
               && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
            nbops = lops;
            res = lout.str();
            std::clog << "# Found algorithm[" << i << "], operations: " << lops << std::endl;
        }
        omp_unset_lock(&writelock);

            // =============================================
            // then trying a directed selection of variables
        lops = BiLinearAlgorithm(oout, pA, pB, pT, true);

        omp_set_lock(&writelock);
        if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
             ( (std::get<0>(lops)==std::get<0>(nbops))
               && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
            nbops = lops;
            res = oout.str();
            std::clog << "# Found oriented [" << i << "], operations: " << lops << std::endl;
        }
        omp_unset_lock(&writelock);
    }

        // Print the chosen algorithm
    out << res << std::flush;
    return nbops;
}
// ===============================================================
