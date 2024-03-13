// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, B. Grenet
// ==========================================================================

/****************************************************************
 * PLinOpt Library, In-place inline implementations
 ****************************************************************/

#include "plinopt_inplace.h"

// ===============================================================
// Tools for atomic operation (ADD or SCA)
struct Atom {
    char _var;
    size_t _src;
    char _ope;
    Givaro::Rational _val;
    size_t _des;

    Atom(char v, size_t s, char o, const Givaro::Rational& r, size_t d=0)
            : _var(v),_src(s),_ope(o),_val(r),_des(d) {}

    friend std::ostream& operator<<(std::ostream& out, const Atom& p) {
        const bool bsca(isSca(p._ope));

        if (bsca && isOne(p._val)) return out;


        if (p._ope == ' ') return out << "# row computed in "
                                      << p._var << p._src;


        out << p._var << p._src << ":="
            << p._var << p._src << ' ';

        if (bsca) {
            out << p._ope << ' ' << VALPAR(p._val);
        } else {
            const auto uval(sign(p._val)<0?-p._val:p._val);
            const auto uope(sign(p._val)<0?SWAPOP(p._ope):p._ope);
            out << uope << ' ';
            if (! isOne(uval)) {
                out << VALPAR(uval) << '*';
            }
            out << p._var << p._des;
        }
        return out << ';';
    }

    bool isinv(const Atom& p) const {
        bool binv(true);
        binv &= (this->_var == p._var);
        binv &= (this->_src == p._src);
        if (isAdd(this->_ope)) {
            if (this->_ope==p._ope) {
                binv &= isZero(this->_val+p._val);
            } else {
                binv &= isZero(this->_val-p._val);
            }
        } else {
            binv &= isOne(this->_val * p._val);
        }
        return binv &= (this->_des == p._des);
    }
};

auto isScaOne {[](const Atom& p){ return (isSca(p._ope) && isOne(p._val));} };

Tricounter complexity(const Program_t& p) {
    Tricounter nops;
    for(const auto& iter: p) {
        if (isAdd(iter._ope)) ++std::get<0>(nops);
        if (isSca(iter._ope)) ++std::get<1>(nops);
        if (iter._ope == ' ') ++std::get<2>(nops);
    }
    return nops;
}
std::ostream& operator<< (std::ostream& out, const Program_t& p) {
    for(const auto& iter: p) {
        out << iter << std::endl;
    }
    return out;
}


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
    for(auto iter=L.begin(); iter!=L.end(); ++iter) {
        vnext.push_back(iter);
    }
    std::shuffle (vnext.begin(), vnext.end(),
                  std::default_random_engine(Givaro::BaseTimer::seed()));
    return vnext.front();
#else
    return orientindex(preci, L, oriented);
#endif
}
// ===============================================================



// ===============================================================
// Removes the first atom followed by its inverse
bool simplify(Program_t& Program) {
    size_t cnt(0);
    for(auto iter = Program.begin(); iter != Program.end(); ++iter) {
        ++ cnt;
        auto next(iter);
        for(++next; next != Program.end(); ++next) {
            if (next->isinv(*iter)) {
#ifdef VERBATIM_PARSING
                std::clog << "# Removing[" << cnt << "] " << *iter << " with " << *next << std::endl;
#endif
                Program.erase(next);
                Program.erase(iter);

                return true;

            }
            if (next->_src == iter->_src) {
                if ( (next->_ope == ' ') || isSca(next->_ope) ) {
                    break;
                }
            }
        }
    }
    return false;
}

// ===============================================================





// ===============================================================
// In-place program realizing a linear function
Program_t& LinearAlgorithm(Program_t& Program, const Matrix& A,
                           const char variable,
                           const bool transposed, const bool oriented) {

    const QRat& QQ = A.field();
    Tricounter opcount; // 0:ADD, 1:SCA, 2:MUL
    const size_t m(A.rowdim());
    size_t preci(A.coldim());

    for(size_t l=0; l<m; ++l) {
        const auto aiter { nextindex(preci, A[l], oriented) }; // choose var
        const size_t i(aiter->first);

            // scale choosen var
        Program.emplace_back(variable,i,(transposed?'/':'*'),aiter->second);

			// Add the rest of the row
       for(auto iter = A[l].begin(); iter != A[l].end(); ++iter) {
            if (iter != aiter) {
                if (transposed) {
                    Program.emplace_back(variable, iter->first,
                                         MONEOP('-',aiter->second),
                                         iter->second, i);
                } else {
                    Program.emplace_back(variable, i,
                                         '+',
                                         iter->second, iter->first);
                }

            }
        }

        Program.emplace_back(variable,i,' ',QQ.zero); // placeholder for stop

			// Sub the rest of the row
        for(auto iter = A[l].begin(); iter != A[l].end(); ++iter) {
            if (iter != aiter) {
                if (transposed) {
                    Program.emplace_back(variable, iter->first,
                                         MONEOP('+',aiter->second),
                                         iter->second, i);
                } else {
                    Program.emplace_back(variable, i,
                                         '-',
                                         iter->second, iter->first);
                }
            }
        }

            // Un-scale choosen var
        Program.emplace_back(variable,i,(transposed?'*':'/'),aiter->second);

        preci = i;
    }

// std::clog << "# P setup : " << Program.size() << std::endl;
// std::clog <<  Program << std::endl;

        // Removing noops
    Program.erase(std::remove_if(Program.begin(),
                                 Program.end(),
                                 isScaOne),
                  Program.end());

// std::clog << "# P pruned: " << Program.size() << std::endl;
// std::clog <<  Program << std::endl;

        // Removes atoms followed by their inverses
    if (! transposed) {
        bool simp; do {
            simp = simplify(Program);
// std::clog << "# P simplf: " << Program.size() << std::endl;
// std::clog <<  Program << std::endl;
        } while(simp);
    }

    return Program;
}
// ===============================================================



// ===============================================================
// Searching the space of in-place linear programs
Tricounter SearchLinearAlgorithm(Program_t& Program, const Matrix& A,
                                 const char variable, size_t randomloops,
                                 const bool transposed) {
    const QRat& QQ = A.field();

    Givaro::Timer elapsed;
    std::ostringstream sout, matout;

    LinearAlgorithm(Program, A, variable, transposed, true);
    Tricounter nbops { complexity(Program) };


    std::string res(sout.str());
    std::clog << "# Oriented number of operations for " << variable
              << ": " << nbops << std::endl;

#pragma omp parallel for shared(A,Program,nbops)
    for(size_t i=0; i<randomloops; ++i) {

            // =============================================
            // random permutation of the rows
        LinBox::Permutation<QRat> P(QQ,A.rowdim());
        Givaro::GivRandom generator;
        P.random(generator.seed());

        Matrix pA(QQ,A.rowdim(), A.coldim());
        permuteRows(pA,P,A,QQ);

        Program_t lProgram; LinearAlgorithm(lProgram, A, variable, transposed);
        Tricounter lops { complexity(lProgram) };

        if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
             ( (std::get<0>(lops)==std::get<0>(nbops))
               && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
            nbops = lops;
            Program = lProgram;
            std::clog << "# Found algorithm[" << i << "] for " << variable
                      << ", operations: " << lops << std::endl;
        }
    }

    return nbops;
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

        if ( (!QQ.isOne(ckter->second)) && (!QQ.isMOne(ckter->second))) {
            if ( (k == preck) && (l>0) && (ckter->second == T.getEntry(l-1,k)) ) {
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
        if ( (!QQ.isOne(ckter->second)) && (!QQ.isMOne(ckter->second))) {
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

        if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
             ( (std::get<0>(lops)==std::get<0>(nbops))
               && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
            nbops = lops;
            res = lout.str();
            std::clog << "# Found algorithm[" << i << "], operations: " << lops << std::endl;
        }

            // =============================================
            // then trying a directed selection of variables
        lops = BiLinearAlgorithm(oout, pA, pB, pT, true);

        if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
             ( (std::get<0>(lops)==std::get<0>(nbops))
               && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
            nbops = lops;
            res = oout.str();
            std::clog << "# Found oriented [" << i << "], operations: " << lops << std::endl;
        }
    }

        // Print the chosen algorithm
    out << res << std::flush;
    return nbops;
}
// ===============================================================




// ===============================================================
// Enables checking a matrix multiplication with Maple
#ifdef INPLACE_CHECKER

    // Setup auxilliary variables
void InitializeVariables(const char L, const size_t m,
                         const char H, const size_t n,
                         const char F, const size_t s) {
    for(size_t h=0; h<m; ++h)
        std::clog << L << h << ":=L[" << (h+1) << "];";
    std::clog << std::endl;
    for(size_t h=0; h<n; ++h)
        std::clog << H << h << ":=H[" << (h+1) << "];";
    std::clog << std::endl;
    for(size_t h=0; h<s; ++h)
        std::clog << F << h << ":=F[" << (h+1) << "];";
    std::clog << std::endl;
    std::clog << std::string(30,'#') << std::endl;
}

    // Collect result of program
void CollectVariables(const char L, const size_t m,
                      const char H, const size_t n,
                      const char F, const size_t s) {
    std::clog << std::string(30,'#') << std::endl;
    for(size_t h=0; h<s; ++h)
        std::clog << "R[" << (h+1) << "]:=simplify(" << F << h << ",symbolic);";
    std::clog << std::endl;
    for(size_t h=0; h<n; ++h)
        std::clog << H << h << "-H[" << (h+1) << "],";
    std::clog << "0;" << std::endl;
    for(size_t h=0; h<m; ++h)
        std::clog << L << h << "-L[" << (h+1) << "],";
    std::clog << "0;" << std::endl;
    std::clog << std::string(30,'#') << std::endl;
}


    // Compare program with a matrix multiplication
void CheckMatrixMultiplication(const Matrix& A, const Matrix& B,
                               const Matrix& C) {
    Tricounter mkn(LRP2MM(A,B,C));
    const size_t& n(std::get<0>(mkn)), t(std::get<1>(mkn)), m(std::get<2>(mkn));
    std::clog <<"# code-checking for "
              << m << 'x' << t << 'x' << n
              << " Matrix-Multiplication" << std::endl;
    std::clog << '<';
    for(size_t i=0; i<m; ++i) {
        if (i!=0) std::clog << ',' << std::endl;
        std::clog << '<';
        for(size_t j=0; j<n; ++j) {
            if (j!=0) std::clog << '|';
            for(size_t k=0; k<t; ++k) {
                if (k!=0) std::clog << '+';
                std::clog << "L[" << (i*t+k+1) << "]*H[" << (k*n+j+1) << ']';
            }
            std::clog << " + F[" << (i*n+j+1) << "] - R[" << (i*n+j+1) << ']';
        }
        std::clog << '>';
    }
    std::clog << ">;" << std::endl;
}

#endif
// ===============================================================
