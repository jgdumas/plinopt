// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
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
                                      << p._var << p._src << VALPAR(p._val);


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
        if (isAdd(this->_ope) && isAdd(p._ope)) {
            if (this->_ope==p._ope) {
                binv &= isZero(this->_val+p._val);
            } else {
                binv &= isZero(this->_val-p._val);
            }
        } else if (isSca(this->_ope) && isSca(p._ope)) {
            binv &= isOne(this->_val / p._val);
        } else {
            // at least one of the two operations is not arithmetic
            return false;
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
    if ( (nexti == L.end()) || (! isOne(nexti->second)) ) {
//         // Only otherwise, prefer a variable with coefficient One
//     if (nexti == L.end()) {
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
bool simplify(Program_t& Program, const bool transposed) {
    size_t cnt(0);
    for(auto iter = Program.begin(); iter != Program.end(); ++iter) {
        ++ cnt; size_t nnt(0);

        if (iter->_ope == ' ') {
            continue;
        }

        auto next(iter);
        for(++next; next != Program.end(); ++next) {
            ++ nnt;
            if (next->isinv(*iter)) {
#ifdef VERBATIM_PARSING
                std::clog << "# Removing[" << cnt << "] " << *iter
                          << " with[" << nnt << "] " << *next << std::endl;

                for(auto rit=iter; rit != next; ++rit) {
                    std::clog << "# %% " << *rit << std::endl;
                }
                std::clog << "# %% " << *next << std::endl;
#endif

                Program.erase(next);
                Program.erase(iter);

                return true;

            }

				// Stop optimizing, if iter is reused from now, in certain cases
            if (transposed) {
                if ( ( (iter->_des == next->_des)
                       &&
                       ( (next->_ope == ' ') || isSca(next->_ope) ) )
                     ||
                     (iter->_des == next->_src)
                     ) {
                    break;
                }
            } else {
                if ( ( (iter->_src == next->_src)
                       &&
                       ( (next->_ope == ' ') || isSca(next->_ope) ) )
                     ||
                     ( (iter->_des == next->_src) && (next->_ope != ' ') )
                     ||
                     (iter->_src == next->_des)
                     ) {
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
Tricounter LinearAlgorithm(Program_t& Program, const Matrix& A,
                           const char variable,
                           const bool transposed, const bool oriented) {

    const QRat& QQ = A.field();
    Tricounter opcount; // 0:ADD, 1:SCA, 2:MUL
    const size_t m(A.rowdim());
    size_t preci(A.coldim());

    for(size_t l=0; l<m; ++l) {
        const auto& CurrentRow(A[l]);
            // choose var
        const auto aiter { nextindex(preci, CurrentRow, oriented) };
        const size_t i(aiter->first);

            // scale choosen var
        if (transposed) {
            if ( (!QQ.isOne(aiter->second)) && (!QQ.isMOne(aiter->second))) {
                Program.emplace_back(variable,i,'/',aiter->second);
            }
        } else {
            Program.emplace_back(variable,i,'*',aiter->second);
        }


			// Add the rest of the row
       for(auto iter = CurrentRow.begin(); iter != CurrentRow.end(); ++iter) {
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

           // placeholder for barrier:
           //   (AXPY will have to performed at this point)
       Program.emplace_back(variable,i,' ', aiter->second);

			// Sub the rest of the row
        for(auto iter = CurrentRow.begin(); iter != CurrentRow.end(); ++iter) {
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
        if (transposed) {
            if ( (!QQ.isOne(aiter->second)) && (!QQ.isMOne(aiter->second))) {
                Program.emplace_back(variable,i,'*',aiter->second);
            }
        } else {
            Program.emplace_back(variable,i,'/',aiter->second);
        }


        if (CurrentRow.size()>1) preci = i;
    }

        // Removing noops
    Program.erase(std::remove_if(Program.begin(), Program.end(), isScaOne),
                  Program.end());

        // Removes atoms followed by their inverses, one at a time
    bool simp; do {
        simp = simplify(Program, transposed);
    } while(simp);

    return complexity(Program);
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

    Tricounter nbops{ LinearAlgorithm(Program, A, variable, transposed, true) };


    std::string res(sout.str());
#ifdef VERBATIM_PARSING
    std::clog << "# Oriented number of operations for " << variable
              << ": " << nbops << std::endl;
#endif

#pragma omp parallel for shared(A,Program,nbops)
    for(size_t i=0; i<randomloops; ++i) {

            // =============================================
            // random permutation of the rows
        LinBox::Permutation<QRat> P(QQ,A.rowdim());
        Givaro::GivRandom generator;
        P.random(generator.seed());

            // =============================================
            // Apply permutation to A
        Matrix pA(QQ,A.rowdim(), A.coldim());
        permuteRows(pA,P,A,QQ);

        Program_t lProgram;
        Tricounter lops { LinearAlgorithm(lProgram, pA, variable, transposed) };

        if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
             ( (std::get<0>(lops)==std::get<0>(nbops))
               && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
            nbops = lops;
            Program = lProgram;
#ifdef VERBATIM_PARSING
            std::clog << "# Found algorithm[" << i << "] for " << variable
                      << ", operations: " << lops << std::endl;
#endif
        }

        lops = LinearAlgorithm(lProgram, pA, variable, transposed, true);

        if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
             ( (std::get<0>(lops)==std::get<0>(nbops))
               && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
            nbops = lops;
            Program = lProgram;
#ifdef VERBATIM_PARSING
            std::clog << "# Found oriented[" << i << "] for " << variable
                      << ", operations: " << lops << std::endl;
#endif
        }
    }

    return nbops;
}
// ===============================================================





// ===============================================================
// In-place program realizing a trilinear function
Tricounter TriLinearAlgorithm(std::ostream& out,
                              const Matrix& A, const Matrix& B,
                              const Matrix& T, const bool oriented=false) {

    const QRat& QQ = T.field();

    Tricounter opcount; // 0:ADD, 1:SCA, 2:MUL
    size_t preci(A.coldim()), precj(B.coldim()), preck(T.coldim());

        // Buffering clauses, to enable optimizing them out
    std::ostringstream saout, sbout, scout;
    std::map<std::pair<size_t,Givaro::Rational>,std::ostringstream>
        maaout, mabout, macout;

    const size_t m(A.rowdim());
    for(size_t l=0; l<m; ++l) {
        const auto& CurrentARow(A[l]), CurrentBRow(B[l]), CurrentTRow(T[l]);

            /******************
             * Left hand side *
             ******************/
            // choose var
        const auto aiter { nextindex(preci, CurrentARow, oriented) };
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

        for(auto iter=CurrentARow.begin(); iter != CurrentARow.end(); ++iter) {
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
        const auto bjter( nextindex(precj, CurrentBRow, oriented) );
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

        for(auto jter = CurrentBRow.begin(); jter != CurrentBRow.end(); ++jter) {
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
        const auto ckter( nextindex(preck, CurrentTRow, oriented) );
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

        for(auto kter = CurrentTRow.begin(); kter != CurrentTRow.end(); ++kter) {
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
        for(auto kter = CurrentTRow.begin(); kter != CurrentTRow.end(); ++kter) {
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
        for(auto jter = CurrentBRow.begin(); jter != CurrentBRow.end(); ++jter) {
            if (jter != bjter) {
                ADD(mabout[*jter], 'b',j,'-',jter->second, jter->first, opcount);
            }
        }
        SCA(sbout, 'b',j,'/',bjter->second, opcount);

            /***************************
             * RESTORE: Left hand side *
             ***************************/
        for(auto iter = CurrentARow.begin(); iter != CurrentARow.end(); ++iter) {
            if (iter != aiter) {
                ADD(maaout[*iter], 'a',i,'-',iter->second, iter->first, opcount);
            }
        }
        SCA(saout, 'a',i,'/',aiter->second, opcount);

        if (CurrentARow.size()>1) preci=i;
        if (CurrentBRow.size()>1) precj=j;
        if (CurrentTRow.size()>1) preck=k;
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


Tricounter& operator+=(Tricounter& l, const Tricounter& r) {
    std::get<0>(l) += std::get<0>(r);
    std::get<1>(l) += std::get<1>(r);
    std::get<2>(l) += std::get<2>(r);
    return l;
}




// ===============================================================
// In-place program realizing a trilinear function
Tricounter TriLinearProgram(std::ostream& out, const Matrix& A, const Matrix& B,
                           const Matrix& T, const bool oriented) {
    const QRat& QQ(T.field());

    Program_t aprog, bprog, cprog;


    Tricounter aops { LinearAlgorithm(aprog, A, 'a', false, oriented) };
    out << "# Found " << aops << " for a" << std::endl;

    Tricounter bops { LinearAlgorithm(bprog, B, 'b', false, oriented) };
    out << "# Found " << bops << " for b" << std::endl;

    Tricounter cops { LinearAlgorithm(cprog, T, 'c', true, oriented) };
    out << "# Found " << cops << " for c" << std::endl;

    Tricounter pty;
        // Synchronizing the programs
    auto aiter(aprog.begin()), bjter(bprog.begin()), ckter(cprog.begin()) ;
    for( ; ckter != cprog.end() ; ) {
            // Print ADD/SCA, left hand sides
        for( ; (aiter != aprog.end()) && (aiter->_ope != ' ') ; ++aiter) {
            out << *aiter << std::endl;
        }
            // Print ADD/SCA, right hand sides
        for( ; (bjter != bprog.end()) && (bjter->_ope != ' ') ; ++bjter) {
            out << *bjter << std::endl;
        }
            // Print ADD/SCA, products
        for( ; (ckter != cprog.end()) && (ckter->_ope != ' ') ; ++ckter) {
            out << *ckter << std::endl;
        }

            // Print (recursive) AXPY
        if ( (aiter != aprog.end()) &&  (bjter != bprog.end()) &&  (ckter != cprog.end()) ) {

            MUL(out, ckter->_var, ckter->_src, MONEOP('+',ckter->_val),
                aiter->_var, aiter->_src, bjter->_var, bjter->_src, pty);

            ++aiter; ++bjter; ++ckter;
        }
    }


        // Last reversions to recover initial state for inputs
    for( ; ckter != cprog.end() ; ++ckter) {
        out << *ckter << std::endl;
    }
    for( ; bjter != bprog.end() ; ++bjter) {
        out << *bjter << std::endl;
    }
    for( ; aiter != aprog.end() ; ++aiter) {
        out << *aiter << std::endl;
    }

        // Number of operations in the program
    Tricounter nops; nops += aops; nops += bops; nops += cops;

    std::get<2>(nops) /= 3; // MUL is counted in each of the three programs

    return nops;
}
// ===============================================================


// ===============================================================
// Searching the space of in-place trilinear programs
Tricounter SearchTriLinearAlgorithm(std::ostream& out,
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
    Tricounter nbops{ TriLinearProgram(sout, A, B, T, true) };
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

            // Apply permutation to A, B, T
        Matrix pA(QQ,A.rowdim(), A.coldim()), pB(QQ,B.rowdim(), B.coldim()),
            pT(QQ,T.rowdim(), T.coldim());
        permuteRows(pA,P,A,QQ); permuteRows(pB,P,B,QQ); permuteRows(pT,P,T,QQ);


            // =============================================
            // random coherent negation of the rows
        for(size_t i=0; i<pA.rowdim(); ++i) {
            const bool negA( generator.brand() ), negB( generator.brand() );
            if (negA) negRow(pA, i, QQ);
            if (negB) negRow(pB, i, QQ);
            if (negA != negB) negRow(pT, i, QQ);
        }

        std::ostringstream lout, sout;

            // =============================================
            // trying first a directed selection of variables
        Tricounter lops{ TriLinearProgram(lout, pA, pB, pT, true) };

        omp_set_lock(&writelock);
        if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
             ( (std::get<0>(lops)==std::get<0>(nbops))
               && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
            nbops = lops;
            res = lout.str();
            std::clog << "# Found oriented [" << i << "], operations: " << lops << std::endl;
        }
        omp_unset_lock(&writelock);

            // =============================================
            // then trying a random selection of variables
        lops = TriLinearProgram(sout, pA, pB, pT);

        omp_set_lock(&writelock);
        if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
             ( (std::get<0>(lops)==std::get<0>(nbops))
               && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
            nbops = lops;
            res = sout.str();
            std::clog << "# Found combined [" << i << "], operations: " << lops << std::endl;
        }
        omp_unset_lock(&writelock);
    }

        // Print the chosen algorithm
    out << res << std::flush;
    return nbops;

//    return  TriLinearProgram(out, A, B, T);



 }
// ===============================================================




// ===============================================================
// Enables checking a matrix multiplication with Maple
#ifdef INPLACE_CHECKER

    // Setup auxilliary variables
void InitializeVariable(const char L, const size_t m, const char M)  {
    for(size_t h=0; h<m; ++h)
        std::clog << L << h << ":=" << M << "[" << (h+1) << "];";
    std::clog << std::endl;
}

void InitializeVariables(const char L, const size_t m,
                         const char H, const size_t n,
                         const char F, const size_t s) {
    InitializeVariable(L, m, 'L');
    InitializeVariable(H, n, 'H');
    InitializeVariable(F, s, 'F');
    std::clog << std::string(30,'#') << std::endl;
}

    // Collect result of program
void CollectVariable(const char L, const size_t m, const char M) {
    for(size_t h=0; h<m; ++h)
        std::clog << L << h << "-" << M << "[" << (h+1) << "],";
    std::clog << "0;" << std::endl;
}

void CollectVariables(const char L, const size_t m,
                      const char H, const size_t n,
                      const char F, const size_t s) {
    std::clog << std::string(30,'#') << std::endl;
    for(size_t h=0; h<s; ++h)
        std::clog << "R[" << (h+1) << "]:=simplify(" << F << h << ",symbolic);";
    std::clog << std::endl;
    CollectVariable(H, n, 'H');
    CollectVariable(L, m, 'L');
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
