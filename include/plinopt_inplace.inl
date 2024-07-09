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
    long _des;

    Atom(char v, size_t s, char o, const Givaro::Rational& r, long d=-1)
            : _var(v),_src(s),_ope(o),_val(r),_des(d) {}

    Atom(const Atom& c)
            : _var(c._var),_src(c._src),_ope(c._ope),_val(c._val),_des(c._des){}

    Atom& operator= (const Atom& c) {
        this->_var = c._var;
        this->_src = c._src;
        this->_ope = c._ope;
        this->_val = c._val;
        this->_des = c._des;
        return *this;
    }

        // Raw print
    friend std::ostream& dprint(std::ostream& out, const Atom& p) {
        return out << p._var << p._src << ' ' << p._ope << ' '
                   << p._val << ' ' << p._var << p._des;
    }

        // print as a program (removes noops and barriers)
    friend std::ostream& operator<<(std::ostream& out, const Atom& p) {
        const bool bsca(isMulDiv(p._ope));
        QRat QQ;

        if (bsca && QQ.isOne(p._val)) return out;


        if (p._ope == ' ') {
#ifdef DEBUG
            std::clog << "# row computed in "
                      << p._var << p._src << VALPAR(p._val) << std::endl;
#endif
            if (QQ.isZero(p._val)) {
                return out << '0' << ';';
            } else {
                return out << p._var << p._src << ';';
            }
        }


        out << p._var << p._src << ":=";

        size_t nbmul(0);

        const auto uval(sign(p._val)<0?-p._val:p._val);

        if (bsca) {
            if (sign(p._val)<0) out << '-';
            printSCA(out, p._var, p._src, p._ope, uval, nbmul, QQ);
        } else {
            const auto uope(sign(p._val)<0?SWAPOP(p._ope):p._ope);
            out << p._var << p._src << uope;

            printmulorjustdiv(out, p._var, p._des, uval, nbmul, QQ);
        }
        return out << ';';
    }

        // Operation is on the same variable(s)
    bool sameops(const Atom& p) const {
        return (this->_var == p._var)
            && (this->_src == p._src)
            && (this->_des == p._des);
    }

        // Operation is a no-op
    bool isnoop() const {
        return (isAddSub(this->_ope) && isZero(this->_val))
            || (isMulDiv(this->_ope) && isOne(this->_val));
    }

        // Combine both operations, if compatible
    bool cumulate(const Atom& p) {
        if (this->sameops(p)) {
            if (isAddSub(this->_ope) && isAddSub(p._ope)) {
                if (this->_ope==p._ope)
                    this->_val += p._val;
                else
                    this->_val -= p._val;
                if (sign(this->_val) < 0) {
                    this->_ope = SWAPOP(this->_ope);
                    QRat QQ; QQ.negin(this->_val);
                }
                return true;
            } else if (isMulDiv(this->_ope) && isMulDiv(p._ope)) {
                if (this->_ope==p._ope)
                    this->_val *= p._val;
                else
                    this->_val /= p._val;
                QRat QQ;
                if (abs(this->_val) < QQ.one) {
                    this->_ope = INVOP(this->_ope);
                    QQ.invin(this->_val);
                }
                return true;
            }
        }
        return false;
    }

};


// Lambda Testing whether atom is '*' or '/' with a 1
auto isMulDivOne {[](const Atom& p){ return (isMulDiv(p._ope)
                                          && Givaro::isOne(p._val));} };


template<typename Ring>
Tricounter complexity(const Ring& R, const AProgram_t& p) {
    Tricounter nops;
    for(const auto& iter: p) {
        if (isAddSub(iter._ope)) {
            ++std::get<0>(nops);
            if (notAbsOne(R,iter._val)) ++std::get<1>(nops);
        }
        if (isMulDiv(iter._ope)) ++std::get<1>(nops);
        if (iter._ope == ' ') ++std::get<2>(nops);
    }
    return nops;
}

// Printing a runable program
std::ostream& operator<< (std::ostream& out, const AProgram_t& p) {
    for(const auto& iter: p) if (iter._ope != ' ') out << iter << std::endl;
    return out;
}

// Raw print program (including barrier)
std::ostream& dprint (std::ostream& out, const AProgram_t& p) {
    for(const auto& iter: p)
        if (iter._ope != ' ')
            out << iter << std::endl;
        else
            dprint(out << "# barrier # ", iter) << std::endl;
    return out;
}

// Raw print program (with explicitly permuted output)
std::ostream& Pprint(std::ostream& out, const char c,
                     const AProgram_t& atomP,
                     const LinBox::Permutation<QRat>& P) {
//     P.write(std::clog << "# Permutation: ") << std::endl;
    size_t numop(0);
    for(const auto& iter: atomP) {
        if (iter._ope == ' ') {
            out << c << P[numop++] << ":=";
        }
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
bool simplify(AProgram_t& Program, const bool transposed) {
    size_t cnt(0);
    for(auto iter = Program.begin(); iter != Program.end(); ++iter) {
        ++ cnt; size_t nnt(0);

        if (iter->_ope == ' ') {
            continue;
        }

        auto next(iter);
        for(++next; next != Program.end(); ++next) {
            ++ nnt;
            if (next->sameops(*iter)) { // same variables
#ifdef VERBATIM_PARSING
                std::clog << "# Found[" << cnt << "] " << *iter
                          << " with[" << cnt+nnt << "] " << *next << std::endl;
                for(auto rit=iter; rit != next; ++rit) {
                    dprint(std::clog << "## %f% ", *rit) << std::endl;
                }
                dprint( std::clog << "## %f% ", *next) << std::endl;
#endif

                Atom c(*iter);
                bool bpos = c.cumulate(*next);

                if (bpos) { // same variables and compatible operations
#ifdef VERBATIM_PARSING
                    std::clog << "## --> replaced by: " << c << std::endl;
#endif
                    Program.erase(next);
                    if (c.isnoop()) { // operations are inverse of one another
                        Program.erase(iter);
                    } else {
                        *iter = c;
                    }
                    return true;
                }
            }

                // Whether or not to Stop optimizing
                // --> if iter is reused from now, in certain cases
            bool needbreak(false);

                // STOP if: same source but different usage
            needbreak |=
                (iter->_src == next->_src)
                &&
                ( (next->_ope == ' ')
                  ||
                  (isAddSub(iter->_ope) && isMulDiv(next->_ope))
                  ||
                  (isMulDiv(iter->_ope) && isAddSub(next->_ope))
                  );

                // STOP if: previous rhs is now modified
                // except: barrier will modify output, but not inputs
            if (transposed)
                needbreak |= (iter->_des == (long)(next->_src));
            else
                needbreak |= ((iter->_des == (long)(next->_src)) && (next->_ope != ' '));

                // STOP if: previous lhs is now used
            needbreak |= ((long)(iter->_src) == next->_des);

            if (needbreak) break;
        }
    }
    return false;
}

// ===============================================================



// ===============================================================
// Push independent variables towards the end
// --> semantically equivalent if: AXPY is not encountered
// --> otherwise until next modification of same variable
// --> makes same variables closer, thus helps simplification
void pushvariables(AProgram_t& Program, size_t numout) {
        // For each variable
    for(size_t i=0; i<numout; ++i) {
        bool found(false);
        auto findex(Program.begin());

        for(auto iter = Program.begin(); iter != Program.end(); ++iter) {
            if (found) { // (findex_src == i) and (findex_ope is Add or Mul)
                bool rotate(false);
                if (isAddSub(findex->_ope)) {
                    if (findex->_des == (long)(iter->_src)) {
                            // f_des is modified, cannot push
                        found = false; continue;
                    } else if (iter->_src == i) {
                        if (findex->_des == iter->_des) {
                                // i_ope must be '+'
                            rotate=true;
                        } else if (isMulDiv(iter->_ope)) {
                            found = false; continue;
                        }
                    }
                } else {
                        // findex_ope is MUL
                    if (iter->_des == i) {
                            // f_des is modified, cannot push
                        found = false; continue;
                    } else {
                        if (iter->_src == i) {
                            if (isMulDiv(iter->_ope))
                                rotate = true;
                            else {
                                found = false; continue;
                            }
                        }
                    }
                }

                if (rotate) {
                    auto nindex(findex); ++nindex;
                    if (nindex != iter) { // Has to rotate
// dprint(std::clog << "######## BEF ROT\n", Program) << std::endl;
// dprint(std::clog << "### ROT " << findex->_var << i << ':', *findex);
// dprint(std::clog << " | ", *nindex);
// dprint(std::clog << " --> ", *iter) << std::endl;
                    std::rotate(findex, nindex, iter);
// dprint(std::clog << "######## AFT ROT\n", Program) << std::endl;
                    }
                    found = false;
                }
            } else {
                    // Try the next subregion of the Program
                if ( (iter->_ope != ' ') && (iter->_src == i) ) {
                    found = true;
                    findex = iter;
                }
            }
        }

            // Can be moved at the end
        if (found) {
            auto nindex(findex); ++nindex;
            if (nindex != Program.end()) {
// dprint(std::clog << "######## BEF ROT\n", Program) << std::endl;
// dprint(std::clog << "### ROT " << findex->_var << i << ':', *findex);
// dprint(std::clog << " | ", *nindex)
//                  << " --> END" << std::endl;
                std::rotate(findex, nindex, Program.end());
// dprint(std::clog << "######## AFT ROT\n", Program) << std::endl;
            }
        }
    }
}
// ===============================================================



// ===============================================================
// In-place program realizing a linear function
Tricounter LinearAlgorithm(AProgram_t& Program, const Matrix& A,
                           const char variable,
                           const bool transposed, const bool oriented) {

    const QRat& QQ = A.field();
    Tricounter opcount; // 0:ADD, 1:SCA, 2:MUL
    const size_t m(A.rowdim());
    size_t preci(A.coldim());
    for(size_t l=0; l<m; ++l) {
        const auto& CurrentRow(A[l]);
        if (CurrentRow.size()>0) {

                // choose var
            const auto aiter { nextindex(preci, CurrentRow, oriented) };
            const size_t i(aiter->first);

                // scale choosen var
            if (!QQ.isOne(aiter->second)) {
                if (transposed) {
                    if (!QQ.isMOne(aiter->second)) {
                        Program.emplace_back(variable,i,'/',aiter->second);
                    }
                } else {
                    Program.emplace_back(variable,i,'*',aiter->second);
                }
            }

                // Add the rest of the row
            for(auto iter=CurrentRow.begin(); iter != CurrentRow.end(); ++iter){
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
            for(auto iter=CurrentRow.begin(); iter != CurrentRow.end(); ++iter){
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
            if (!QQ.isOne(aiter->second)) {
                if (transposed) {
                    if (!QQ.isMOne(aiter->second)) {
                        Program.emplace_back(variable,i,'*',aiter->second);
                    }
                } else {
                    Program.emplace_back(variable,i,'/',aiter->second);
                }
            }

            if (CurrentRow.size()>1) preci = i;
        } else {
            Program.emplace_back(' ',l,' ',QQ.zero);
        }
    }

        // Removing noops
    Program.erase(std::remove_if(Program.begin(), Program.end(), isMulDivOne),
                  Program.end());



//     std::clog << std::string(20,'#') << std::endl;
//     std::clog << std::string(10,'#') << " BEF TRSP SIMP" << std::endl;
//     dprint(std::clog, Program) << std::endl;

    bool simp; do {
        // Moving variables towards themselves whenever possible
        if (transposed) pushvariables(Program, A.coldim());
        // Removes atoms followed by their inverses, one at a time
        simp = simplify(Program, transposed);
    } while(simp);

//     std::clog << std::string(20,'#') << std::endl;
//     std::clog << std::string(10,'#') << " AFT TRSP SIMP" << std::endl;
//     dprint(std::clog, Program) << std::endl;
//     std::clog << std::string(20,'#') << std::endl;

    return complexity(QQ,Program);
}
// ===============================================================



Tricounter TransposedDoubleAlgorithm(AProgram_t& Program, const Matrix& T,
                                     const char variable) {

// T.write(std::clog << "TDA:\n",FileFormat::Pretty) << std::endl;
    const QRat& QQ = T.field();
    Tricounter opcount;			// 0:ADD, 1:SCA, 2:MUL
    const size_t m(T.rowdim());		// should be even
    for(size_t l=0; l<m; ++l) {
        const auto& UpperRow(T[l]);
        const auto& LowerRow(T[++l]);
        if (UpperRow.size()>0) {

                // choose matrix <<a|c>,<0|a>>
            const auto aiter { UpperRow.begin() };
            const auto i { aiter->first };
            const auto cindex(i+1);

                // compute inverse matrix  <<y|z>,<0|y>>=<<a|c>,<0|a>>^{-1}
            QRat::Element c(QQ.zero), z(QQ.zero), y; QQ.init(y);
            const QRat::Element& a(aiter->second);
            QQ.inv(y,a);
            if ((UpperRow.size()>1) && ((aiter+1)->first == cindex)){
                c = (aiter+1)->second;
                QQ.mul(z,y,c);
                QQ.mulin(z,y);
                QQ.negin(z);   // z = - a^{-1} c a^{-1}
            }

                // scale choosen var
            if (notAbsOne(QQ,y))
                Program.emplace_back(variable, cindex, '*', y);
            if (!QQ.isZero(z))
                Program.emplace_back(variable, cindex, MONEOP('+',y), z, i);
            if (notAbsOne(QQ,y))
                Program.emplace_back(variable, i, '*', y);

                // Add the rest of the row
            for(auto iter=UpperRow.begin(); iter != UpperRow.end(); ++iter){
                if (iter != aiter) {
                    if (iter->first != cindex) {
                        Program.emplace_back(variable, iter->first,
                                             MONEOP('-',y), iter->second, i);
                    }
                    Program.emplace_back(variable, iter->first+1,
                                         MONEOP('-',y), iter->second, cindex);
                }
            }

                // placeholder for barrier:
                //   (AXPY will have to performed at this point)
            Program.emplace_back(variable,i,' ', aiter->second);
            Program.emplace_back(variable,cindex,' ', aiter->second);

                // Sub the rest of the row
            for(auto iter=UpperRow.begin(); iter != UpperRow.end(); ++iter){
                if (iter != aiter) {
                    if (iter->first != cindex) {
                        Program.emplace_back(variable, iter->first,
                                             MONEOP('+',a), iter->second, i);
                    }
                    Program.emplace_back(variable, iter->first+1,
                                         MONEOP('+',a), iter->second, cindex);
                }
            }

                // Un-scale choosen var
            if (notAbsOne(QQ,a))
                Program.emplace_back(variable, cindex, '*', a);
            if (!QQ.isZero(c))
                Program.emplace_back(variable, cindex, MONEOP('+',a), c, i);
            if (notAbsOne(QQ,a))
                Program.emplace_back(variable, i, '*', a);

        } else {
            Program.emplace_back(' ',l,' ',QQ.zero);
        }
    }
// dprint(std::clog << "Program:\n", Program) << std::endl;

        // Removing noops
    Program.erase(std::remove_if(Program.begin(), Program.end(), isMulDivOne),
                  Program.end());

    bool simp; do {
        // Moving variables towards themselves whenever possible
        pushvariables(Program, T.coldim());
        // Removes atoms followed by their inverses, one at a time
        simp = simplify(Program, true);
    } while(simp);

    return complexity(QQ,Program);
}
// ===============================================================



// ===============================================================
// Searching the space of in-place linear programs
Tricounter SearchLinearAlgorithm(AProgram_t& Program,
                                 LinBox::Permutation<QRat>& P,
                                 const Matrix& A, const char variable,
                                 size_t randomloops, const bool transposed) {
    const QRat& QQ = A.field();

    Givaro::Timer elapsed;
    std::ostringstream sout, matout;

    Tricounter nbops{ LinearAlgorithm(Program, A, variable, transposed, true) };

    std::string res(sout.str());
#ifdef VERBATIM_PARSING
    std::clog << "# Oriented number of operations for " << variable
              << ": " << nbops << std::endl;
#endif

#pragma omp parallel for shared(A,Program,P,nbops)
    for(size_t i=0; i<randomloops; ++i) {

            // =============================================
            // random permutation of the rows
        LinBox::Permutation<QRat> lP(QQ,A.rowdim());
        Givaro::GivRandom generator;
        lP.random(generator.seed());

            // =============================================
            // Apply permutation to A
        Matrix pA(QQ,A.rowdim(), A.coldim());
        permuteRows(pA,lP,A,QQ);

        AProgram_t lProgram;
        Tricounter lops { LinearAlgorithm(lProgram, pA, variable, transposed) };

#pragma omp critical
        {
            if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
                 ( (std::get<0>(lops)==std::get<0>(nbops))
                   && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
                nbops = lops;
                Program = lProgram;
                P.getStorage() = std::move(lP.getStorage());
#ifdef VERBATIM_PARSING
                std::clog << "# Found algorithm[" << i << "] for " << variable
                          << ", operations: " << lops << std::endl;
#endif
            }
        }


        lops = LinearAlgorithm(lProgram, pA, variable, transposed, true);

#pragma omp critical
        {
            if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
                 ( (std::get<0>(lops)==std::get<0>(nbops))
                   && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
                nbops = lops;
                Program = lProgram;
                P.getStorage() = std::move(lP.getStorage());
#ifdef VERBATIM_PARSING
                std::clog << "# Found oriented[" << i << "] for " << variable
                          << ", operations: " << lops << std::endl;
#endif
            }
        }
    }

    return nbops;
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
                            const Matrix& T, const bool oriented,
                            const bool expanded) {
    const QRat& QQ(T.field());

    AProgram_t aprog, bprog, cprog;


    Tricounter aops { LinearAlgorithm(aprog, A, 'a', false, oriented) };
    out << "# Found " << aops << " for a" << std::endl;

    Tricounter bops { LinearAlgorithm(bprog, B, 'b', false, oriented) };
    out << "# Found " << bops << " for b" << std::endl;

    Tricounter cops;
    if (expanded) {
        cops = TransposedDoubleAlgorithm(cprog, T, 'c');
    } else {
        cops = LinearAlgorithm(cprog, T, 'c', true, oriented);
    }
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

            if (expanded) {
                MULTD(out, ckter->_var, ckter->_src, (ckter+1)->_src,
                      MONEOP('+',ckter->_val),
                      aiter->_var, aiter->_src, bjter->_var, bjter->_src, pty);
                ++ckter;
            } else {
                MUL(out, ckter->_var, ckter->_src, MONEOP('+',ckter->_val),
                    aiter->_var, aiter->_src, bjter->_var, bjter->_src, pty);
            }
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

    if (expanded) { std::get<2>(cops) >>= 1;  } // MUL2D is counted twice
    Tricounter nops; nops += aops; nops += bops; nops += cops;
    std::get<2>(nops) /= 3; // MUL is counted in each of the three programs

    return nops;
}
// ===============================================================


// ===============================================================
// Searching the space of in-place trilinear programs
Tricounter SearchTriLinearAlgorithm(std::ostream& out,
                                    const Matrix& A, const Matrix& B,
                                    const Matrix& T,
                                    size_t randomloops, bool expanded) {

    const size_t expArowdim(expanded?A.rowdim()<<1:A.rowdim());
    if ( (A.rowdim() != B.rowdim()) || (expArowdim != T.rowdim()) ) {
        std::cerr << "Incorrect dimensions :" << std::endl;
		std::cerr << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;
		std::cerr << "B is " << B.rowdim() << " by " << B.coldim() << std::endl;
		std::cerr << "C is " << T.coldim() << " by " << T.rowdim() << std::endl;
        return Tricounter{0,0,0};
    }

    const QRat& QQ = T.field();

    Givaro::Timer elapsed;
    std::ostringstream sout, matout;
    Tricounter nbops{ TriLinearProgram(sout, A, B, T, true, expanded) };
    std::string res(sout.str());
    std::clog << "# Oriented number of operations: " << nbops << std::endl;

//     std::clog << res << std::endl;
//     std::clog << std::string(40,'#') << std::endl;

#pragma omp parallel for shared(A,B,T,res,nbops)
    for(size_t i=0; i<randomloops; ++i) {

            // =============================================
            // random permutation of the rows
        LinBox::Permutation<QRat> P(QQ,A.rowdim());
        Givaro::GivRandom generator;
        P.random(generator.seed());



            // =============================================
            // Apply permutation to A, B, T
        Matrix pA(QQ,A.rowdim(), A.coldim()), pB(QQ,B.rowdim(), B.coldim()),
            pT(QQ,T.rowdim(), T.coldim());

        permuteRows(pA,P,A,QQ); permuteRows(pB,P,B,QQ);

        if (expanded) {
            LinBox::Permutation<QRat> PP(QQ,T.rowdim());
            auto & PPindices(PP.getStorage());
            const auto& Pindices(P.getStorage());
            for(size_t i(0); i<Pindices.size(); ++i) {
                PPindices[i<<1]     = Pindices[i]<<1;
                PPindices[1+(i<<1)] = 1 + (Pindices[i]<<1);
            }
            permuteRows(pT,PP,T,QQ);
        } else {
            permuteRows(pT,P,T,QQ);
        }



            // =============================================
            // random coherent negation of the rows
        for(size_t i=0; i<pA.rowdim(); ++i) {
            const bool negA( generator.brand() ), negB( generator.brand() );
            if (negA) negRow(pA, i);
            if (negB) negRow(pB, i);
            if (negA != negB) {
                if (expanded) {
                    negRow(pT, i<<1);
                    negRow(pT, 1+(i<<1));
                } else {
                    negRow(pT, i);
                }
            }

        }

        std::ostringstream lout, sout;

            // =============================================
            // trying first a directed selection of variables
        Tricounter lops{ TriLinearProgram(lout, pA, pB, pT, true, expanded) };

#pragma omp critical
        {
            if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
                 ( (std::get<0>(lops)==std::get<0>(nbops))
                   && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
                nbops = lops;
                res = lout.str();
                std::clog << "# Found algorithm[" << i << "], operations: "
                          << lops << std::endl;


//     std::clog << res << std::endl;
//     std::clog << std::string(40,'#') << std::endl;

            }
        }

            // =============================================
            // then trying a random selection of variables
        lops = TriLinearProgram(sout, pA, pB, pT, false, expanded);

#pragma omp critical
        {
            if ( (std::get<0>(lops)<std::get<0>(nbops)) ||
                 ( (std::get<0>(lops)==std::get<0>(nbops))
                   && (std::get<1>(lops)<std::get<1>(nbops)) ) ) {
                nbops = lops;
                res = sout.str();
                std::clog << "# Found oriented [" << i << "], operations: "
                          << lops << std::endl;


//     std::clog << res << std::endl;
//     std::clog << std::string(40,'#') << std::endl;

            }
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
void InitializeVariable(const char L, const size_t m, const char M)  {
    for(size_t h=0; h<m; ++h)
        std::clog << L << h << ":=" << M << "[" << (h+1) << "]:";
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
    std::clog << '<';
    for(size_t h=0; h<m; ++h) {
        if (h>0) std::clog << '|';
        std::clog << L << h << '-' << M << '[' << (h+1) << ']';
    }
    std::clog << ">;" << std::endl;
}

void CollectVariables(const char L, const size_t m,
                      const char H, const size_t n,
                      const char F, const size_t s) {
    std::clog << std::string(30,'#') << std::endl;
    for(size_t h=0; h<s; ++h)
        std::clog << "R[" << (h+1) << "]:=simplify(" << F << h << ",symbolic):";
    std::clog << std::endl;
    std::clog << L << "Input:="; CollectVariable(L, m, 'L');
    std::clog << H << "Input:="; CollectVariable(H, n, 'H');
    std::clog << std::string(30,'#') << std::endl;
}


    // Compare program with a matrix multiplication
void CheckMatrixMultiplication(const char L, const Matrix& A,
                               const char H, const Matrix& B,
                               const char F, const Matrix& C) {
    Tricounter mkn(LRP2MM(A,B,C));
    const size_t& n(std::get<0>(mkn)), t(std::get<1>(mkn)), m(std::get<2>(mkn));
    std::clog <<"# code-checking for "
              << m << 'x' << t << 'x' << n
              << " Matrix-Multiplication" << std::endl;
    std::clog << F << "Output:=<";
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

    std::clog << "errors:=map(M->nops(select(x->x<>0,convert(M,list))),["
              << L << "Input," << H << "Input," << F << "Output]);"
              << std::endl;
}

// Produce code for the application of a matrix to a variable
std::ostream& ApplyMatrix(std::ostream& out,
                          const char O, const Matrix& A, const char I) {
//     A.write(std::clog << "AppMat:\n",FileFormat::Pretty) << std::endl;

    for(size_t i(0); i<A.rowdim(); ++i) {
        const auto row(A[i]);
        out << O << '[' << (i+1) << "]:=";
        for(auto iter(row.begin()); iter != row.end(); ++iter) {
            if (iter != row.begin()) out << '+';
            out << '(';
            if (!A.field().isOne(iter->second)) out << iter->second << '*';
            out << I << '[' << (iter->first+1) << "])";
        }
        out << ':';
    }
    return out << std::endl;
}


    // Compare program with direct linear applications
void CheckTriLinearProgram(const char L, const Matrix& AA,
                           const char H, const Matrix& BB,
                           const char F, const Matrix& CC, bool expanded) {

    std::clog <<"# code-checking for "
              << AA.coldim() << ':' << BB.coldim() << ' '
              << CC.rowdim() << 'x' << CC.coldim()
              << " trilinear operator" << std::endl;

    CollectVariables(L, AA.coldim(), H, BB.coldim(), F, CC.rowdim());

    ApplyMatrix(std::clog, 'X', AA, 'L');
    ApplyMatrix(std::clog, 'Y', BB, 'H');

    std::string zone[2] { "hig", "low" };
    for(size_t i(1); i<=AA.rowdim(); ++i) {
        std::clog << "T[" << i << "]:=X[" << i << "]*Y[" << i << ']';
        if (expanded) std::clog << '*' << zone[i & 0x1];
        std::clog << ':';
    }
    std::clog << std::endl;

    ApplyMatrix(std::clog, 'Z', CC, 'T');

    std::clog << F << "Output:=map(expand,<";
    for(size_t i(1); i<=CC.rowdim(); ++i) {
        if (i>1) std::clog << '|';
        std::clog << "F[" << i << "]+Z[" << i << "]-R[" << i << ']';
    }
    std::clog << ">);" << std::endl;

    std::clog << "errors:=map(M->nops(select(x->x<>0,convert(M,list))),["
              << L << "Input," << H << "Input," << F << "Output]);"
              << std::endl;
}

#endif
// ===============================================================
