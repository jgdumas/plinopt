// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, B. Grenet
// ==========================================================================

/****************************************************************
 * PLinOpt Library, In-place inline implementations
 ****************************************************************/

#include "plinopt_inplace.h"


Matrix::Row::const_iterator nextindex(const size_t preci, const Matrix::Row& L) {
        // Try to reuse previous variable
    auto nexti(std::find_if(L.begin(),L.end(),
                            [preci](const auto&a) { return a.first == preci; } )
               );

        // Otherwise, prefer a variable with coefficient One
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

std::string rmnl(const std::string& str) {
    std::string s(str);
    s.erase(std::remove(s.begin(), s.end(), '\n'), s.cend());
    return s;
}


Tricounter BiLinearAlgorithm(std::ostream& out,
                             const Matrix& A, const Matrix& B,
                             const Matrix& T) {


    if ( (A.rowdim() != B.rowdim())
         ||
         (A.rowdim() != T.rowdim())
         ) {

        std::cerr << "Incorrect dimensions :" << std::endl;
		std::cerr << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;
		std::cerr << "B is " << B.rowdim() << " by " << B.coldim() << std::endl;
		std::cerr << "C is " << T.coldim() << " by " << T.rowdim() << std::endl;

        return Tricounter{0,0,0};
    }


    const QRat& QQ = T.field();


#ifdef INPLACE_CHECKER
    const size_t tt(A.coldim()), n(B.coldim()), s(T.coldim());
    for(size_t h=0; h<tt; ++h)
        std::clog << 'a' << h << ":=L[" << (h+1) << "];";
    std::clog << std::endl;
    for(size_t h=0; h<n; ++h)
        std::clog << 'b' << h << ":=H[" << (h+1) << "];";
    std::clog << std::endl;
    for(size_t h=0; h<s; ++h)
        std::clog << 'c' << h << ":=F[" << (h+1) << "];";
    std::clog << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

    Tricounter opcount; // 0:ADD, 1:SCA, 2:MUL
    size_t preci(A.coldim()), precj(B.coldim()), preck(T.coldim());

    std::ostringstream saout, sbout, scout;

    std::map<std::pair<size_t,Givaro::Rational>,std::ostringstream> maaout, mabout, macout;

    const size_t m(A.rowdim());
    for(size_t l=0; l<m; ++l) {
// std::clog << "# BEG row: " << l << ", astr: " << rmnl(saout.str()) << ", bstr: " <<  rmnl(sbout.str()) << std::endl;

           /* Left hand side */
        const auto aiter { nextindex(preci, A[l]) };
        const size_t i(aiter->first);
        if ( (i == preci) && (l>0) && (aiter->second == A.getEntry(l-1,i)) ) {
            if (! QQ.isOne(aiter->second)) {
#ifdef VERBATIM_PARSING
                std::clog << "# Optimized out: scalar (div;mul), "
                          << 'a' << i << " /* " << VALPAR(aiter->second) << std::endl;
#endif
                --std::get<1>(opcount);
                saout.clear(); saout.str(std::string());
            }
        } else {
            SCA(saout, 'a',i,'*',aiter->second, opcount);
        }
        const bool aSCA(saout.tellp() != std::streampos(0));

        for(auto iter = A[l].begin(); iter != A[l].end(); ++iter) {
            if (iter != aiter) {
                if ( (i == preci) && (l>0) && (! aSCA) &&
                     (std::find(A[l-1].begin(), A[l-1].end(), *iter ) != A[l-1].end()) ) {
#ifdef VERBATIM_PARSING
                    std::clog << "# Optimized out: (add;sub), "
                              << 'a' << i << " +- "
                              << VALPAR(iter->second) << 'a' << iter->first << std::endl;
#endif
                    maaout.erase(maaout.find(*iter));
                    --std::get<0>(opcount);
                } else {
                    ADD(saout, 'a',i,'+',iter->second, iter->first, opcount);
                }
            }
        }
//         std::clog << "### AFT maaout" << std::endl;
        for(const auto& [key,value]: maaout) {
            out << value.str() << std::flush;
        }
        maaout.clear();
//         std::clog << "### AFT saout" << std::endl;
        out << saout.str() << std::flush;
        saout.clear(); saout.str(std::string());

            /* Right hand side */
        const auto bjter( nextindex(precj, B[l]) );
        const size_t j(bjter->first);
        if ( (j == precj) && (l>0) && (bjter->second == B.getEntry(l-1,j)) ) {
            if (! QQ.isOne(bjter->second)) {
#ifdef VERBATIM_PARSING
                std::clog << "# Optimized out, scalar (div;mul): "
                          << 'b' << j << " /* " << VALPAR(bjter->second) << std::endl;
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
                if ( (j == precj) && (l>0) && (! bSCA) &&
                     (std::find(B[l-1].begin(), B[l-1].end(), *jter ) != B[l-1].end()) ) {
#ifdef VERBATIM_PARSING
                    std::clog << "# Optimized out: (add;sub), "
                              << 'b' << j << " +- "
                              << VALPAR(jter->second) << 'b' << jter->first << std::endl;
#endif
                    mabout.erase(mabout.find(*jter));
                    --std::get<0>(opcount);
                } else {
                    ADD(sbout, 'b',j,'+',jter->second, jter->first, opcount);
                }
            }
        }
//         std::clog << "### AFT mabout" << std::endl;
        for(const auto& [key,value]: mabout) {
            out << value.str() << std::flush;
        }
        mabout.clear();
//         std::clog << "### AFT sbout" << std::endl;
        out << sbout.str() << std::flush;
        sbout.clear(); sbout.str(std::string());


            /* Product */
        const auto ckter( nextindex(preck, T[l]) );
        const size_t k(ckter->first);

        if ( (!QQ.isOne(ckter->second)) && (!QQ.isMOne(ckter->second))) {
            if ( (k == preck) && (l>0) && (ckter->second == T.getEntry(l-1,k)) ) {
                if (! QQ.isOne(ckter->second)) {
#ifdef VERBATIM_PARSING
                    std::clog << "# Optimized out, scalar (div;mul): "
                              << 'c' << k << " /* " << VALPAR(ckter->second) << std::endl;
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
                if ( (k == preck) && (l>0) && (! cSCA) &&
                     (std::find(T[l-1].begin(), T[l-1].end(), *kter ) != T[l-1].end()) ) {
#ifdef VERBATIM_PARSING
                    std::clog << "# Optimized out: (add;sub), "
                              << 'c' << kter->first << " +- "
                              << VALPAR(kter->second) << 'c' << k << std::endl;
#endif
                    macout.erase(macout.find(*kter));
                    --std::get<0>(opcount);
                } else {
                    ADD(scout, 'c', kter->first, MONEOP('-',ckter->second),
                        kter->second, k, opcount);
                }
            }
        }
//         std::clog << "### AFT macout" << std::endl;
        for(const auto& [key,value]: macout) {
            out << value.str() << std::flush;
        }
        macout.clear();
//         std::clog << "### AFT scout" << std::endl;
        out << scout.str() << std::flush;
        scout.clear(); scout.str(std::string());

            /* the recursive (multiplicative) calls */
        MUL(out, 'c',k,MONEOP('+',ckter->second),'a',i,'b',j, opcount);


            /* RESTORE: Product */
        for(auto kter = T[l].begin(); kter != T[l].end(); ++kter) {
            if (kter != ckter) {
                ADD(macout[*kter], 'c', kter->first, MONEOP('+',ckter->second),
                    kter->second, k, opcount);
            }
        }
        if ( (!QQ.isOne(ckter->second)) && (!QQ.isMOne(ckter->second))) {
            SCA(scout, 'c', ckter->first, '*', ckter->second, opcount);
        }

            /* RESTORE: Right hand side */
        for(auto jter = B[l].begin(); jter != B[l].end(); ++jter) {
            if (jter != bjter) {
                ADD(mabout[*jter], 'b',j,'-',jter->second, jter->first, opcount);
            }
        }
        SCA(sbout, 'b',j,'/',bjter->second, opcount);

            /* RESTORE: Left hand side */
        for(auto iter = A[l].begin(); iter != A[l].end(); ++iter) {
            if (iter != aiter) {
                ADD(maaout[*iter], 'a',i,'-',iter->second, iter->first, opcount);
            }
        }
        SCA(saout, 'a',i,'/',aiter->second, opcount);

// std::clog << "# END row: " << l << std::endl;

        preci=i; precj=j; preck=k;
    }

//     std::clog << "### END macout" << std::endl;
    for(const auto& [key,value]: macout) {
        out << value.str() << std::flush;
    }
//     std::clog << "### END scout" << std::endl;
    out << scout.str() << std::flush;


//     std::clog << "### END mabout" << std::endl;
    for(const auto& [key,value]: mabout) {
        out << value.str() << std::flush;
    }
//     std::clog << "### END sbout" << std::endl;
    out << sbout.str() << std::flush;


//     std::clog << "### END maaout" << std::endl;
    for(const auto& [key,value]: maaout) {
        out << value.str() << std::flush;
    }
//     std::clog << "### END saout" << std::endl;
    out << saout.str() << std::flush;


#ifdef INPLACE_CHECKER
    std::clog << std::string(30,'#') << std::endl;
    for(size_t h=0; h<s; ++h)
        std::clog << "R[" << (h+1) << "]:=simplify(" << 'c' << h << ",symbolic);";
    std::clog << std::endl;
    for(size_t h=0; h<n; ++h)
        std::clog << 'b' << h << "-H[" << (h+1) << "],";
    std::clog << "0;" << std::endl;
    for(size_t h=0; h<tt; ++h)
        std::clog << 'a' << h << "-L[" << (h+1) << "],";
    std::clog << "0;" << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif


    return opcount;
}


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
    AA.write(std::clog << "A:=",LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
    BB.write(std::clog << "B:=",LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
    TT.write(std::clog << "T:=",LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

}
