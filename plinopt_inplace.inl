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

        // Otherwise, prefer One
    if (nexti == L.end()) {
        nexti = std::find_if(L.begin(),L.end(),
                             [](const auto&a) { return isOne(a.second); } );
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
//     std::clog << std::endl << "## RMNL1, str: " << str << ", s: " << s  << std::endl << std::endl;
    s.erase(std::remove(s.begin(), s.end(), '\n'), s.end());
//     std::clog << std::endl << "## RMNL2, str: " << str << ", s: " << s  << std::endl << std::endl;
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
    size_t preci(0), precj(0), preck(0);

    std::ostringstream saout, sbout;

//     std::map<Matrix::Row::const_iterator,std::ostringstream> maaout;

    const size_t m(A.rowdim());
    for(size_t l=0; l<m; ++l) {
// std::clog << "# BEG row: " << l << ", astr: " << rmnl(saout.str()) << ", bstr: " <<  rmnl(sbout.str()) << std::endl;

            /* Left hand side */
        const auto aiter { nextindex(preci, A[l]) };
        const size_t i(aiter->first);

        if ( (i == preci) && (l>0) && (aiter->second == A.getEntry(l-1,i)) ) {
            if (! QQ.isOne(aiter->second)) {
#ifdef VERBATIM_PARSING
                std::clog << "# scalar div then mul optimized out: "
                          << 'a' << i << "/*" << aiter->second << std::endl;
#endif
                --std::get<1>(opcount);
            }

        } else {
            out << saout.str(); // print previous DIV
            SCA(out, 'a',i,'*',aiter->second, opcount);
        }
        saout.clear(); saout.str(std::string());

        for(auto iter = A[l].begin(); iter != A[l].end(); ++iter) {
            if (iter != aiter) {

                if ( (i == preci) && (l>0) &&
                     (std::find(A[l-1].begin(), A[l-1].end(), *iter ) != A[l-1].end()) ) {
#ifdef VERBATIM_PARSING
            std::clog << "# add then sub could be optimized out: "
                      << 'a' << i << "+-"
                      << VALPAR(iter->second) << 'a' << iter->first << std::endl;
#endif
                }

                ADD(out, 'a',i,'+',iter->second, iter->first, opcount);
            }
        }

            /* Right hand side */
        const auto bjter( nextindex(precj, B[l]) );
        const size_t j(bjter->first);

        if ( (j == precj) && (l>0) && (bjter->second == B.getEntry(l-1,j)) ) {
            if (! QQ.isOne(bjter->second)) {
#ifdef VERBATIM_PARSING
                std::clog << "# scalar div then mul optimized out: "
                          << 'b' << j << "/*" << bjter->second << std::endl;
#endif
                --std::get<1>(opcount);
            }
        } else {
            out << sbout.str(); // print previous DIV
            SCA(out, 'b',j,'*',bjter->second, opcount);
        }
        sbout.clear(); sbout.str(std::string());

        for(auto jter = B[l].begin(); jter != B[l].end(); ++jter) {
            if (jter != bjter) {
                if ( (j == precj) && (l>0) &&
                     (std::find(B[l-1].begin(), B[l-1].end(), *jter ) != B[l-1].end()) ) {
#ifdef VERBATIM_PARSING
            std::clog << "# add then sub could be optimized out: "
                      << 'b' << j << "+-"
                      << VALPAR(jter->second) << 'b' << jter->first << std::endl;
#endif
                }
                ADD(out, 'b',j,'+',jter->second, jter->first, opcount);
            }
        }

            /* Product */
        const auto ckter( nextindex(preck, T[l]) );
        const size_t k(ckter->first);

        if ( (!QQ.isOne(ckter->second)) && (!QQ.isMOne(ckter->second))) {
            SCA(out, 'c', k, '/', ckter->second, opcount);
        }
        for(auto kter = T[l].begin(); kter != T[l].end(); ++kter) {
            if (kter != ckter) {
                ADD(out, 'c', kter->first, MONEOP('-',ckter->second),
                    kter->second, k, opcount);
            }
        }

            /* the recursive (multiplicative) calls */
        MUL(out, 'c',k,MONEOP('+',ckter->second),'a',i,'b',j, opcount);


            /* RESTORE: Product */
        for(auto kter = T[l].begin(); kter != T[l].end(); ++kter) {
            if (kter != ckter) {
                ADD(out, 'c', kter->first, MONEOP('+',ckter->second),
                    kter->second, k, opcount);
            }
        }
        if ( (!QQ.isOne(ckter->second)) && (!QQ.isMOne(ckter->second))) {
            SCA(out, 'c', ckter->first, '*', ckter->second, opcount);
        }

            /* RESTORE: Right hand side */
        for(auto jter = B[l].begin(); jter != B[l].end(); ++jter) {
            if (jter != bjter) {
                ADD(out, 'b',j,'-',jter->second, jter->first, opcount);
            }
        }
        SCA(sbout, 'b',j,'/',bjter->second, opcount);

            /* RESTORE: Left hand side */
        for(auto iter = A[l].begin(); iter != A[l].end(); ++iter) {
            if (iter != aiter) {
                ADD(out, 'a',i,'-',iter->second, iter->first, opcount);
            }
        }
        SCA(saout, 'a',i,'/',aiter->second, opcount);

// std::clog << "# END row: " << l << std::endl;

        preci=i; precj=j; preck=k;
    }


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
