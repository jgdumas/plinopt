// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear programs
// Authors: J-G. Dumas, B. Grenet
// ==========================================================================

/****************************************************************
 * PLinOpt Library, In-place inline implementations
 ****************************************************************/

#include "plinopt_inplace.h"

bool BiLinearAlgorithm(const Matrix& A, const Matrix& B,
                       const Matrix& T) {


    if ( (A.rowdim() != B.rowdim())
         ||
         (A.rowdim() != T.rowdim())
         ) {

        std::cerr << "Incorrect dimensions :" << std::endl;
		std::cerr << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;
		std::cerr << "B is " << B.rowdim() << " by " << B.coldim() << std::endl;
		std::cerr << "C is " << T.coldim() << " by " << T.rowdim() << std::endl;

        return false;
    }


    const QRat& QQ = T.field();
    const size_t t(A.rowdim()), m(A.coldim()), n(B.coldim()), s(T.coldim());


#ifdef VERBATIM_PARSING
    for(size_t h=0; h<m; ++h)
        std::clog << 'a' << h << ":=L[" << (h+1) << "];";
    std::clog << std::endl;
    for(size_t h=0; h<n; ++h)
        std::clog << 'b' << h << ":=H[" << (h+1) << "];";
    std::clog << std::endl;
    for(size_t h=0; h<s; ++h)
        std::clog << 'c' << h << ":=F[" << (h+1) << "];";
    std::clog << std::endl;
    std::clog << std::string(20,'#') << std::endl;
#endif



    for(size_t l=0; l<t; ++l) {
            /* Left hand side */
        auto iter(A[l].begin());
        const size_t i(iter->first);
        SCA('a',i,'*',iter->second);
        for(++iter; iter != A[l].end(); ++iter) {
            ADD('a',i,'+',iter->second, iter->first);
        }

            /* Right hand side */
        auto jter(B[l].begin());
        const size_t j(jter->first);
        SCA('b',j,'*',jter->second);
        for(++jter; jter != B[l].end(); ++jter) {
            ADD('b',j,'+',jter->second, jter->first);
        }

            /* Product */
        const auto ckter(T[l].begin());
        const size_t k(ckter->first);

        if ( (!QQ.isOne(ckter->second)) && (!QQ.isMOne(ckter->second)))
            SCA('c', k, '/', ckter->second);
        auto kter(T[l].begin());
        for(++kter; kter != T[l].end(); ++kter) {
            ADD('c', kter->first, MONEOP('-',kter->second), kter->second, k);
        }

            /* the recursive (multiplicative) calls */
        MUL('c',k,MONEOP('+',ckter->second),'a',i,'b',j);


            /* RESTORE: Product */
        kter = T[l].begin();
        for(++kter; kter != T[l].end(); ++kter) {
            ADD('c', kter->first, MONEOP('+',kter->second),  kter->second, k);
        }
        if ( (!QQ.isOne(ckter->second)) && (!QQ.isMOne(ckter->second)))
            SCA('c', ckter->first, '*', ckter->second);

            /* RESTORE: Right hand side */
        jter = B[l].begin();
        for(++jter; jter != B[l].end(); ++jter) {
            ADD('b',j,'-',jter->second, jter->first);
        }
        SCA('b',j,'/',B[l].begin()->second);

            /* RESTORE: Left hand side */
        iter = A[l].begin();
        for(++iter; iter != A[l].end(); ++iter) {
            ADD('a',i,'-',iter->second, iter->first);
        }
        SCA('a',i,'/',A[l].begin()->second);

    }

#ifdef VERBATIM_PARSING
    std::clog << std::string(20,'#') << std::endl;
    for(size_t h=0; h<s; ++h)
        std::clog << "R[" << (h+1) << "]:=simplify(" << 'c' << h << ",symbolic);";
    std::clog << std::endl;
    for(size_t h=0; h<n; ++h)
        std::clog << 'b' << h << "-H[" << (h+1) << "],";
    std::clog << "0;" << std::endl;
    for(size_t h=0; h<m; ++h)
        std::clog << 'a' << h << "-L[" << (h+1) << "],";
    std::clog << "0;" << std::endl;
    std::clog << std::string(20,'#') << std::endl;
#endif


    return true;
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
    AA.write(std::clog << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;
    BB.write(std::clog << "B:=",Tag::FileFormat::Maple) << ';' << std::endl;
    TT.write(std::clog << "T:=",Tag::FileFormat::Maple) << ';' << std::endl;
    std::clog << std::string(20,'#') << std::endl;
#endif
    
}


