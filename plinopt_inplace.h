// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet
// ==========================================================================

/****************************************************************
 * PLinOpt Library, In-place definitions
 * Reference:
 *   [ J-G. Dumas, B. Grenet; Jul. 2023
 *     In-place accumulation of fast multiplication formulae
 *     (https://hal.science/hal-04167499) ]
 ****************************************************************/

#include "plinopt_library.h"

#ifndef _PLINOPT_LIBRARY_INPLACE_H_
#define _PLINOPT_LIBRARY_INPLACE_H_

// ===============================================================
// In-place program realizing a bilinear function
Tricounter BiLinearAlgorithm(std::ostream& out,
                             const Matrix& A, const Matrix& B,
                             const Matrix& T, const bool oriented=false);

// ===============================================================
// Searching the space of in-place bilinear programs
Tricounter SearchBiLinearAlgorithm(std::ostream& out,
                                   const Matrix& A, const Matrix& B,
                                   const Matrix& T, size_t randomloops);

// ===============================================================
// Double expansion of an HM representation
void DoubleExpand(Matrix& AA, Matrix& BB, Matrix& TT,
                  const Matrix& A, const Matrix& B, const Matrix& T);



// ===============================================================
// Directed selection of accumulation i/o variables
Matrix::Row::const_iterator orientindex(const size_t preci,
                                        const Matrix::Row& L,
                                        const bool oriented);

// ===============================================================
// Random selection of accumulation i/o variables
Matrix::Row::const_iterator nextindex(const size_t preci,
                                      const Matrix::Row& L,
                                      const bool oriented);


// ===============================================================
// Definition of ADD, SCA, MUL atomic operations
// Switch (__INPLOP__) between explicit/in-place operator notation

#define VALPAR(v) '(' << v << ')'
#define SWAPOP(op) (op=='+'?'-':'+')
#define MONEOP(op,val) (QQ.isMOne(val)?SWAPOP(op):op)

// n is triple of operation count, 0:ADD, 1:SCA, 2:MUL
#ifdef __INPLOP__
#define INOPSYNT(op) ' ' << op << '=' << ' '
#define UNITSOP(out,op,val,n) out << INOPSYNT(MONEOP(op,val)); ++std::get<0>(n); if ((!QQ.isOne(val)) && (!QQ.isMOne(val))) { out << val << '*'; ++std::get<1>(n); }


#define SCA(out,c,i,op,val,n) if (! QQ.isOne(val)) { out << c << i << INOPSYNT(op) << val << ';' << std::endl; ++std::get<1>(n); }
#define ADD(out,c,i,op,val,h,n) out << c << i; UNITSOP(out,op,val,n); out << c << h << ';' << std::endl;
#define MUL(out,c,k,op,a,i,b,j,n) ++std::get<2>(n); out << c << k << INOPSYNT(op) << a << i << " . " << b << j << ';' << " ### AXPYIN ###" << std::endl;

#else
#define OPSYNT(o) ' ' << o << ' '
#define UNITSOP(out,op,val,n) out << OPSYNT(MONEOP(op,val)); ++std::get<0>(n); if ((!QQ.isOne(val)) && (!QQ.isMOne(val))) { out << VALPAR(val) << '*'; ++std::get<1>(n); }

#define SCA(out,c,i,op,val,n) if (! QQ.isOne(val)) { out << c << i << ":=" << c << i << OPSYNT(op) << VALPAR(val) << ';' << std::endl; ++std::get<1>(n); }
#define ADD(out,c,i,op,val,h,n) out << c << i << ":=" << c << i; UNITSOP(out,op,val,n); out << c << h << ';' << std::endl;
#define MUL(out,c,k,op,a,i,b,j,n) ++std::get<2>(n); out << c << k << ":=" << c << k << OPSYNT(op) << a << i << " * " << b << j << ';' << " ### AXPY ###" << std::endl;
#endif


// ===============================================================
// Tools

    // Removing newlines in strings
std::string rmnl(const std::string& str) {
    std::string s(str);
    s.erase(std::remove(s.begin(), s.end(), '\n'), s.cend());
    return s;
}


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


#include "plinopt_inplace.inl"

#endif
