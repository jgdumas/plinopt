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
// Inplace Program is a vector of atomic operations (ADD or SCA)
struct Atom;
typedef std::vector<Atom> AProgram_t;

// Printing a runable program
std::ostream& operator<< (std::ostream&, const AProgram_t&);
// Print a runable program with explicitly permuted output
std::ostream& Pprint(std::ostream& out, const char c, const AProgram_t& atomP,
                     const LinBox::Permutation<QRat>& P);
// Raw print program (including barrier)
std::ostream& dprint (std::ostream& out, const AProgram_t& p);

Tricounter complexity(const AProgram_t& p);

// ===============================================================
// In-place program realizing a linear function
Tricounter LinearAlgorithm(AProgram_t& Program, const Matrix& A,
                           const char a, const bool transposed=false,
                           const bool oriented=false);

// ===============================================================
// Searching the space of in-place linear programs
Tricounter SearchLinearAlgorithm(AProgram_t& Program, LinBox::Permutation<QRat>& P, const Matrix& A,
                                 const char a, size_t randomloops, const bool transposed=false);

// ===============================================================
// In-place optimized program realizing a trilinear function
Tricounter TriLinearProgram(std::ostream& out, const Matrix& A, const Matrix& B,
                            const Matrix& T, const bool oriented=false,
                            const bool expanded=false);

// ===============================================================
// Searching the space of in-place trilinear programs
Tricounter SearchTriLinearAlgorithm(std::ostream& out,
                                    const Matrix& A, const Matrix& B,
                                    const Matrix& T, size_t randomloops,
                                    const bool expanded=false);

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

auto isAddSub { [](const char c) { return ( (c=='+') || (c=='-') ); } };
auto isMulDiv { [](const char c) { return ( (c=='*') || (c=='/') ); } };

// ===============================================================
// Definition of ADD, SCA, MUL atomic operations
// Switch (__INPLOP__) between explicit/in-place operator notation

#define VALPAR(v) '(' << v << ')'
#define SWAPOP(op) (op=='+'?'-':'+')
#define INVOP(op) (op=='*'?'/':'*')
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
#define MULTD(out,c,k,l,op,a,i,b,j,n) ++std::get<2>(n); out << c << k << ":=" << c << k << OPSYNT(op) << '(' << a << i << " * " << b << j << ")*low;" << " ### AXPY low  ###" << std::endl; out << c << l << ":=" << c << l << OPSYNT(op) << '(' << a << i << " * " << b << j << ")*hig;" << " ### AXPY high ###" << std::endl;
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
    // Compare program with a matrix multiplication
void CheckMatrixMultiplication(const Matrix& A, const Matrix& B,
                               const Matrix& C);

    // Compare program with direct linear applications
void CheckTriLinearProgram(const char L, const Matrix& AA,
                           const char H, const Matrix& BB,
                           const char F, const Matrix& CC,
                           bool expanded = false);
#endif


#include "plinopt_inplace.inl"

#endif
