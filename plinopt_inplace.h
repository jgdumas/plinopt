// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
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
// Definition of ADD, SCA, MUL atomic operations
// Switch (__INPLOP__) between explicit/in-place operator notation

#define VALPAR(v) '(' << v << ')'
#define SWAPOP(op) (op=='+'?'-':'+')
#define MONEOP(op,val) (QQ.isMOne(val)?SWAPOP(op):op)

// n is operation count, 0:ADD, 1:SCA, 2:MUL
#ifdef __INPLOP__
#define INOPSYNT(op) ' ' << op << '=' << ' '
#define UNITSOP(out,op,val,n) out << INOPSYNT(MONEOP(op,val)); ++std::get<0>(n); if ((!QQ.isOne(val)) && (!QQ.isMOne(val))) { out << val << '*'; ++std::get<1>(n); }


#define SCA(out,c,i,op,val,n) if (! QQ.isOne(val)) { out << c << i << INOPSYNT(op) << val << ';' << std::endl; ++std::get<1>(n); }
#define ADD(out,c,i,op,val,h,n) out << c << i; UNITSOP(out,op,val,n); out << c << h << ';' << std::endl;
#define MUL(out,c,k,op,a,i,b,j,n) out << c << k << INOPSYNT(op) << a << i << " . " << b << j << ';' << std::endl; ++std::get<2>(n);

#else
#define OPSYNT(o) ' ' << o << ' '
#define UNITSOP(out,op,val,n) out << OPSYNT(MONEOP(op,val)); ++std::get<0>(n); if ((!QQ.isOne(val)) && (!QQ.isMOne(val))) { out << VALPAR(val) << '*'; ++std::get<1>(n); }

#define SCA(out,c,i,op,val,n) if (! QQ.isOne(val)) { out << c << i << ":=" << c << i << OPSYNT(op) << VALPAR(val) << ';' << std::endl; ++std::get<1>(n); }
#define ADD(out,c,i,op,val,h,n) out << c << i << ":=" << c << i; UNITSOP(out,op,val,n); out << c << h << ';' << std::endl;
#define MUL(out,c,k,op,a,i,b,j,n) out << c << k << ":=" << c << k << OPSYNT(op) << a << i << " * " << b << j << ';' << std::endl; ++std::get<2>(n);
#endif


// ===============================================================


#include "plinopt_inplace.inl"

#endif
