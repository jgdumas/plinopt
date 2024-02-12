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
bool BiLinearAlgorithm(const Matrix& A, const Matrix& B,
                       const Matrix& T);


// ===============================================================
// Double expansion of an HM representation
void DoubleExpand(Matrix& AA, Matrix& BB, Matrix& TT,
                  const Matrix& A, const Matrix& B, const Matrix& T);


// ===============================================================
// Definition of ADD, SCA, MUL atomic operations
// Switch (__INPLOP__) between explicit/in-place operator notation

#define SWAPOP(op) (op=='+'?'-':'+')
#define MONEOP(op,val) (QQ.isMOne(val)?SWAPOP(op):op)

#ifdef __INPLOP__
#define INOPSYNT(op) ' ' << op << '=' << ' '
#define UNITSOP(op,val) std::cout << INOPSYNT(MONEOP(op,val)); if ((!QQ.isOne(val)) && (!QQ.isMOne(val))) std::cout << val << '*';


#define SCA(c,i,op,val) if (! QQ.isOne(val)) std::cout << c << i << INOPSYNT(op) << val << ';' << std::endl;
#define ADD(c,i,op,val,h) std::cout << c << i; UNITSOP(op,val); std::cout << c << h << ';' << std::endl;
#define MUL(c,k,op,a,i,b,j) std::cout << c << k << INOPSYNT(op) << a << i << " . " << b << j << ';' << std::endl;

#else
#define OPSYNT(o) ' ' << o << ' '
#define VALPAR(v) '(' << v << ')'
#define UNITSOP(op,val) std::cout << OPSYNT(MONEOP(op,val)); if ((!QQ.isOne(val)) && (!QQ.isMOne(val))) std::cout << VALPAR(val) << '*';


#define SCA(c,i,op,val) if (! QQ.isOne(val)) std::cout << c << i << ":=" << c << i << OPSYNT(op) << VALPAR(val) << ';' << std::endl;
#define ADD(c,i,op,val,h) std::cout << c << i << ":=" << c << i; UNITSOP(op,val); std::cout << c << h << ';' << std::endl;
#define MUL(c,k,op,a,i,b,j) std::cout << c << k << ":=" << c << k << OPSYNT(op) << a << i << " * " << b << j << ';' << std::endl;
#endif


// ===============================================================

#include "plinopt_inplace.inl"

#endif
