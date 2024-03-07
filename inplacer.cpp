// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, B. Grenet
// ==========================================================================

/****************************************************************
 * Returns an in-place program
 *          computing the bilinear function in HM representation
 *
 * Usage: L.sms R.sms P.sms [expansion]
 *          L.sms/R.sms/P.sms the 3 HM matrices
 *          expansion: if present intermediate results grouped by 2
 *
 * For an HM representation of a matrix product,
 *     the vectorization is supposed row-major:
 *        [ a11 a12 ]
 *        [ a21 a22 ] is vectorized as [a11 a12 a21 a22]
 *
 * Reference:
 *      [ In-place accumulation of fast multiplication formulae
 *        J-G. Dumas, B. Grenet
 *        https://hal.science/hal-04167499 ]
 *
 * Matrix syntax: SMS format, see:
 *      [Sparse Integer Matrix Collection](https://hpac.imag.fr)
 *		- Starts with: `m n 'R'`
 *		- then: `i j value`
 *		- ends with: `0 0 0`
 ****************************************************************/

// ============================================================
// Define to switch between explicit/in-place operator notation
//#define __INPLOP__
// ============================================================

// ============================================================
// Define to print comments during parsing
//#define VERBATIM_PARSING
// ============================================================

// ============================================================
// Define to print verification codes
//#define INPLACE_CHECKER
// ============================================================

#include "plinopt_inplace.h"

// ===============================================================
// argv[1-3]: L.sms R.sms P.sms
// argv[4]: if present, double expand intermediate result
int main(int argc, char ** argv) {


    if ((argc <=3) || (std::string(argv[1]) == "-h")) {
        std::clog << "Usage: " << argv[0] << " L.sms R.sms P.sms [expand]\n";
        exit(-1);
    }

        // ================================
        // Reading matrices
	std::ifstream left (argv[1]), right (argv[2]), product(argv[3]);

    QRat QQ;
    QMstream ls(QQ, left), rs(QQ, right), ps(QQ, product);
    Matrix A(ls), B(rs), C(ps);

    Matrix T(QQ); Transpose(T, C);


#ifdef VERBATIM_PARSING
    A.write(std::clog << "A:=Matrix(",FileFormat::Maple) << ");" << std::endl;
    B.write(std::clog << "B:=Matrix(",FileFormat::Maple) << ");" << std::endl;
    C.write(std::clog << "C:=Matrix(",FileFormat::Maple) << ");" << std::endl;
    T.write(std::clog << "T:=Matrix(",FileFormat::Maple) << ");" << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif
    Tricounter opcount; // 0:ADD, 1:SCA, 2:MUL

    if (argc>4) {
            // ================================
            // Duplicate intermediate products
            // to group them 2 by 2
        Matrix AA(QQ), BB(QQ), TT(QQ);
        DoubleExpand(AA,BB,TT, A,B,T);
        opcount = BiLinearAlgorithm(std::cout, AA, BB, TT);
    } else {
            // ================================
            // Direct computation
        opcount = BiLinearAlgorithm(std::cout, A, B, T);
    }

    std::clog << std::string(40,'#') << std::endl;
    std::clog << "# \033[1;32m" << std::get<0>(opcount) << "\tADD" << "\033[0m" << std::endl;
    std::clog << "# \033[1;32m" << std::get<1>(opcount)  << "\tSCA" << "\033[0m" << std::endl;
    std::clog << "# \033[1;32m" << std::get<2>(opcount)  << "\tAXPY" << "\033[0m" << std::endl;
    std::clog << std::string(40,'#') << std::endl;

#ifdef INPLACE_CHECKER
        // =============================================
        // Checking, when algorithm is a matrix product
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
#endif


    return 0;
}
