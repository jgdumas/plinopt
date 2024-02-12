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
    A.write(std::clog << "A:=",FileFormat::Maple) << ';' << std::endl;
    B.write(std::clog << "B:=",FileFormat::Maple) << ';' << std::endl;
    C.write(std::clog << "C:=",FileFormat::Maple) << ';' << std::endl;
    T.write(std::clog << "T:=",FileFormat::Maple) << ';' << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

    if (argc>4) {
            // ================================
            // Duplicate intermediate products
            // to group them 2 by 2
        Matrix AA(QQ), BB(QQ), TT(QQ);
        DoubleExpand(AA,BB,TT, A,B,T);
        BiLinearAlgorithm(AA, BB, TT);
    } else {
            // ================================
            // Direct computation
        BiLinearAlgorithm(A, B, T);
    }

#ifdef INPLACE_CHECKER
        // Enables checking If algorithm is 2x2 matrix product
    std::clog << "<H[1]*L[1] + H[3]*L[2] + F[1] - R[1]," << std::endl;
    std::clog << "H[2]*L[1] + H[4]*L[2] + F[2] - R[2]," << std::endl;
    std::clog << "H[1]*L[3] + H[3]*L[4] + F[3] - R[3]," << std::endl;
    std::clog << "H[2]*L[3] + H[4]*L[4] + F[4] - R[4]>;" << std::endl;
#endif


    return 0;
}
