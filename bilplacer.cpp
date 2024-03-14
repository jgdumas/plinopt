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


void usage(int argc, char ** argv) {
    std::clog << "Usage: " << argv[0] << " L.sms R.sms P.sms [-e] [-O #]\n";

    std::clog << "  -e: double expands the intermediate result\n"
              << "  -O #: randomized search with that many loops\n";
    exit(-1);
}


// ===============================================================
// argv[1-3]: L.sms R.sms P.sms
//   -e: double expand the intermediate result
//   -O # randomized search of that many loops
//        looks for reduced number of additions, then multiplications
int main(int argc, char ** argv) {

    size_t randomloops(DORANDOMSEARCH?ceil(std::sqrt(DEFAULT_RANDOM_LOOPS)):1);
    bool doexpand(false);

    if (argc<4) usage(argc,argv);

    for (int i = 4; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") usage(argc,argv);
        else if (args == "-e") { doexpand = true; }
        else if (args == "-O") {
            randomloops = atoi(argv[++i]);
            if ( (randomloops>1) && (!DORANDOMSEARCH) ) {
                randomloops = 1;
                std::cerr << "#  \033[1;36mWARNING: RANDOM_TIES not defined,"
                          << " random loops disabled\033[0m." << std::endl;
            }
        }
    }
        // =================================
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

    if (doexpand) {
            // =================================
            // Duplicate intermediate products
            // to group them 2 by 2
        Matrix AA(QQ), BB(QQ), TT(QQ);
        DoubleExpand(AA,BB,TT, A,B,T);
        opcount = BiLinearAlgorithm(std::cout, AA, BB, TT);
    } else {
            // =================================
            // Direct computation
#ifdef INPLACE_CHECKER
        InitializeVariables('a',A.coldim(), 'b', B.coldim(), 'c', T.coldim());
#endif

        opcount = SearchBiLinearAlgorithm(std::cout, A, B, T, randomloops);

#ifdef INPLACE_CHECKER
        CollectVariables('a',A.coldim(), 'b', B.coldim(), 'c', T.coldim());
        CheckMatrixMultiplication(A,B,C);
#endif

    }

        // =================================
        // Resulting program operation count
    std::clog << std::string(40,'#') << std::endl;
    std::clog << "# \033[1;32m" << std::get<0>(opcount) << "\tADD\033[0m\n";
    std::clog << "# \033[1;32m" << std::get<1>(opcount) << "\tSCA\033[0m\n";
    std::clog << "# \033[1;32m" << std::get<2>(opcount) << "\tAXPY\033[0m\n";
    std::clog << std::string(40,'#') << std::endl;

    return 0;
}
