// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet
// ==========================================================================

/****************************************************************
 * Returns an in-place program
 *          computing the trilinear function in HM representation
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


void usage(char ** argv) {
    std::clog << "Usage: " << argv[0] << " L.sms R.sms P.sms [-e|-m] [-O #]\n";

    std::clog << "  -m: check for a matrix multiplication\n"
              << "  -e: double expands the intermediate result\n"
              << "  -O #: randomized search with that many loops\n";
    exit(-1);
}


// ===============================================================
// argv[1-3]: L.sms R.sms P.sms
//   -e: double expand the intermediate result
//   -O # randomized search of that many loops
//        looks for reduced number of additions, then multiplications
int main(int argc, char ** argv) {
    using PLinOpt::FileFormat;

    size_t randomloops(DORANDOMSEARCH?DEFAULT_RANDOM_LOOPS:1);
    bool doexpand(false), checkmat(false);

    if (argc<4) usage(argv);

    std::string files[3]; size_t numfil(0);

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") usage(argv);
        else if (args == "-m") { checkmat = true; }
        else if (args == "-e") { doexpand = true; }
        else if (args == "-O") {
            randomloops = atoi(argv[++i]);
            if ( (randomloops>1) && (!DORANDOMSEARCH) ) {
                randomloops = 1;
                std::cerr << "#  \033[1;36mWARNING: RANDOM_TIES not defined,"
                          << " random loops disabled\033[0m." << std::endl;
            }
        } else {
            files[numfil] = std::move(args); ++numfil;
        }
    }
        // =================================
        // Reading matrices
	std::ifstream left (files[0]), right (files[1]), product(files[2]);

    PLinOpt::QRat QQ;
    PLinOpt::QMstream ls(QQ, left), rs(QQ, right), ps(QQ, product);
    PLinOpt::Matrix A(ls), B(rs), C(ps);

    PLinOpt::Matrix T(QQ); PLinOpt::Transpose(T, C);


#ifdef VERBATIM_PARSING
    A.write(std::clog << "A:=Matrix(",FileFormat::Maple) << ");" << std::endl;
    B.write(std::clog << "B:=Matrix(",FileFormat::Maple) << ");" << std::endl;
    C.write(std::clog << "C:=Matrix(",FileFormat::Maple) << ");" << std::endl;
    T.write(std::clog << "T:=Matrix(",FileFormat::Maple) << ");" << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif


    PLinOpt::Tricounter opcount; // 0:ADD, 1:SCA, 2:MUL

    if (doexpand) {
            // =================================
            // Duplicate intermediate products
            // to group them 2 by 2
        PLinOpt::Matrix AA(QQ), BB(QQ), TT(QQ);
        PLinOpt::DoubleExpand(AA,BB,TT, A,B,T);

#ifdef VERBATIM_PARSING
    A.write(std::clog << "A:\n",FileFormat::Pretty) << std::endl;
    B.write(std::clog << "B:\n",FileFormat::Pretty) << std::endl;
    TT.write(std::clog << "TT:\n",FileFormat::Pretty) << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

#ifdef INPLACE_CHECKER
        PLinOpt::InitializeVariables('a',AA.coldim(), 'b', BB.coldim(),
                                     'c', TT.coldim());
#endif
            // TODO: a SearchTriLinearAlgorithm preserving 2 by 2 products ...
        opcount = PLinOpt::SearchTriLinearAlgorithm(std::cout, A, B, TT,
                                                    randomloops, true);

#ifdef INPLACE_CHECKER
        PLinOpt::Matrix CC(QQ); PLinOpt::Transpose(CC, TT);
        PLinOpt::CheckTriLinearProgram('a', AA, 'b', BB, 'c', CC, true);
#endif

    } else {
            // =================================
            // Direct computation
#ifdef INPLACE_CHECKER
        PLinOpt::InitializeVariables('a',A.coldim(), 'b', B.coldim(),
                                     'c', T.coldim());
#endif

        opcount = PLinOpt::SearchTriLinearAlgorithm(std::cout, A, B, T,
                                                    randomloops);

#ifdef INPLACE_CHECKER
        PLinOpt::CollectVariables('a',A.coldim(), 'b', B.coldim(),
                                  'c', T.coldim());
        if (checkmat)
            PLinOpt::CheckMatrixMultiplication('a', A, 'b', B, 'c', C);
        else
            PLinOpt::CheckTriLinearProgram('a', A, 'b', B, 'c', C, false);
#endif

    }

        // =================================
        // Resulting program operation count
    std::clog << std::string(40,'#') << std::endl;
    std::clog << "# \033[1;32m" << std::get<0>(opcount) << "\tADD\033[0m\n";
    std::clog << "# \033[1;32m" << std::get<1>(opcount) << "\tSCA\033[0m\n";
    std::clog << "# \033[1;32m" << std::get<2>(opcount);
    if (doexpand)
        std::clog << "\tAXPY (double size)\033[0m\n";
    else
        std::clog << "\tAXPY\033[0m\n";
    std::clog << std::string(40,'#') << std::endl;

    return 0;
}
