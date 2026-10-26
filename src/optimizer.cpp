// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/*************************************************************************
 * Optimization of linear programs via Common Subexpression Elimination &
 *                                     Kernel Output Factorization
 * Matrix syntax: SMS format, see:
 *                [Sparse Integer Matrix Collection](https://hpac.imag.fr)
 *				- Starts with: `m n 'R'`
 *				- then: `i j value`
 *				- ends with: `0 0 0`
 * Reference:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
 *************************************************************************/

// ============================================================
// Define to print comments during parsing
//#define VERBATIM_PARSING
// ============================================================

// ============================================================
//#define RANDOM_TIES
// Acts on:
//   - Common Subexpressions selection:
//       will determine the maximal size of CSE
//       then among CSE with the maximal size, select:
//         * by default (if not defined), via a score
//         * otherwise define the following for random choice
//   - Choice of independent rows for the Kernel method;
//         * by default (not defined), start with sparsest
//         * otherwise random order
// ============================================================

#include <sstream>

#include "plinopt_optimize.h"

namespace PLinOpt {


// ============================================================
// Optimizing a linear program over a field
template<typename Field>
int FOptimiser(std::istream& input, std::ostream& output,
               const size_t randomloops,
               const bool printMaple, const bool printPretty,
               const bool tryDirect, const bool tryKernel, const bool tryLU,
               const bool tryAB, const bool mostCSE, const bool allkernels,
               const bool KFI, const Field& F, const int verbose) {

        // ============================================================
        // Read Matrix of Linear Transformation
    QRat QQ;
    QMstream ms(QQ,input);
    Matrix M(ms); M.resize(M.rowdim(),M.coldim());

    Givaro::Timer global; global.start();
    Matrix T(QQ,M.coldim(),M.rowdim()); Transpose(T,M);
    global.stop();

    if (printPretty) {
        M.write(std::clog,FileFormat::Pretty) << ';' << std::endl;
        std::clog << std::string(40,'#') << std::endl;
    }

    if (printMaple) {
        M.write(std::clog << "M:=",FileFormat::Maple) << ';' << std::endl;
        std::clog << std::string(40,'#') << std::endl;
    }

    std::clog << std::string(40,'#') << std::endl;
    if (verbose>0 && DORANDOMSEARCH) {
        std::clog << "# random seed " << randomSeedBase() << std::endl;
    }

    const Pair<size_t> opsinit(naiveOps(M));
    std::ostringstream ssout;

        // ============================================================
        // Try the different optimization methods
    Pair<size_t> nbops( OptMethods(opsinit, ssout, F, M, T, global, randomloops,
                                   printMaple, printPretty,
                                   tryDirect, tryKernel, tryLU, tryAB,
                                   mostCSE, allkernels, KFI, verbose) );

    std::stringstream ssin; ssin << ssout.str();
    reseed('C', 0);     // the compacter's draws do not depend on the loops
    PLinOpt::Compacter(output,ssin);

    if ((nbops.first !=0 || nbops.second != 0)) {
        std::clog << std::string(40,'#') << std::endl;
        std::clog << "# \033[1;32m" << nbops.first
                  << "\tadditions\tinstead of " << opsinit.first
                  << "\033[0m \t" << global << std::endl;
        std::clog << "# \033[1;32m" << nbops.second
                  << "\tmultiplications\tinstead of " << opsinit.second
                  << "\033[0m" << std::endl;
        std::clog << std::string(40,'#') << std::endl;
    }

        // ============================================================
        // Maple syntax checker
#ifdef INPLACE_CHECKER
    std::clog << std::string(40,'#') << std::endl;
    std::clog << '<';
    for(size_t i=0; i< M.rowdim(); ++i) {
        if (i != 0) std::clog << '|' << std::endl;
        std::clog << '<';
        for(size_t j=0; j<M.coldim(); ++j) {
            if (j != 0) std::clog << '+';
            std::clog << 'i' << j << "*(" << M.getEntry(i,j) << ')';
        }
        std::clog << " - o" << i << '>';
    }
    std::clog << '>';
    if (F.characteristic() > 0) std::clog << " mod " << F.characteristic();
    std::clog << ';' << std::endl;
#endif

    return 0;
}


// ============================================================
// Does the program compute the matrix, over the rationals?
bool checkProgram(const std::string& matrix, const std::string& program) {
    QRat QQ;
    std::stringstream mss(matrix);
    QMstream ms(QQ,mss);
    Matrix M(ms); M.resize(M.rowdim(),M.coldim());

    std::stringstream pss(program);
    VProgram_t PV; programParser(PV, pss);
    char next(0); parenthesisExpand(PV, next);
    Matrix A(QQ); matrixBuilder(A, PV);

    const size_t m(std::max(M.rowdim(),A.rowdim()));
    const size_t n(std::max(M.coldim(),A.coldim()));
    M.resize(m,n); A.resize(m,n);
    DenseMatrix dA(QQ,m,n), dM(QQ,m,n), R(QQ,m,n);
    any2dense(dA,A); any2dense(dM,M);
    LinBox::MatrixDomain<QRat> BMD(QQ);
    BMD.sub(R,dM,dA);
    return BMD.isZero(R);
}

// ============================================================
// Choice between modular computation, lifted computation, or rationals
int Selector(std::istream& input, const size_t randomloops,
             const bool printMaple, const bool printPretty,
             const bool tryDirect, const bool tryKernel, const bool tryLU,
             const bool tryAB, const bool mostCSE, const bool allkernels, const bool KFI,
             const Givaro::Integer& q, const bool lift, const int verbose) {
    if (! Givaro::isZero(q)) {
        Givaro::Modular<Givaro::Integer> FF(q);
        return FOptimiser(input, std::cout, randomloops, printMaple, printPretty, tryDirect,
                          tryKernel, tryLU, tryAB, mostCSE, allkernels, KFI, FF, verbose);
    }
    if (lift) {
            // Search over the lifting field, check over the rationals
        std::stringstream matrix; matrix << input.rdbuf();
        std::stringstream mss(matrix.str());
        std::ostringstream program;
        LiftField FF(LiftPrime);
        const int status( FOptimiser(mss, program, randomloops, printMaple, printPretty,
                                     tryDirect, tryKernel, tryLU, tryAB, mostCSE,
                                     allkernels, KFI, FF, verbose) );
        if (checkProgram(matrix.str(), program.str())) {
            std::cout << program.str();
            return status;
        }
        std::clog << "# \033[1;36mLifted program fails the rational check,"
                  << " redoing the search over the rationals.\033[0m" << std::endl;
        std::stringstream again(matrix.str());
        QRat QQ;
        return FOptimiser(again, std::cout, randomloops, printMaple, printPretty, tryDirect,
                          tryKernel, tryLU, tryAB, mostCSE, allkernels, KFI, QQ, verbose);
    }
    QRat QQ;
    return FOptimiser(input, std::cout, randomloops, printMaple, printPretty, tryDirect,
                      tryKernel, tryLU, tryAB, mostCSE, allkernels, KFI, QQ, verbose);
}

} // End of namespace PLinOpt
// ============================================



// ============================================================
// Main: select between file / std::cin
// -O # search for reduced number of additions and/then multiplications
//      i.e. min of random # tries (requires definition of RANDOM_TIES)
int main(int argc, char** argv) {

    bool printMaple(false), printPretty(false), mostCSE(false),
        tryKernel(false), tryLU(false), tryAB(false),
        tryDirect(false), allkernels(false), lift(false);
#if defined(KERNEL_FULL_IDENTITY)
    bool KFI(true);
#else
    bool KFI(false);
#endif

    int verbose(1);
    size_t randomloops(DORANDOMSEARCH?DEFAULT_RANDOM_LOOPS:1);
    Givaro::Integer prime(0u);

    std::string filename;

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << " [-h|-M|-P|-K|-D|-G|-E|-N|-A|-Q|-q #|-O #|-s #] [stdin|matrixfile.sms]\n";

            std::clog
                << "  -D/-K/-G: direct/kernel/LU methods (default is all)\n"
                << "  -F: Kernel Full Identity (default is " << (KFI?"yes":"no") << ")\n"
                << "  -A: factoring method (default is not)\n"
                << "  -E: exhaustive greedy CSE search (default is not)\n"
                << "  -N: exhaustive nullspace permutations (default is not)\n"
                << "  -P/-M: choose the printing format\n"
                << "  -q #: search modulo (default is Rationals)\n"
                << "  -Q: search over a word-size prime field, lift the program to the"
                << " rationals and check it (falls back to the rational search)\n"
                << "  -O #: randomized search with that many loops (default "
                << randomloops << " loops)\n"
                << "  -s #: seed of the randomized search (default is the clock;"
                << " the seed used is printed)\n";

            exit(-1);
        }
        else if (args == "-M") { printMaple = true; }
        else if (args == "-P") { printPretty = true; }
        else if (args == "-G") { tryLU = true; }
        else if (args == "-A") { tryAB = true; }
        else if (args == "-Q") { lift = true; }
        else if (args == "-D") { tryDirect = true; }
        else if (args == "-K") { tryKernel = true; }
        else if (args == "-F") { KFI = true; }
        else if (args == "-E") { mostCSE = true; }
        else if (args == "-N") { allkernels = true; }
        else if (args == "-q") { prime = Givaro::Integer(argv[++i]); }
        else if (args == "-v") { verbose = atoi(argv[++i]); }
        else if (args == "-s") { PLinOpt::RandomSeed = atoll(argv[++i]); }
        else if (args == "-O") {
            randomloops = atoi(argv[++i]);
            if ( (randomloops>1) && (!DORANDOMSEARCH) ) {
                randomloops = 1;
                std::cerr << "#  \033[1;36mWARNING: RANDOM_TIES not defined,"
                          << " random loops disabled\033[0m." << std::endl;
            }
        }
        else { filename = args; }
    }

    if (!tryKernel && !tryDirect && !tryLU) {
            // Try all
        tryLU = tryDirect = tryKernel = true;
    }

    if (filename == "") {
        return PLinOpt::Selector(std::cin, randomloops, printMaple,
                                 printPretty, tryDirect, tryKernel, tryLU,
                                 tryAB, mostCSE, allkernels, KFI, prime, lift, verbose);
    } else {
        std::ifstream inputmatrix(filename);
        if ( inputmatrix ) {
            int s=PLinOpt::Selector(inputmatrix, randomloops, printMaple,
                                    printPretty, tryDirect, tryKernel, tryLU,
                                    tryAB, mostCSE, allkernels, KFI, prime, lift, verbose);
            inputmatrix.close();
            return s;
        }
    }
    return -1;
}
// ============================================================
