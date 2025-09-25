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
int FOptimiser(std::istream& input, const size_t randomloops,
               const bool printMaple, const bool printPretty,
               const bool tryDirect, const bool tryKernel, const bool tryLU,
               const bool mostCSE, const bool allkernels, const Field& F) {

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

    const Pair<size_t> opsinit(naiveOps(M));
    std::ostringstream ssout;

        // ============================================================
        // Try the different optimization methods
    Pair<size_t> nbops( OptMethods(opsinit, ssout, F, M, T, global, randomloops,
                                   printMaple, printPretty,
                                   tryDirect, tryKernel, tryLU,
                                   mostCSE, allkernels) );

    std::cout << ssout.str() << std::flush;

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
// Choice between modular computation or over the rationals
int Selector(std::istream& input, const size_t randomloops,
             const bool printMaple, const bool printPretty,
             const bool tryDirect, const bool tryKernel, const bool tryLU,
             const bool mostCSE, const bool allkernels,
             const Givaro::Integer& q) {
    if (! Givaro::isZero(q)) {
        Givaro::Modular<Givaro::Integer> FF(q);
        return FOptimiser(input, randomloops, printMaple, printPretty,
                          tryDirect, tryKernel, tryLU, mostCSE, allkernels, FF);
    } else {
        QRat QQ;
        return FOptimiser(input, randomloops, printMaple, printPretty,
                          tryDirect, tryKernel, tryLU, mostCSE, allkernels, QQ);
    }
}

} // End of namespace PLinOpt
// ============================================



// ============================================================
// Main: select between file / std::cin
// -D/-K options select direct/kernel methods only (default is both)
// -P/-M option choose the printing format
// -E option for exhaustive greedy CSE search
// -N option for exhaustive nullspace permutation search
// -q # search for linear map modulo # (default is over the rationals)
// -O # search for reduced number of additions and/then multiplications
//      i.e. min of random # tries (requires definition of RANDOM_TIES)
int main(int argc, char** argv) {

        // ============================================================
        // Linear Transformation
    bool printMaple(false), printPretty(false), mostCSE(false),
        tryKernel(false), tryLU(false), tryDirect(false), allkernels(false);

    size_t randomloops(DORANDOMSEARCH?DEFAULT_RANDOM_LOOPS:1);
    Givaro::Integer prime(0u);

    std::string filename;

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << " [-h|-M|-P|-K|-D|-G|-E|-N|-q #|-O #] [stdin|matrixfile.sms]\n";

            std::clog
                << "  -D/-K/-G: direct/kernel/LU methods (default is all)\n"
                << "  -E: exhaustive greedy CSE search (default is not)\n"
                << "  -N: exhaustive nullspace permutations (default is not)\n"
                << "  -P/-M: choose the printing format\n"
                << "  -q #: search modulo (default is Rationals)\n"
                << "  -O #: search for reduced number of additions and/then multiplications\n";

            exit(-1);
        }
        else if (args == "-M") { printMaple = true; }
        else if (args == "-P") { printPretty = true; }
        else if (args == "-G") { tryLU = true; }
        else if (args == "-D") { tryDirect = true; }
        else if (args == "-K") { tryKernel = true; }
        else if (args == "-E") { mostCSE = true; }
        else if (args == "-N") { allkernels = true; }
        else if (args == "-q") { prime = Givaro::Integer(argv[++i]); }
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
                                 mostCSE, allkernels, prime);
    } else {
        std::ifstream inputmatrix(filename);
        if ( inputmatrix ) {
            int s=PLinOpt::Selector(inputmatrix, randomloops, printMaple,
                                    printPretty, tryDirect, tryKernel, tryLU,
                                    mostCSE, allkernels, prime);
            inputmatrix.close();
            return s;
        }
    }
    return -1;
}
// ============================================================
