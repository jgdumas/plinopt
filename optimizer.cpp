// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
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
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *     Strassen's algorithm is not optimally accurate
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
//         * by default (not defined), via a score
//         * otherwise define the following for random choice
//   - Choice of independent rows for the Kernel method;
//         * by default (not defined), start with sparsest
//         * otherwise random order
// ============================================================

#include <sstream>
#include "plinopt_optimize.h"


// ============================================================
// Optimizing a linear program
int Selector(std::istream& input, const size_t randomloops,
             const bool printMaple, const bool printPretty,
             const bool tryDirect, const bool tryKernel) {

        // ============================================================
        // Read Matrix of Linear Transformation
    QRat QQ;
    QMstream ms(QQ,input);
    Matrix M(ms); M.resize(M.rowdim(),M.coldim());

    Givaro::Timer chrono; chrono.start();

    Matrix T(QQ); Transpose(T,M);

    if (printPretty) {
        M.write(std::clog,FileFormat::Pretty) << ';' << std::endl;
        std::clog << std::string(40,'#') << std::endl;
    }

    if (printMaple) {
        M.write(std::clog << "M:=",FileFormat::Maple) << ';' << std::endl;
        std::clog << std::string(40,'#') << std::endl;
    }



        // ============================================================
        // Compute naive number of operations
    int addinit(0), mulinit(0);
    for(auto iter=M.rowBegin(); iter != M.rowEnd(); ++iter)
        addinit += std::max((int)iter->size()-1,0);
    for(auto it = M.IndexedBegin(); it != M.IndexedEnd(); ++it)
        if (!isOne(abs(it.value()))) ++mulinit;

    std::clog << std::string(40,'#') << std::endl;

    std::ostringstream ssout;
    std::pair<size_t,size_t> nbops{addinit,mulinit};

        // ============================================================
        // First try: optimize the whole matrix
    if (tryDirect) {
        for(size_t i=0; i<randomloops; ++i) {
            Matrix lM(QQ,M.rowdim(),M.coldim()); sparse2sparse(lM,M);
            Matrix lT(QQ,T.rowdim(),T.coldim()); sparse2sparse(lT,T);
            std::ostringstream lssout;
                // Cancellation-free optimization
            input2Temps(lssout, lM.coldim(), 'i', 't', lT);
            auto lnbops( Optimizer(lssout, lM, 'i', 'o', 't', 'r') );
#ifdef VERBATIM_PARSING
            std::clog << "# Found, direct: " << lnbops.first << "\tadditions, "
                       << lnbops.second << "\tmultiplications." << std::endl;
#endif
            if ( (ssout.tellp() == std::streampos(0)) ||
                 (lnbops.first<nbops.first) ||
                 ( (lnbops.first==nbops.first) && (lnbops.second<nbops.second) ) ) {
                ssout.clear(); ssout.str(std::string());
                ssout << lssout.str();
                nbops = lnbops;
            }
        }
    }

        // ============================================================
        // Second try: separate independent and dependent rows
    if (tryKernel) {
        for(size_t i=0; i<randomloops; ++i) {
            Matrix lT(QQ,T.rowdim(),T.coldim()); sparse2sparse(lT,T);
            std::ostringstream lssout;
            Matrix NullSpace(QQ,lT.coldim(),T.coldim());
            auto lnbops( nullspacedecomp(lssout, NullSpace, lT) );
#ifdef VERBATIM_PARSING
            std::clog << "# Found, kernel: " << lnbops.first << "\tadditions, "
                       << lnbops.second << "\tmultiplications." << std::endl;
#endif
            if ( (ssout.tellp() == std::streampos(0)) ||
                 (lnbops.first<nbops.first) ||
                 ( (lnbops.first==nbops.first) && (lnbops.second<nbops.second) ) ) {
                ssout.clear(); ssout.str(std::string());
                ssout << lssout.str();
                nbops = lnbops;
            }
        }

        if (nbops == std::pair<size_t,size_t>{-1,-1}) {
                // Zero dimensional kernel
            std::cerr << "# Fail: zero dimensional kernel.\n"
                      << "# --> try direct or transposing." << std::endl;
            nbops = std::pair<size_t,size_t>{0,0};
        }
    }

    chrono.stop();

    std::cout << ssout.str() << std::flush;

    if ((nbops.first !=0 || nbops.second != 0)) {
        std::clog << std::string(40,'#') << std::endl;
        std::clog << "# " << nbops.first << "\tadditions\tinstead of " << addinit
                  << ' ' << chrono << std::endl;
        std::clog << "# " << nbops.second << "\tmultiplications\tinstead of " << mulinit << std::endl;
        std::clog << std::string(40,'#') << std::endl;
    }

    return 0;
}

#ifdef RANDOM_TIES
#  define DORANDOMSEARCH true
#else
#  define DORANDOMSEARCH false
#endif


// ============================================================
// Main: select between file / std::cin
// -D/-K options seect direct/kernel methods only (default is both)
// -P/-M option choose the printing format
// -O # search for reduced number of additions, then multiplications
//      i.e. min of random # tries (requires RANDOM_TIES)
int main(int argc, char** argv) {

        // ============================================================
        // Linear Transformation
    bool printMaple(false),
        printPretty(false),
        directOnly(false),
        kernelOnly(false);
    size_t randomloops(DEFAULT_RANDOM_LOOPS);

    std::string filename;

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0] << " [-h|-M|-P|-K|-D] [stdin|matrixfile.sms]\n";
            exit(-1);
        } else if (args == "-M") { printMaple = true; }
        else if (args == "-P") { printPretty = true; }
        else if (args == "-D") { directOnly = true; }
        else if (args == "-K") { kernelOnly = true; }
        else if (args == "-O") {
            randomloops = atoi(argv[++i]);
            if ( (randomloops>1) && (!DORANDOMSEARCH) ) {
                randomloops = 1;
                std::cerr << "# WARNING: RANDOM_TIES not defined,"
                          << " random loops disabled." << std::endl;
            }
        }
        else { filename = args; }
    }

    const bool tryDirect(directOnly || !kernelOnly);
    const bool tryKernel(kernelOnly || !directOnly);

    QRat QQ;

    if (filename == "") {
        return Selector(std::cin, randomloops, printMaple, printPretty, tryDirect, tryKernel);
    } else {
        std::ifstream inputmatrix(filename);
        if ( inputmatrix ) {
            int rt=Selector(inputmatrix, randomloops, printMaple, printPretty, tryDirect, tryKernel);
            inputmatrix.close();
            return rt;
        }
    }
    return -1;
}
// ============================================================
