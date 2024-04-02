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
template<typename Field>
int DKOptimiser(std::istream& input, const size_t randomloops,
                const bool printMaple, const bool printPretty,
                const bool tryDirect, const bool tryKernel,
                const bool mostCSE, const bool allkernels, const Field& F) {

        // ============================================================
        // Read Matrix of Linear Transformation
    QRat QQ;
    QMstream ms(QQ,input);
    Matrix M(ms); M.resize(M.rowdim(),M.coldim());

    Givaro::Timer chrono, global;
    chrono.start();

    Matrix T(QQ,M.coldim(),M.rowdim()); Transpose(T,M);

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
    Pair<size_t> nbops{addinit,mulinit};

        // ============================================================
        // Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;

        // ============================================================
        // First try: optimize the whole matrix
    if (tryDirect) {
#pragma omp parallel for shared(M,T,ssout,nbops)
        for(size_t i=0; i<randomloops; ++i) {
            FMatrix lM(M, F);
            FMatrix lT(T, F);

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
#pragma omp parallel for shared(T,ssout,nbops)
        for(size_t i=0; i<randomloops; ++i) {
            FMatrix lT(T, F);
//             Matrix lT(QQ,T.rowdim(),T.coldim()); matrixCopy(lT,T,QQ);
            std::ostringstream lssout;
            FMatrix NullSpace(F,lT.coldim(),T.coldim());
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

        if (nbops == Pair<size_t>{-1,-1}) {
                // Zero dimensional kernel
            std::cerr << "# \033[1;36mFail: zero dimensional kernel.\033[0m\n"
                      << "# --> try direct or transposing." << std::endl;
            nbops = Pair<size_t>{0,0};
        }
    }

    chrono.stop(); global = chrono;
// std::clog << ssout.str() << std::flush;

        // ============================================================
        // Exhaustive nullspace permutation search (if # is <= 12!)
    if (allkernels && (M.rowdim() < 13)) {
        chrono.clear(); chrono.start();

        const size_t m(M.rowdim());
        std::vector<long> Fm { factorial(m) };
        std::ostringstream kout;

        auto knbops(nbops);

#pragma omp parallel for shared(Fm,T,kout,knbops)
        for(size_t i=0; i<Fm.back(); ++i) {
            std::vector<size_t> l{kthpermutation(i,m,Fm)};
            FMatrix lT(T, F);
            std::ostringstream lkout;
            FMatrix NullSpace(F,lT.coldim(),T.coldim());
            auto lkops( nullspacedecomp(lkout, NullSpace, lT, l, mostCSE) );
#ifdef VERBATIM_PARSING
            std::clog << "# Found, kernel: " << lkops.first << "\tadditions, "
                       << lkops.second << "\tmultiplications." << std::endl;
#endif
            if ( (kout.tellp() == std::streampos(0)) ||
                 (lkops.first<knbops.first) ||
                 ((lkops.first==knbops.first)&&(lkops.second<knbops.second))) {
                kout.clear(); kout.str(std::string());
                kout << lkout.str();
                knbops = lkops;
            }
        }

        chrono.stop(); global += chrono;

        if ( (knbops.first < nbops.first) ||
             ( (knbops.first == nbops.first) && (knbops.second < nbops.second) ) ) {
            nbops = knbops;
            std::clog << "# \033[1;36m"
                      << "Exhaustive kernel permutation, found:"
                      << "\033[0m" << std::endl;

// std::clog << kout.str() << std::flush;

            ssout.str(kout.str());

            if ((knbops.first !=0 || knbops.second != 0)) {
                std::clog << std::string(40,'#') << std::endl;
                std::clog << "# \033[1;32m" << knbops.first << "\tadditions\tinstead of " << addinit
                          << "\033[0m \t" << chrono << std::endl;
                std::clog << "# \033[1;32m" << knbops.second << "\tmultiplications\tinstead of " << mulinit << "\033[0m" << std::endl;
                std::clog << std::string(40,'#') << std::endl;
            }
        } else {
                // Optimal was already found
            std::clog << "# \033[1;36m"
                      << "No kernel permutation has less additions.\033[0m \t"
                      << chrono << std::endl;
        }
    }

        // ============================================================
        // MostCSE Greedy CSE search
    if (mostCSE) {
        chrono.clear(); chrono.start();

        std::vector<std::ostringstream> out;
        FMatrix lM(M, F);
        FMatrix lT(T, F);
        std::ostringstream iout;
            // Cancellation-free optimization
        input2Temps(iout, lM.coldim(), 'i', 't', lT);
        auto rnbops( RecOptimizer(iout, lM, 'i', 'o', 't', 'r') );

        chrono.stop(); global += chrono;

        if ( (rnbops.first < nbops.first) ||
             ( (rnbops.first == nbops.first) && (rnbops.second < nbops.second) ) ) {
            std::clog << "# \033[1;36m"
                      << "Exhaustive greedy CSE search, found:"
                      << "\033[0m" << std::endl;

            ssout.str(iout.str());

// std::clog << iout.str() << std::flush;


            if ((rnbops.first !=0 || rnbops.second != 0)) {
                std::clog << std::string(40,'#') << std::endl;
                std::clog << "# \033[1;32m" << rnbops.first << "\tadditions\tinstead of " << addinit
                          << "\033[0m \t" << chrono << std::endl;
                std::clog << "# \033[1;32m" << rnbops.second << "\tmultiplications\tinstead of " << mulinit << "\033[0m" << std::endl;
                std::clog << std::string(40,'#') << std::endl;
            }
        } else {
                // Optimal was already found
            std::clog << "# \033[1;36m"
                      << "No CSE scheduling has less additions.\033[0m \t"
                      << chrono << std::endl;
        }
    }


    std::cout << ssout.str() << std::flush;


    if ((nbops.first !=0 || nbops.second != 0)) {
        std::clog << std::string(40,'#') << std::endl;
        std::clog << "# \033[1;32m" << nbops.first << "\tadditions\tinstead of " << addinit
                  << "\033[0m \t" << global << std::endl;
        std::clog << "# \033[1;32m" << nbops.second << "\tmultiplications\tinstead of " << mulinit << "\033[0m" << std::endl;
        std::clog << std::string(40,'#') << std::endl;
    }


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
// ====================================================================



// ============================================================
// Choice between modular computation or over the rationals
int Selector(std::istream& input, const size_t randomloops,
             const bool printMaple, const bool printPretty,
             const bool tryDirect, const bool tryKernel,
             const bool mostCSE, const bool allkernels,
             const Givaro::Integer& q) {
    if (! Givaro::isZero(q)) {
        Givaro::Modular<Givaro::Integer> FF(q);
        return DKOptimiser(input, randomloops, printMaple, printPretty,
                           tryDirect, tryKernel, mostCSE, allkernels, FF);
    } else {
        QRat QQ;
        return DKOptimiser(input, randomloops, printMaple, printPretty,
                           tryDirect, tryKernel, mostCSE, allkernels, QQ);
    }
}
// ============================================================


// ============================================================
// Main: select between file / std::cin
// -D/-K options select direct/kernel methods only (default is both)
// -P/-M option choose the printing format
// -E option for exhaustive greedy CSE search
// -N option for exhaustive nullspace permutation search
// -q # search for linear map modulo # (default is over the rationals)
// -O # search for reduced number of additions, then multiplications
//      i.e. min of random # tries (requires definition of RANDOM_TIES)
int main(int argc, char** argv) {

        // ============================================================
        // Linear Transformation
    bool printMaple(false), printPretty(false), mostCSE(false),
        directOnly(false), kernelOnly(false), allkernels(false);
    size_t randomloops(DORANDOMSEARCH?DEFAULT_RANDOM_LOOPS:1);
    Givaro::Integer prime(0u);

    std::string filename;

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << " [-h|-M|-P|-K|-D|-E|-N|-q #|-O #] [stdin|matrixfile.sms]\n";

            std::clog
                << "  -D/-K: select direct/kernel methods only (default is both)\n"
                << "  -E: exhaustive greedy CSE search (default is not)\n"
                << "  -N: exhaustive nullspace permutations (default is not)\n"
                << "  -P/-M: choose the printing format\n"
                << "  -q #: search modulo (default is Rationals)\n"
                << "  -O #: search for reduced number of additions, then multiplications\n";

            exit(-1);
        }
        else if (args == "-M") { printMaple = true; }
        else if (args == "-P") { printPretty = true; }
        else if (args == "-D") { directOnly = true; }
        else if (args == "-K") { kernelOnly = true; }
        else if (args == "-E") { mostCSE = true; }
        else if (args == "-N") { allkernels = true; }
        else if (args == "-q") { prime = Givaro::Integer(argv[++i]); }
        else if (args == "-O") { randomloops = atoi(argv[++i]);
            if ( (randomloops>1) && (!DORANDOMSEARCH) ) {
                randomloops = 1;
                std::cerr << "#  \033[1;36mWARNING: RANDOM_TIES not defined,"
                          << " random loops disabled\033[0m." << std::endl;
            }
        }
        else { filename = args; }
    }

    const bool tryDirect(directOnly || !kernelOnly);
    const bool tryKernel(kernelOnly || !directOnly);

    if (filename == "") {
        return Selector(std::cin, randomloops, printMaple, printPretty,
                        tryDirect, tryKernel, mostCSE, allkernels, prime);
    } else {
        std::ifstream inputmatrix(filename);
        if ( inputmatrix ) {
            int s=Selector(inputmatrix, randomloops, printMaple, printPretty,
                           tryDirect, tryKernel, mostCSE, allkernels, prime);
            inputmatrix.close();
            return s;
        }
    }
    return -1;
}
// ============================================================
