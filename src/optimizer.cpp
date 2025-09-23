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

#ifdef OPTIMIZE_ADDITIONS
// Optimize for additions first, then multiplications
auto cmpOpCount {[](const auto& a, const auto& b) { return (a.first<b.first) || ( (a.first==b.first) && (a.second<b.second) ); } };
#else
// Optimize for sum of additions and multiplications
auto cmpOpCount {[](const auto& a, const auto& b) { return (a.first+a.second<b.first+b.second); } };
#endif


namespace PLinOpt {


// ============================================================
// Optimizing a linear program (Direct CSE or Kernel methods)
template<typename Field>
Pair<size_t>& DKOptimiser(Pair<size_t>& nbops, std::ostringstream& ssout,
                          const Field& F, const Matrix& M, const Matrix& T,
                          Givaro::Timer& global, const size_t randomloops,
                          const bool printMaple, const bool printPretty,
                          const bool tryDirect, const bool tryKernel,
                          const bool mostCSE, const bool allkernels) {

    Givaro::Timer chrono; chrono.start();
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
            auto lnbops( Optimizer(lssout, lM, 'o', 't', 'r') );


#pragma omp critical
            {
#ifdef VERBATIM_PARSING
                std::clog << "# Found, direct: "
                          << lnbops.first << "\tadditions, "
                          << lnbops.second << "\tmultiplications." << std::endl;
#endif
                if ( (ssout.tellp() == std::streampos(0))
                     || cmpOpCount(lnbops,nbops) ) {
                    std::clog << "# Found D: "
                              << lnbops.first << '|' << lnbops.second
                              << " instead of "
                              << nbops.first << '|' << nbops.second
                              << std::endl;
                    ssout.clear(); ssout.str(std::string());
                    ssout << lssout.str();
                    nbops = lnbops;
                }
            }
        }
    }

        // ============================================================
        // Second try: separate independent and dependent rows
    if (tryKernel) {
#pragma omp parallel for shared(T,ssout,nbops)
        for(size_t i=0; i<randomloops; ++i) {
            FMatrix lT(T, F);
            std::ostringstream lssout;
            FMatrix NullSpace(F,lT.coldim(),T.coldim());
            auto lnbops( nullspacedecomp(lssout, NullSpace, lT) );

#pragma omp critical
            {
#ifdef VERBATIM_PARSING
                std::clog << "# Found, kernel: "
                          << lnbops.first << "\tadditions, "
                          << lnbops.second << "\tmultiplications." << std::endl;
#endif
                if ( (ssout.tellp() == std::streampos(0) )
                     || cmpOpCount(lnbops,nbops) ) {
                    ssout.clear(); ssout.str(std::string());
                    ssout << lssout.str();
                    if (cmpOpCount(lnbops,nbops)) {
                        std::clog << "# Found K: "
                                  << lnbops.first << '|' << lnbops.second
                                  << " instead of "
                                  << nbops.first << '|' << nbops.second
                                  << std::endl;
                        nbops = lnbops;
                    }
                }
            }
        }

        if (nbops == Pair<size_t>{-1,-1}) {
                // Zero dimensional kernel
            std::cerr << "# \033[1;36mFail: zero dimensional kernel.\033[0m\n"
                      << "# --> try direct or transposing." << std::endl;
            nbops = Pair<size_t>{0,0};
        }
    }

    chrono.stop(); global += chrono;

        // ============================================================
        // Exhaustive nullspace permutation search (if # is <= 12!)
    if (allkernels && (M.rowdim() < 13)) {
        chrono.clear(); chrono.start();

        const size_t m(M.rowdim());
        std::vector<size_t> Fm { factorial(m) };
        std::ostringstream kout;

        auto knbops(nbops);

#pragma omp parallel for shared(Fm,T,kout,knbops)
        for(size_t i=0; i<Fm.back(); ++i) {
            std::vector<size_t> l{kthpermutation(i,m,Fm)};
            FMatrix lT(T, F);
            std::ostringstream lkout;
            FMatrix NullSpace(F,lT.coldim(),T.coldim());
            auto lkops( nullspacedecomp(lkout, NullSpace, lT, l, mostCSE) );

#pragma omp critical
            {
#ifdef VERBATIM_PARSING
                std::clog << "# Found, kernel: "
                          << lkops.first << "\tadditions, "
                          << lkops.second << "\tmultiplications." << std::endl;
#endif
                if ( (kout.tellp() == std::streampos(0))
                     || cmpOpCount(lkops,knbops) ) {
                    std::clog << "# Found E: "
                              << lkops.first << '|' << lkops.second
                              << " instead of "
                              << knbops.first << '|' << knbops.second
                              << std::endl;
                    kout.clear(); kout.str(std::string());
                    kout << lkout.str();
                    knbops = lkops;
                }
            }
        }

        chrono.stop(); global += chrono;

        if (cmpOpCount(knbops,nbops)) {
            std::clog << "# \033[1;36m"
                      << "Exhaustive kernel permutation, found:"
                      << "\033[0m" << std::endl;

// std::clog << kout.str() << std::flush;

            ssout.str(kout.str());

            if ((knbops.first !=0 || knbops.second != 0)) {
                std::clog << std::string(40,'#') << std::endl
                          << "# \033[1;32m" << knbops.first
                          << "\tadditions\tinstead of " << nbops.first
                          << "\033[0m \t" << chrono << std::endl
                          << "# \033[1;32m" << knbops.second
                          << "\tmultiplications\tinstead of " << nbops.second
                          << "\033[0m" << std::endl
                          << std::string(40,'#') << std::endl;
            }
            nbops = knbops;
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
        auto rnbops( RecOptimizer(iout, lM, 'o', 't', 'r') );

        chrono.stop(); global += chrono;

        if (cmpOpCount(rnbops, nbops)) {
            std::clog << "# \033[1;36m"
                      << "Exhaustive greedy CSE search, found:"
                      << "\033[0m" << std::endl;

            ssout.str(iout.str());

            if ((rnbops.first !=0 || rnbops.second != 0)) {
                std::clog << std::string(40,'#') << std::endl;
                std::clog << "# \033[1;32m" << rnbops.first << "\tadditions\tinstead of " << nbops.first
                          << "\033[0m \t" << chrono << std::endl;
                std::clog << "# \033[1;32m" << rnbops.second << "\tmultiplications\tinstead of " << nbops.second << "\033[0m" << std::endl;
                std::clog << std::string(40,'#') << std::endl;
                nbops = rnbops;
            }
        } else {
                // Optimal was already found
            std::clog << "# \033[1;36m"
                      << "No greedy CSE schedule has less additions.\033[0m \t"
                      << chrono << std::endl;
        }
    }

    return nbops;
}

#include <linbox/algorithms/gauss.h>

// ============================================================
// Optimizing a linear program (Gaussian elimination method)
template<typename outstream, typename Field>
Pair<size_t>& LUOptimiser(Pair<size_t>& nbops, outstream& gout,
                          const Field& F, const Matrix& M, const Matrix& T,
                          const size_t randomloops, const bool mostCSE) {

        // ============================================================
        // Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;

        // ============================================================

    FMatrix U(M, F);

#ifdef DEBUG
    U.write(std::clog << "## M:=Matrix(", FileFormat(8)) << ");" << std::endl;
#endif


    typename Field::Element Det;
    size_t Rank;
    FMatrix L(F, U.rowdim(), U.rowdim());
    LinBox::Permutation<Field> Q(F,U.rowdim());
    LinBox::Permutation<Field> P(F,U.coldim());
    LinBox::GaussDomain<Field> GD(F);
    GD.QLUPin(Rank, Det, Q, L, U, P, U.rowdim(), U.coldim() );

#ifdef DEBUG
    Q.write(std::clog << "## Q:=Matrix(", FileFormat(8)) << ");" << std::endl;
    L.write(std::clog << "## L:=Matrix(", FileFormat(8)) << ");" << std::endl;
    U.write(std::clog << "## U:=Matrix(", FileFormat(8)) << ");" << std::endl;
    P.write(std::clog << "## P:=Matrix(", FileFormat(8)) << ");" << std::endl;
#endif

// #pragma omp parallel for shared(Q,L,U,P,gout,nbops)
//     for(size_t i=0; i<randomloops; ++i) {
        std::ostringstream luout;

            // ============================================
            // Applying permutation P to 'i' variables
        input2Temps(luout, P.rowdim(), 'i', 't', P.getStorage());

            // ============================================
            // Applying Upper matrix to 't' variables
        auto Uops = Optimizer(luout, U, 'v', 't', 'r');

            // ============================================
            // Applying Lower matrix to 'v' variables
        auto Lops = Optimizer(luout, L, 'x', 'v', 'g');
            // ============================================
            // Applying permutation Q back into 'o' variables
        input2Temps(luout, Q.rowdim(), 'x', 'o', Q.getStorage());

        Uops.first += Lops.first;
        Uops.second += Lops.second;

// #pragma omp critical
//         {
            if ( (gout.tellp() == std::streampos(0))
                 || cmpOpCount(Uops,nbops) ) {
                std::clog << "# Found G: "
                          << Uops.first << '|' << Uops.second
                          << " instead of "
                          << nbops.first << '|' << nbops.second
                          << std::endl;
                gout.clear(); gout.str(std::string());
                gout.str(luout.str());
                nbops = Uops;
            }
//         }
//     }

    return nbops;
}

// ============================================================
// Optimizing a linear program over a field
template<typename Field>
int FOptimiser(std::istream& input, const size_t randomloops,
               const bool printMaple, const bool printPretty,
               const bool tryDirect, const bool tryKernel,
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


        // ============================================================
        // Compute naive number of operations

    std::ostringstream ssout;
    const Pair<size_t> opsinit(naiveOps(M));
    Pair<size_t> nbops(opsinit);

    DKOptimiser(nbops, ssout, F, M, T, global, randomloops, printMaple,
                printPretty, tryDirect, tryKernel, mostCSE, allkernels);

    std::ostringstream luout;
    LUOptimiser(nbops, luout, F, M, T, randomloops, mostCSE);

    std::cout << ssout.str() << std::flush;

    if ((nbops.first !=0 || nbops.second != 0)) {
        std::clog << std::string(40,'#') << std::endl;
        std::clog << "# \033[1;32m" << nbops.first << "\tadditions\tinstead of " << opsinit.first
                  << "\033[0m \t" << global << std::endl;
        std::clog << "# \033[1;32m" << nbops.second << "\tmultiplications\tinstead of " << opsinit.second << "\033[0m" << std::endl;
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








// ============================================================
// Choice between modular computation or over the rationals
int Selector(std::istream& input, const size_t randomloops,
             const bool printMaple, const bool printPretty,
             const bool tryDirect, const bool tryKernel,
             const bool mostCSE, const bool allkernels,
             const Givaro::Integer& q) {
    if (! Givaro::isZero(q)) {
        Givaro::Modular<Givaro::Integer> FF(q);
        return FOptimiser(input, randomloops, printMaple, printPretty,
                          tryDirect, tryKernel, mostCSE, allkernels, FF);
    } else {
        QRat QQ;
        return FOptimiser(input, randomloops, printMaple, printPretty,
                          tryDirect, tryKernel, mostCSE, allkernels, QQ);
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
                << "  -O #: search for reduced number of additions and/then multiplications\n";

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
        return PLinOpt::Selector(std::cin, randomloops, printMaple, printPretty,
                                 tryDirect, tryKernel, mostCSE, allkernels, prime);
    } else {
        std::ifstream inputmatrix(filename);
        if ( inputmatrix ) {
            int s=PLinOpt::Selector(inputmatrix, randomloops, printMaple, printPretty,
                           tryDirect, tryKernel, mostCSE, allkernels, prime);
            inputmatrix.close();
            return s;
        }
    }
    return -1;
}
// ============================================================
