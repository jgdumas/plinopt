// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library
 * Reference:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *     Strassen's algorithm is not optimally accurate
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/

#include "plinopt_sparsify.h"

/****************************************************************
 * Sparse decomposition of matrices
 * Reference:
 *   [ Gal Beniamini, Oded Schwartz;
 *     Faster Matrix Multiplication via Sparse Decomposition.
 *     SPAA 2019: 11-22]
 ****************************************************************/

// ============================================================
int Factorizer(std::istream& input, const FileFormat& matformat,
               const size_t selectinnerdim, const size_t randomloops) {

        // ============================================================
        // Read Matrix of Linear Transformation
    QRat QQ;
    QMstream ms(QQ,input);
    Matrix M(ms); M.resize(M.rowdim(),M.coldim());
    size_t sc,sb,sr;
    densityProfile(std::clog << "# Initial profile: ", sc, M)
                             << std::endl;

#ifdef VERBATIM_PARSING
    M.write(std::clog,FileFormat::Pretty) << std::endl;
#endif
    size_t innerdim(selectinnerdim == 0 ? M.coldim() : selectinnerdim);
    if ( (innerdim >= M.rowdim()) ||
         (innerdim < M.coldim()) ) {
        std::cerr << "# \033[1;36mFail: inner dimension has to be between " << M.coldim() << " and " << (M.rowdim()-1) << ".\033[0m\n";
        return -1;
    }

    Givaro::Timer elapsed;
    Matrix CoB(QQ, innerdim, M.coldim());
    Matrix Res(QQ, M.rowdim(), innerdim);
    Pair<size_t> nbops{backSolver(CoB, Res, M, QQ)};

    elapsed.start();
#pragma omp parallel for shared(Res,CoB,M,QQ,nbops,innerdim)
    for(size_t i=0; i<randomloops; ++i) {
        Matrix lCoB(QQ, innerdim, M.coldim());
        Matrix lRes(QQ, M.rowdim(), innerdim);
        auto bSops{backSolver(lCoB, lRes, M, QQ)};
#ifdef VERBATIM_PARSING
        std::clog << "# Res/CoB profile[" << i << "]: " << bSops << std::endl;
#endif
        if ( (bSops.first<nbops.first) ||
             ( (bSops.first==nbops.first)
               && (bSops.second<nbops.second) ) ) {
            nbops = bSops;
            std::clog << "# Found [" << i << "], R/CB profile: " << bSops << std::endl;
            matrixCopy(CoB, lCoB, QQ);
            matrixCopy(Res, lRes, QQ);
        }
    }
    elapsed.stop();

        // ============================================================
        // Print resulting matrices

        // change of basis to stdout
    densityProfile(std::clog << "# Alternate basis profile: \033[1;36m",
                   sb, CoB) << "\033[0m" << std::endl;
    CoB.write(std::cout, matformat) << std::endl;


        // residuum sparse matrix to stdlog
    densityProfile(std::clog << "# Sparse residuum profile: \033[1;36m",
                   sr, Res) << "\033[0m" << std::endl;
    Res.write(std::clog, matformat)<< std::endl;


        // Final check that we computed a factorization M=Res.CoB
    std::clog << std::string(30,'#') << std::endl;
    consistency(std::clog, M, Res, CoB)
        << " \033[1;36m" << Res.rowdim() << 'x' << Res.coldim()
        << " by " << CoB.rowdim() << 'x' << CoB.coldim() << " with "
        << sr << " non-zeroes (" << sb << " alt.) instead of " << sc
        << "\033[0m:" << ' ' << elapsed << std::endl;

    return 0;
}


// ============================================================
// Main: select between file / std::cin
//       -k #: sets the inner dimension (default is column dimension)
//       -M/-P/-S/-L: selects the ouput format
//       -O # search for reduced randomized sparsity
//      i.e. min of random # tries (requires definition of RANDOM_TIES)
int main(int argc, char** argv) {

    FileFormat matformat = FileFormat::Pretty;
    std::string filename;
    size_t innerdim(0);					// will default to columndimension
    size_t randomloops(DORANDOMSEARCH?DEFAULT_RANDOM_LOOPS:1);

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << " [-h|-M|-P|-S|-L|-k #|-O #] [stdin|matrixfile.sms]\n";

            std::clog
                << "  -k #: inner dimension (default is column dimension)\n"
                << "  -M/-P/-S/-L: selects the ouput format\n"
                << "  -O #: search for reduced randomized sparsity\n";

            exit(-1);
        }
        else if (args == "-M") { matformat = FileFormat(1); } // Maple
        else if (args == "-S") { matformat = FileFormat(5); } // SMS
        else if (args == "-P") { matformat = FileFormat(8); } // Pretty
        else if (args == "-L") { matformat = FileFormat(12); }// Linalg
        else if (args == "-k") { innerdim = atoi(argv[++i]); }
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

    if (filename == "") {
        return Factorizer(std::cin, matformat, innerdim, randomloops);
    } else {
        std::ifstream inputmatrix(filename);
        return Factorizer(inputmatrix, matformat, innerdim, randomloops);
    }

    return -1;
}
// ============================================================
