// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library
 * Reference:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
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

namespace PLinOpt {
// ============================================================
template<typename _Mat>
int TFactorizer(const _Mat& M, const FileFormat& matformat,
                const size_t selectinnerdim, const size_t randomloops) {



    using FMatrix = _Mat;
    using Field = typename _Mat::Field;
    const Field& F(M.field());

    size_t sc,sb,sr;
    densityProfile(std::clog << "# Initial profile: ", sc, M)
                             << std::endl;

    FMatrix Alt(F, M.rowdim(), selectinnerdim);
    FMatrix CoB(F, selectinnerdim, M.coldim());

    Givaro::Timer elapsed; elapsed.start();
    Factorizer(Alt, CoB, M, randomloops, selectinnerdim);
    elapsed.stop();

        // ============================================================
        // Print resulting matrices

        // change of basis to stdout
    densityProfile(std::clog << "# Alternate basis profile: \033[1;36m",
                   sb, CoB) << "\033[0m" << std::endl;
    CoB.write(std::cout, matformat) << std::endl;

        // residuum sparse matrix to stdlog
    densityProfile(std::clog << "# Sparse residuum profile: \033[1;36m",
                   sr, Alt) << "\033[0m" << std::endl;
    Alt.write(std::clog, matformat)<< std::endl;

        // Final check that we computed a factorization M=Alt.CoB
    std::clog << std::string(30,'#') << std::endl;
    consistency(std::clog, M, Alt, CoB)
        << " \033[1;36m" << Alt.rowdim() << 'x' << Alt.coldim()
        << " by " << CoB.rowdim() << 'x' << CoB.coldim() << " with "
        << sr << " non-zeroes (" << sb << " alt.) instead of " << sc
        << "\033[0m:" << ' ' << elapsed << std::endl;

    return 0;
}



} // End of namespace PLinOpt
// ============================================


// ============================================================
// Checking a linear program with a matrix
int Fmain(std::istream& input, const Givaro::Integer& q,
          const PLinOpt::FileFormat& mformat,
          const size_t innerdim, const size_t loops){
        // ============================================================
        // Read Matrix of Linear Transformation
    PLinOpt::QRat QQ;
    PLinOpt::QMstream ms(QQ,input);
    PLinOpt::Matrix rM(ms); rM.resize(rM.rowdim(),rM.coldim());

    if (! Givaro::isZero(q)) {
            // ============================================================
            // Rebind matrix type over sub field matrix type
        using Field = Givaro::Modular<Givaro::Integer>;
        using FMatrix=typename PLinOpt::Matrix::template rebind<Field>::other;

        const Field FF(q);
        FMatrix fM(rM, FF);
        return PLinOpt::TFactorizer(fM, mformat, innerdim, loops);
    } else {
        return PLinOpt::TFactorizer(rM, mformat, innerdim, loops);
    }
}

// ============================================================
// Main: select between file / std::cin
//       -k #: sets the inner dimension (default is column dimension)
//       -M/-P/-S/-L: selects the ouput format
//       -O # search for reduced randomized sparsity
//      i.e. min of random # tries (requires definition of RANDOM_TIES)
int main(int argc, char** argv) {
    using PLinOpt::FileFormat;

    FileFormat matformat = FileFormat::Pretty;
    std::string filename;
    size_t innerdim(0);					// will default to columndimension
    size_t randomloops(DORANDOMSEARCH?DEFAULT_RANDOM_LOOPS:1);
    Givaro::Integer q(0u);

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog
                << "Usage: " << argv[0]
                << " [-h|-M|-P|-S|-L|-k #|-O #] [stdin|matrixfile.sms]\n"
                << "  -k #: inner dimension (default is column dimension)\n"
                << "  -M/-P/-S/-L: selects the ouput format\n"
                << "  -q #: search modulo (default is Rationals)\n"
                << "  -O #: search for reduced randomized sparsity (default "
                << randomloops << " loops)\n";

            exit(-1);
        }
        else if (args == "-M") { matformat = FileFormat(1); } // Maple
        else if (args == "-S") { matformat = FileFormat(5); } // SMS
        else if (args == "-P") { matformat = FileFormat(8); } // Pretty
        else if (args == "-L") { matformat = FileFormat(12); }// Linalg
        else if (args == "-k") { innerdim = atoi(argv[++i]); }
        else if (args == "-q") { q = Givaro::Integer(argv[++i]); }
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
        return Fmain(std::cin, q, matformat, innerdim, randomloops);
    } else {
        std::ifstream inputmatrix(filename);
        return Fmain(inputmatrix, q, matformat, innerdim, randomloops);
    }
}
// ============================================================
