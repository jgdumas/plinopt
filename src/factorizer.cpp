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
int TFactorizer(const _Mat& A, const FileFormat& matformat,
                const size_t selectinnerdim, const size_t randomloops,
                const size_t blocksize, const size_t maxnumcoeff,
                const bool initialElimination,
                const bool initialSparsification) {



    using FMatrix = _Mat;
    using Field = typename _Mat::Field;
    using DenseFMatrix = LinBox::DenseMatrix<Field>;
    const Field& F(A.field());

    size_t sc;
    densityProfile(std::clog << "# [FCTZ] Initial profile: ", sc, A)
                             << std::endl;

    Givaro::Timer elapsed;
    FMatrix Alt(F, A.rowdim(), selectinnerdim);
    FMatrix CoB(F, selectinnerdim, A.coldim());

    if (initialSparsification) {
            // ============================================================
            // Sparsifier
        FMatrix Cs(F,A.coldim(), A.coldim());
        FMatrix M(F,A.rowdim(), A.coldim());
        blockSparsifier(elapsed, Cs, M, A, blocksize,
                        matformat, maxnumcoeff, initialElimination);

        const std::string blsp{"[SPRB]"};
        size_t sm(profileConsistency(matformat, elapsed,
                                     A, sc, blsp, M, 0, Cs, 0));
        std::clog << std::string(30,'#') << std::endl;

            // ============================================================
            // Factorizer
        FMatrix Ca(F, selectinnerdim, M.coldim());
        Givaro::Timer Felapse; Felapse.start();
        Factorizer(Alt, Ca, M, randomloops, selectinnerdim);
        Felapse.stop();
        elapsed += Felapse;

        const std::string sprf{"[FCTA]"};
        profileConsistency(matformat, elapsed, M, sm, sprf, Alt, 0, Ca, 0);
        std::clog << std::string(30,'#') << std::endl;

            // ============================================================
            // Recombination, from A = M.Cs ; M = Alt.Ca:
            // CoB = Ca.Cs, so that A = Alt.CoB
        DenseFMatrix Do(F, Ca.rowdim(), Cs.coldim());
        LinBox::MatrixDomain<Field> BMD(F);

        Felapse.clear(); Felapse.start();
        BMD.mul(Do, Ca, Cs);
        CoB.resize(Do.rowdim(), Do.coldim());
        dense2sparse(CoB, Do);
        Felapse.stop();
        elapsed += Felapse;
    } else {
        elapsed.start();
        Factorizer(Alt, CoB, A, randomloops, selectinnerdim);
        elapsed.stop();
    }

        // ============================================================
        // Print resulting matrices
    const std::string fctz{"[FCTZ]"};
    profileConsistency(matformat, elapsed, A, sc, fctz, Alt, -1, CoB, 1);
    return 0;
}



} // End of namespace PLinOpt
// ============================================


// ============================================================
// Checking a linear program with a matrix
int Fmain(std::istream& input, const Givaro::Integer& q,
          const PLinOpt::FileFormat& mformat,
          const size_t innerdim, const size_t loops,
          const size_t blocksize, const size_t maxnumcoeff,
          const bool initialElimination, const bool initialSparsification) {
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
        return PLinOpt::TFactorizer(fM, mformat, innerdim, loops,
                                    blocksize, maxnumcoeff,
                                    initialElimination, initialSparsification);
    } else {
        return PLinOpt::TFactorizer(rM, mformat, innerdim, loops,
                                    blocksize, maxnumcoeff,
                                    initialElimination, initialSparsification);
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
    size_t maxnumcoeff(COEFFICIENT_SEARCH); // default max coefficients
    size_t blocksize(4u);                   // default column block size
    bool initialSparsification(false);
    bool initialElimination(true);
    Givaro::Integer q(0u);

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog
                << "Usage: " << argv[0]
                << " [-h|-M|-P|-S|-L|[-k|-O|-c|-b|-U|-V #]] [stdin|matrixfile.sms]\n"
                << "  -k #: inner dimension (default is column dimension)\n"
                << "  -M/-P/-S/-L: selects the ouput format\n"
                << "  -V [1|0]: initial sparsification or not (default 0)\n"
                << "  -b #: states the blocking dimension (default "
                << blocksize << ")\n"
                << "  -c #: max number of coefficients per iteration (default "
                << maxnumcoeff << ")\n"
                << "  -U [1|0]: initial LU factorization or not (default 1) \n"
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
        else if (args == "-V") { initialSparsification = atoi(argv[++i]); }
        else if (args == "-b") { blocksize = atoi(argv[++i]); }
        else if (args == "-c") { maxnumcoeff = atoi(argv[++i]); }
        else if (args == "-U") { initialElimination = atoi(argv[++i]); }
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
        return Fmain(std::cin, q, matformat, innerdim, randomloops,
                     blocksize, maxnumcoeff,
                     initialElimination, initialSparsification);
    } else {
        std::ifstream inputmatrix(filename);
        return Fmain(inputmatrix, q, matformat, innerdim, randomloops,
                     blocksize, maxnumcoeff,
                     initialElimination, initialSparsification);
    }
}
// ============================================================
