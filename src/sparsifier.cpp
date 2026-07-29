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

namespace PLinOpt {
// ============================================================
// Sparsifying and reducing coefficient diversity of a matrix
int Selector(std::istream& input, const FileFormat& matformat,
             const size_t blocksize, const size_t maxnumcoeff,
             const bool initialElimination) {

        // ============================================================
        // Read Matrix of Linear Transformation
    QRat QQ;
    QMstream ms(QQ,input);
    Matrix M(ms); M.resize(M.rowdim(),M.coldim());
    size_t sc;
    densityProfile(std::clog << "# [SPRF] Initial profile: ", sc, M)
                             << std::endl;

#ifdef VERBATIM_PARSING
    M.write(std::clog,FileFormat::Pretty);
#endif

    Givaro::Timer elapsed;
    Matrix CoB(QQ,M.coldim(), M.coldim());
    Matrix Res(QQ,M.rowdim(), M.coldim());
    blockSparsifier(elapsed, CoB, Res, M, blocksize,
                    matformat, maxnumcoeff, initialElimination);

        // ============================================================
        // Print resulting matrices
    const std::string fctz{"[SPRF]"};
    profileConsistency(matformat, elapsed, M, sc, fctz, Res, -1, CoB, 1);

    return 0;
}


} // End of namespace PLinOpt
// ============================================


// ============================================================
// Main: select between file / std::cin
//       -c #: sets the max number of coefficients per iteration
//       -M/-P/-S: selects the ouput format
//       -b #: states the blocking dimension
//       -U [1|0]: uses an initial LU factorization | or not
int main(int argc, char** argv) {
    using PLinOpt::FileFormat;

    FileFormat matformat = FileFormat::Pretty;
    std::string filename;
    size_t maxnumcoeff(COEFFICIENT_SEARCH); // default max coefficients
    size_t blocksize(4u);                   // default column block size
    bool initialElimination(true);

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog
                << "Usage: " << argv[0]
                << " [-h|-M|-P|-S|-L|-c #|-U [1|0]] [stdin|matfile.sms]\n"
                << "  -c #: max number of coefficients per iteration (default "
                << maxnumcoeff << ")\n"
                << "  -b #: states the blocking dimension (default "
                << blocksize << ")\n"
                << "  -U [1|0]: initial LU factorization (default) or not\n"
                << "  -M/-P/-S/-L: selects the ouput format\n";

            exit(-1);
        }
        else if (args == "-M") { matformat = FileFormat(1); } // Maple
        else if (args == "-S") { matformat = FileFormat(5); } // SMS
        else if (args == "-P") { matformat = FileFormat(8); } // Pretty
        else if (args == "-L") { matformat = FileFormat(12); }// Linalg
        else if (args == "-c") { maxnumcoeff = atoi(argv[++i]); }
        else if (args == "-b") { blocksize = atoi(argv[++i]); }
        else if (args == "-U") { initialElimination = atoi(argv[++i]); }
        else { filename = args; }
    }

    if (filename == "") {
        return PLinOpt::Selector(std::cin, matformat,
                                 blocksize, maxnumcoeff, initialElimination);
    } else {
        std::ifstream inputmatrix(filename);
        return PLinOpt::Selector(inputmatrix, matformat,
                                 blocksize, maxnumcoeff, initialElimination);
    }

    return -1;
}
// ============================================================
