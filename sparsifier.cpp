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
    size_t sc,sb,sr;
    densityProfile(std::clog << "# Initial profile: ", sc, M)
                             << std::endl;

#ifdef VERBATIM_PARSING
    M.write(std::clog,FileFormat::Pretty);
#endif

    Givaro::Timer elapsed;
    Matrix CoB(QQ,M.coldim(), M.coldim());
    Matrix Res(QQ,M.rowdim(), M.coldim());
    blockSparsifier(elapsed, CoB, Res, M, blocksize,
                    QQ, matformat, maxnumcoeff, initialElimination);

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
        << " \033[1;36m"
        << sr << " non-zeroes (" << sb << " alt.) instead of " << sc
        << "\033[0m:" << ' ' << elapsed << std::endl;

    return 0;
}


// ============================================================
// Main: select between file / std::cin
//       -c #: sets the max number of coefficients per iteration
//       -M/-P/-S: selects the ouput format
//       -b #: states the blocking dimension
//       -U [1|0]: uses an initial LU factorization | or not
int main(int argc, char** argv) {

    FileFormat matformat = FileFormat::Pretty;
    std::string filename;
    size_t maxnumcoeff(COEFFICIENT_SEARCH); // default max coefficients
    size_t blocksize(4u);                   // default column block size
    bool initialElimination(true);

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0] << " [-h|-M|-P|-S|-c #|-U [1|0]] [stdin|matrixfile.sms]\n";

            std::clog << "  -c #: max number of coefficients per iteration\n"
                      << "  -b #: states the blocking dimension\n"
                      << "  -U [1|0]: initial LU factorization | or not\n"
                      << "  -M/-P/-S: selects the ouput format\n";

            exit(-1);
        }
        else if (args == "-M") { matformat = FileFormat(1); } // Maple
        else if (args == "-S") { matformat = FileFormat(5); } // SMS
        else if (args == "-P") { matformat = FileFormat(8); } // Pretty
        else if (args == "-c") { maxnumcoeff = atoi(argv[++i]); }
        else if (args == "-b") { blocksize = atoi(argv[++i]); }
        else if (args == "-U") { initialElimination = atoi(argv[++i]); }
        else { filename = args; }
    }

    if (filename == "") {
        return Selector(std::cin, matformat,
                        blocksize, maxnumcoeff, initialElimination);
    } else {
        std::ifstream inputmatrix(filename);
        return Selector(inputmatrix, matformat,
                        blocksize, maxnumcoeff, initialElimination);
    }

    return -1;
}
// ============================================================
