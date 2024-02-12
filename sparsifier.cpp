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
int Selector(std::istream& input,
             const FileFormat& matformat,
             const size_t maxnumcoeff) {

        // ============================================================
        // Read Matrix of Linear Transformation
    QRat QQ;
    QMstream ms(QQ,input);
    Matrix M(ms); M.resize(M.rowdim(),M.coldim());

#ifdef VERBATIM_PARSING
    M.write(std::clog,FileFormat::Pretty);
#endif

    Givaro::Timer elapsed,chrono;

    Matrix CoB(QQ,M.coldim(), M.coldim());
    Matrix Res(QQ,M.rowdim(), M.coldim());

//         // ============================================================
//         // Sparsify matrix as a whole
//     sparseAlternate(chrono, CoB, Res, M, matformat, maxnumcoeff);

        // ============================================================
        // Deal with blocks of columns
    std::vector<Matrix> vC, vR, vM;
    separateColumnBlocks(vM, M, 4);
    for(const auto& mat: vM) {
        vC.emplace_back(QQ,M.coldim(), M.coldim());
        vR.emplace_back(QQ,M.rowdim(), M.coldim());

        sparseAlternate(chrono, vC.back(), vR.back(), mat,
                        matformat, maxnumcoeff);

        elapsed += chrono;
#ifdef DEBUG
        std::clog << std::string(30,'#') << std::endl;
        consistency(std::clog, mat, vR.back(), vC.back()) << ' ' << chrono << std::endl;
#endif
    }

        // Build resulting matrices
    diagonalMatrix(CoB, vC);
    augmentedMatrix(Res, vR);


        // ============================================================
        // Print resulting matrices

        // change of basis to stdout
    size_t sc;
    densityProfile(std::clog << "# Alternate basis profile: ", sc, CoB)
                             << std::endl;
    CoB.write(std::cout, matformat) << std::endl;


        // residuum sparse matrix to stdlog
    densityProfile(std::clog << "# Sparse residuum profile: ", sc, CoB)
                             << std::endl;
    Res.write(std::clog, matformat)<< std::endl;


        // Final check that we computed a factorization M=Res.CoB
    std::clog << std::string(30,'#') << std::endl;
    consistency(std::clog, M, Res, CoB) << ' ' << elapsed << std::endl;

    return 0;
}


// ============================================================
// Main: select between file / std::cin
//       -c # : sets the max number of coefficients per iteration
//       -f
int main(int argc, char** argv) {

    FileFormat matformat = FileFormat::Pretty;
    std::string filename;
    size_t maxnumcoeff(COEFFICIENT_SEARCH); // default max coefficients

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0] << " [-h|-M|-P|-S|-c #] [stdin|matrixfile.sms]\n";
            exit(-1);
        }
        else if (args == "-M") { matformat = FileFormat(1); }
        else if (args == "-P") { matformat = FileFormat(8); }
        else if (args == "-S") { matformat = FileFormat(5); }
        else if (args == "-c") { maxnumcoeff = atoi(argv[++i]); }
        else { filename = args; }
    }

    if (filename == "") {
        return Selector(std::cin, matformat, maxnumcoeff);
    } else {
        std::ifstream inputmatrix(filename);
        return Selector(inputmatrix, matformat, maxnumcoeff);
    }

    return -1;
}
// ============================================================
