// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

#include "plinopt_sparsify.h"


// ============================================================
// Sparsifying and reducing coefficient diversity of a matrix
int Selector(std::istream& input,
             const FileFormat& matformat,
             const size_t maxnumcoeff) {

    Givaro::Timer chrono;
    chrono.start();

        // ============================================================
        // Read Matrix of Linear Transformation
    QRat QQ;
    QMstream ms(QQ,input);
    Matrix M(ms); M.resize(M.rowdim(),M.coldim());

    if (M.coldim() > 4) {
        std::cerr << "# \033[1;31m******   ERROR coldim > 4   ******\033[0m\n"
                  << "# \033[1;31mColumn dimension must be at most 4.\033[0m\n";
        exit(-1);
    }
#ifdef VERBATIM_PARSING
    M.write(std::clog,FileFormat::Pretty);
#endif


    Matrix TM(QQ,M.coldim(),M.rowdim()); Transpose(TM, M);

        // ============================================================
        // Initialize TICoB to identity
    Matrix TICoB(QQ,TM.rowdim(), TM.rowdim());
    for(size_t i=0; i<TICoB.coldim(); ++i) TICoB.setEntry(i,i,QQ.one);

        // ============================================================
        // Alternating sparsification and column factoring
        //    start by diagonals
    FactorDiagonals(TICoB, TM);
        //    default alternate to sparsify/factor simple things first
    SparseFactor(TICoB, TM);
        //    now try harder (with more potential combination coeffs)
    SparseFactor(TICoB, TM, maxnumcoeff, 1u, maxnumcoeff);

        // ============================================================
        // CoB = TICoB^{-T}, transposed inverse
    Matrix CoB(QQ,TICoB.coldim(), TICoB.rowdim());
    inverseTranspose(CoB, TICoB);
    chrono.stop();

        // ============================================================
        // Print resulting matrices

        // Transposed Inverse change of basis to stdlog
    TICoB.write(std::clog, matformat)<< std::endl;
    std::clog << std::string(30,'#') << std::endl;

        // change of basis to stdout
    size_t sc;
    densityProfile(std::clog << "# CoBasis profile: ", sc, CoB) << std::endl;
    std::clog << "CoB:=" << std::flush;
    CoB.write(std::cout, matformat) << std::flush;
    std::clog << ';' << std::endl;


        // sparse matrix to stdlog
    Matrix TTM(QQ,TM.coldim(), TM.rowdim()); Transpose(TTM, TM);
    TTM.write(std::clog, matformat)<< std::endl;

        // Final check that we computed a factorization M=TTM.CoB
    std::clog << std::string(30,'#') << std::endl;
    consistency(std::clog, M, TTM, CoB) << ' ' << chrono << std::endl;

    return 0;
}


// ============================================================
// Main: select between file / std::cin
//       -c # : sets the max number of coefficients per iteration
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
