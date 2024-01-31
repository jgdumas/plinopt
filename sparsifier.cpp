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

    LinBox::DenseMatrix<QRat> R(M);

    size_t sm(0);
    M.write(std::clog,FileFormat::Pretty);

        // ============================================================
        // Prints and computes density profile of M
    size_t s2;
    densityProfile(std::clog, s2, M) << std::endl;

        // ============================================================
        // Initialize CoB to identity
    Matrix CoB(QQ,M.coldim(), M.coldim());
    for(size_t i=0; i<CoB.coldim(); ++i) CoB.setEntry(i,i,QQ.one);

        // ============================================================
        // Main loop, alternating sparsification and column factoring
    size_t ss(s2);
    do {
        ss = s2;
        FactorDiagonals(CoB, M);
        Sparsifier(CoB, M, maxnumcoeff);
        densityProfile(std::clog, s2, M) << std::endl;

#ifdef DEBUG
        consistency(std::clog, M, R, CoB) << std::endl;
#endif
    } while ( s2 < ss );
    FactorDiagonals(CoB, M);

    chrono.stop();

        // ============================================================
        // Print resulting matrices
    std::clog << std::string(20,'#') << std::endl;
    std::clog << "ICoB:=";
        // change of basis to stdout
    CoB.write(std::cout, matformat);
    std::clog << ';' << std::endl;
        // sparse matrix to stdlog
    M.write(std::clog, matformat)<< std::endl;

    consistency(std::clog, M, R, CoB) << ' ' << chrono << std::endl;

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
