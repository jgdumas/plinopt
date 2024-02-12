// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
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
 *     Strassen's algorithm isnot optimally accurate
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

#include "plinopt_optimize.h"


// ============================================================
// Optimizing a linear program
int Selector(std::istream& input,
             const bool printMaple, const bool printPretty,
             const bool tryDirect, const bool tryKernel) {

        // ============================================================
        // Read Matrix of Linear Transformation
    QRat QQ;
    QMstream ms(QQ,input);
    Matrix M(ms); M.resize(M.rowdim(),M.coldim());
    Matrix T(QQ);
    Transpose(T,M);


    if (printPretty) {
        M.write(std::clog,FileFormat::Pretty) << ';' << std::endl;
        std::clog << std::string(30,'#') << std::endl;
    }

    if (printMaple) {
        M.write(std::clog << "M:=",FileFormat::Maple) << ';' << std::endl;
        std::clog << std::string(30,'#') << std::endl;
    }



        // ============================================================
        // Compute naive number of operations
    int addinit(0), mulinit(0);
    for(auto iter=M.rowBegin(); iter != M.rowEnd(); ++iter)
        addinit += std::max((int)iter->size()-1,0);
    for(auto it = M.IndexedBegin(); it != M.IndexedEnd(); ++it)
        if (!isOne(abs(it.value()))) ++mulinit;

    std::clog << std::string(30,'#') << std::endl;



        // ============================================================
        // First try: optimize the whole matrix

    if (tryDirect) {
        input2Temps(M.coldim(), 'i', 't', T);

            // Cancellation-free optimization
        auto nbops( Optimizer(M, 'i', 'o', 't', 'r') );

        std::clog << std::string(30,'#') << std::endl;
        std::clog << "# " << nbops.first << "\tadditions\tinstead of " << addinit << std::endl;
        std::clog << "# " << nbops.second << "\tmultiplications\tinstead of " << mulinit << std::endl;
        std::clog << std::string(30,'#') << std::endl;

    }



        // ============================================================
        // Second try: separate independent and dependent rows

    if (tryKernel) {
        Matrix NullSpace(QQ,T.coldim(),T.coldim());
        auto Nops( nullspacedecomp(NullSpace, T) );

        if ((Nops.first !=0 || Nops.second != 0)) {
            std::clog << std::string(30,'#') << std::endl;
            std::clog << "# " << Nops.first << "\tadditions\tinstead of " << addinit << std::endl;
            std::clog << "# " << Nops.second << "\tmultiplications\tinstead of " << mulinit << std::endl;
            std::clog << std::string(30,'#') << std::endl;
        }
    }

    return 0;
}



// ============================================================
// Main: select between file / std::cin
// -D/-K options seect direct/kernel methods only (default is both)
// -P/-M option choose the printing format
int main(int argc, char** argv) {

        // ============================================================
        // Linear Transformation
    bool printMaple(false),
        printPretty(false),
        directOnly(false),
        kernelOnly(false);

    std::string filename;

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0] << " [-h|-M|-P|-K|-D] [stdin|matrixfile.sms]\n";
            exit(-1);
        } else if (args == "-M") { printMaple = true; }
        else if (args == "-P") { printPretty = true; }
        else if (args == "-D") { directOnly = true; }
        else if (args == "-K") { kernelOnly = true; }
        else { filename = args; }
    }

    const bool tryDirect(directOnly || !kernelOnly);
    const bool tryKernel(kernelOnly || !directOnly);

    QRat QQ;

    if (filename == "") {
        return Selector(std::cin, printMaple, printPretty, tryDirect, tryKernel);
    } else {
        std::ifstream inputmatrix(filename);
        if ( inputmatrix ) {
            int rt=Selector(inputmatrix, printMaple, printPretty, tryDirect, tryKernel);
            inputmatrix.close();
            return rt;
        }
    }
    return -1;
}
// ============================================================
