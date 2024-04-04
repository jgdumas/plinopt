// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * Transposition of programs via Tellegen's principle
 * Program syntax: see main Tellegen function below
 * References:
 *   [ Tellegen's principle into practice.
 *     ISSAC'03:37-44, A. Bostan, G. Lecerf, É. Schost
 *     https://doi.org/10.1145/860854.860870 ]
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *     Strassen's algorithm is not optimally accurate
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/


#include "plinopt_programs.h"

// ============================================================
// Parsing a program and building the associated matrix
template<typename _Mat>
_Mat& PMbuilder(_Mat& A, std::istream& input) {
    std::stringstream ssin; ssin << input.rdbuf();

        // Line by line parsing
    VProgram_t ProgramVector; programParser(ProgramVector, ssin);
    const size_t PVs { progSize(ProgramVector) };
    std::clog << std::string(40,'#') << std::endl;
    parenthesisExpand(ProgramVector);
    std::clog << std::string(40,'#') << std::endl;

    return matrixBuilder(A, ProgramVector);
}
// ============================================================



// ============================================================
// Checking a linear program with a matrix
template<typename Field>
int PMcheck(const std::string& prgname, const std::string& matname,
            const Field& F) {

        // ============================================================
        // Rebind matrix type over sub field matrix type
    using FMatrix=typename Matrix::template rebind<Field>::other;

    FMatrix A(F);

    if (prgname == "") {
        PMbuilder(A, std::cin);
    } else {
        std::ifstream ifile(prgname);
        if ( ifile ) {
            PMbuilder(A, ifile);
            ifile.close();
        }
    }

    if (matname != "") {
        std::ifstream matfile(matname);
        QRat QQ;
        QMstream ms(QQ, matfile);
        Matrix M(ms); M.resize(M.rowdim(),M.coldim());
        FMatrix B(M, F);
        const size_t m(B.rowdim()), n(B.coldim());
        FMatrix R(F,m,n);
        LinBox::MatrixDomain<Field> BMD(F);
        BMD.sub(R,B,A);

        if (BMD.isZero (R))
            std::clog <<"# \033[1;32mSUCCESS: correct "
                      << m << 'x' << n
                      << " Matrix-Vector multiplication!\033[0m" << std::endl;
        else {
            std::cerr << "# \033[1;31m****** ERROR, not a "
                      << m << 'x' << n
                      << " m-v algorithm******\033[0m"
                      << std::endl;

            A.write(std::clog<<"# Program:\n", FileFormat::Pretty);
            B.write(std::clog<<"\n# Matrix:\n", FileFormat::Pretty);
            R.write(std::clog<<"\n# diff:\n", FileFormat::Pretty) << std::endl;
        }
    } else
        A.write(std::cout, FileFormat::SMS) << std::endl;

    return 0;
}




// ============================================================
// Main: select between file / std::cin
int main(int argc, char** argv) {
    std::string prgname, matname;
    Givaro::Integer q(0u);

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0] << "[-q #] [-M file.sms] [stdin|file.prg] \n"
                      << "  -M f: Matrix file (compared to program file)\n"
                      << "  -q #: search modulo (default is Rationals)\n";
            exit(-1);
        }
        else if (args == "-M") { matname = std::string(argv[++i]); }
        else if (args == "-q") { q = Givaro::Integer(argv[++i]); }
        else { prgname = args; }
    }

    if (! Givaro::isZero(q)) {
        Givaro::Modular<Givaro::Integer> FF(q);
        return PMcheck(prgname, matname, FF);
    } else {
        QRat QQ;
        return PMcheck(prgname, matname, QQ);
    }
}
// ============================================================
