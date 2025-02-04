// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * Checking consistency between a program for a linear function
 *                              and the associated matrix
 * Reference:
 *   [ J-G. Dumas, B. Grenet;
 *     In-place accumulation of fast multiplication formulae
 *     ISSAC 2024, Raleigh, NC USA, pp. 16-25.
 *     (https://hal.science/hal-04167499) ]
 ****************************************************************/


#include "plinopt_programs.h"

// ============================================================
// Parsing a program and building the associated matrix
template<typename _Mat>
_Mat& PMbuilder(_Mat& A, std::istream& input) {
    std::stringstream ssin; ssin << input.rdbuf();

        // Line by line parsing
    VProgram_t ProgramVector; programParser(ProgramVector, ssin);
//     const size_t PVs { progSize(ProgramVector) };
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
        Matrix M(ms);
        const size_t m(std::max(M.rowdim(),A.rowdim()));
        const size_t n(std::max(M.coldim(),A.coldim()));
        M.resize(m,n); A.resize(m,n);
        FMatrix B(M,F);
        LinBox::DenseMatrix<Field> dB(F,m,n),dA(F,m,n),R(F,m,n);
        any2dense(dA,A); any2dense(dB,B);
        LinBox::MatrixDomain<Field> BMD(F);
        BMD.sub(R,dB,dA);

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

            return 1;
        }
    } else
        A.write(std::cout, FileFormat(5)) << std::endl;

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
                      << "        produces the matrix associated to the given program\n"
                      << "  -M f: or compares the program with the given matrix in file f\n"
                      << "  -q #: modular generation/check (default is Rationals)\n";
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
