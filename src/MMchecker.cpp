// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * Randomly checks a matrix multiplication algorithm
 *        given as a row-major HM representation:
 *        [ a11 a12 ]
 *        [ a21 a22 ] is vectorized as [a11 a12 a21 a22]
 *
 * Examples:
 * > ./bin/MMchecker data/2x2x2_7_Strassen_{L,R,P}.sms
 * > ./bin/MMchecker data/2x2x2_7_DPS-accurate_{L,R,P}.sms -m 513083
 * > ./bin/MMchecker data/2x2x2_7_DPS-accurate_{L,R,P}.sms -r 1013 2 3
 * > ./bin/MMchecker data/2x2x2_7_DPS-accurate-X_{L,R,P}.sms -P "X^2-3"
 *
 * References:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
 *   [ J-G. Dumas, B. Grenet;
 *     In-place accumulation of fast multiplication formulae
 *     ISSAC 2024, Raleigh, NC USA, pp. 16-25.
 *     (https://hal.science/hal-04167499) ]
 ****************************************************************/

#include "plinopt_library.h"
#include "plinopt_polynomial.h"

// ===============================================================
void usage(const char* prgname) {
    std::clog << "Usage:" << prgname
              << "  [-h|-b #|-m/-q #|-r # # #|-I x|-P P(x)] L.sms R.sms P.sms\n"
              << "  [-b b]: random check with values of size 'bitsize'\n"
              << "  [-m/-q m]: check is modulo (mod) or (mod/2^k) (default no)\n"
              << "  [-r r e s]: check is modulo (r^e-s) or ((r^e-s)/2^k) (default no)\n"
              << "  [-I x]: indeterminate (default is X)\n"
              << "  [-P P(x)]: modular polynomial (default is none)\n";

    exit(-1);
}


// ===============================================================
// Reading matrices over Base
template<typename Base, typename Field>
int fMMchecker(const Base& BB, const Field& FF, const size_t bitsize,
              const std::vector<std::string>& filenames) {
    typedef LinBox::MatrixStream<Base> BMstream;
    typedef LinBox::SparseMatrix<Base,
        LinBox::SparseMatrixFormat::SparseSeq > BMatrix;
    using FMatrix=typename BMatrix::template rebind<Field>::other;
    using PLinOpt::FileFormat;

        // =============================================
        // Reading matrices
	std::ifstream left (filenames[0]), right (filenames[1]), post(filenames[2]);

    BMstream ls(BB, left), rs(BB, right), ss(BB, post);
    BMatrix BL(ls), BR(rs), BP(ss);
    FMatrix L(BL,FF), R(BR,FF), P(BP,FF);

    if ( (L.rowdim() != R.rowdim()) || (L.rowdim() != P.coldim()) ) {
        std::cerr << "# \033[1;31m****** ERROR, inner dimension mismatch: "
                  << L.rowdim() << "(.)" << R.rowdim() << '|' << P.coldim()
                  << " ******\033[0m"
                  << std::endl;
        return 2;
    }

#ifdef VERBATIM_PARSING
    L.write(std::clog << "L:=",FileFormat::Maple) << ';' << std::endl;
    R.write(std::clog << "R:=",FileFormat::Maple) << ';' << std::endl;
    P.write(std::clog << "P:=",FileFormat::Maple) << ';' << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

    return PLinOpt::MMchecker(FF, bitsize, L, R, P);
}


// ===============================================================
int main(int argc, char ** argv) {

    size_t bitsize(32u);
    Givaro::Integer modulus(0u);
    PLinOpt::QRat QQ;
    using QPol = PLinOpt::PRing<PLinOpt::QRat>;
    using QQuo = Givaro::QuotientDom<QPol>;
    QPol QQX(QQ, 'X');
    QPol::Element QP;
    std::vector<std::string> filenames;

        // Parse arguments

    for(int i=1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args[0] == '-') {
            if (args[1] == 'h') { usage(argv[0]); }
            else if (args[1] == 'b') { bitsize = atoi(argv[++i]); }
            else if ((args[1] == 'm') || (args[1] == 'q')) {
                modulus = Givaro::Integer(argv[++i]);
            } else if (args[1] == 'r') {
                Givaro::Integer r(argv[++i]);
                size_t e(atoi(argv[++i]));
                Givaro::Integer s(argv[++i]);
                modulus = pow(r,e)-s;
            } else if (args[1] == 'I') {
                QQX.setIndeter(Givaro::Indeter(argv[++i][0]));
            } else if (args[1] == 'P') {
                std::stringstream sargs( argv[++i] );
                QQX.read(sargs, QP);
            }
        } else { filenames.push_back(args); }
    }
    if (filenames.size() < 3) { usage(argv[0]); }


        // Select verification

    if (modulus>0) {
            // Remove powers of 2 for better probabilities
        while( (modulus % 2) == 0 ) { modulus >>=1; };
        if (modulus == 1) modulus = 2;
            // Compute modular verification of rational matrices
        Givaro::Modular<Givaro::Integer> F(modulus);
        return fMMchecker(QQ, F, bitsize, filenames);
    }
    else if (QP.size()>0) {
            // Compute modular verification of polynomial matrices
        QQuo QQXm(QQX, QP);
        return fMMchecker(QQX, QQXm, bitsize, filenames);
    } else
            // Compute rational verification of rational matrices
        return fMMchecker(QQ, QQ, bitsize, filenames);
}
