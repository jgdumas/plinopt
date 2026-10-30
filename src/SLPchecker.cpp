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

namespace PLinOpt {

// ============================================================
// Checking a linear program with a matrix
template<typename Field>
int SLPcheck(const std::string& prgname, const std::string& matname,
	    const Field& F) {

	// ============================================================
	// Rebind matrix type over sub field matrix type
    using FMatrix=typename Matrix::template rebind<Field>::other;

    FMatrix A(F);
    Pair<size_t> Ops;

    if (prgname == "") {
      Ops = SLPbuilder(A, std::cin, true);
    } else {
        std::ifstream ifile(prgname);
        if ( ifile ) {
	  Ops = SLPbuilder(A, ifile, true);
            ifile.close();
        }
    }

    if (matname != "") {
      Pair<size_t> nbelts;
      std::ifstream matfile(matname);
      if (! areEqual(A, matfile, nbelts, true)) return 1;
      std::clog <<"# \033[1;32mSUCCESS: correct SLP for "
		<< A.rowdim() << 'x' << A.coldim()
		<< " (" << nbelts.first << '+'
		<< '|' << nbelts.second << 'x'
		<< ") : " << Ops
		<< " Matrix-Vector multiplication!\033[0m" << std::endl;
    } else {
      A.write(std::cout, FileFormat(5));
    }
    return 0;
}

} // End of namespace PLinOpt
// ============================================




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

    if (Givaro::isZero(q)) {
        PLinOpt::QRat QQ;
        return PLinOpt::SLPcheck(prgname, matname, QQ);
    } else {
        Givaro::Modular<Givaro::Integer> FF(q);
        return PLinOpt::SLPcheck(prgname, matname, FF);
    }
}
// ============================================================
