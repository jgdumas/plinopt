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

#include <givaro/modular.h>
#include <givaro/givpoly1.h>
#include "plinopt_library.h"
#include "plinopt_polynomial.h"


// ===============================================================
void usage(const char* prgname) {
    std::clog << "Usage:" << prgname
	      << "  [-h|-b #|-m/-q #] L.sms R.sms P.sms\n"
	      << "  [-b b]: random check with values of size 'bitsize'\n"
	      << "  [-m/-q m]: check is modulo (mod) or /2^k (default no)\n"
	      << "  [-I x]: indeterminate (default is X)\n"
	      << "  [-P P(x)]: irreducible polynomial (default is none)\n";

    exit(-1);
}


// ===============================================================
// Generic random element with bitsize
template<typename Domain>
typename Domain::Element& RandomElt(typename Domain::Element& e,
				    const Domain& D,
				    Givaro::GivRandom& generator,
				    const size_t bitsize) {
    return D.random(generator, e, bitsize);
}

// Specialization for Modular, as bitsize is meanigless for its random
template<>
typename Givaro::Integer& RandomElt(Givaro::Integer& e,
				    const Givaro::Modular<Givaro::Integer>& D,
				    Givaro::GivRandom& generator,
				    const size_t bitsize) {
    return D.random(generator, e);
}
// ===============================================================

using QPol = PLinOpt::PRing<PLinOpt::QRat>;

// ===============================================================
// Generic verification, within Field, of polynomials over Base
template<typename Base, typename Field>
int PMchecker(const Base& BB, const Field& FF, const size_t bitsize, const QPol& QQX,
	      const QPol::Element& QP, const std::vector<std::string>& filenames) {
    typedef LinBox::MatrixStream<Base> BMstream;
    typedef LinBox::SparseMatrix<Base,
	LinBox::SparseMatrixFormat::SparseSeq > BMatrix;
    using FMatrix=typename BMatrix::template rebind<Field>::other;
    typedef LinBox::DenseVector<Field> FVector;
    using FPol = PLinOpt::PRing<Field>;

	// =============================================
	// Reading matrices
	std::ifstream left (filenames[0]), right (filenames[1]), post(filenames[2]);

    using PLinOpt::FileFormat;

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

    PLinOpt::Tricounter mkn(PLinOpt::LRP2MM(L,R,P));
    const size_t& m(std::get<0>(mkn)), k(std::get<1>(mkn)), n(std::get<2>(mkn));

	// =============================================
	// Random inputs
    FVector va(FF,L.rowdim()), ua(FF,L.coldim());
    FVector vb(FF,R.rowdim()), ub(FF,R.coldim());
    FVector wc(FF,P.rowdim()), vc(FF,P.coldim());

    Givaro::GivRandom generator;
    Givaro::Integer::seeding(generator.seed());
    for(auto &iter:ua) RandomElt(iter, FF, generator, bitsize);
    for(auto &iter:ub) RandomElt(iter, FF, generator, bitsize);

	// =============================================
	// Compute matrix product via the HM algorithm
    L.apply(va,ua);
    R.apply(vb,ub);

    for(size_t i=0; i<vc.size(); ++i) FF.mul(vc[i],va[i],vb[i]);

    P.apply(wc,vc);

	// =============================================
	// Compute the polynomial product directly

    FPol FX(FF, QQX.getIndeter());
    typename FPol::Element Pa, Pb, Pc, Delta;
    FX.init(Pa, Givaro::Degree(L.coldim()-1));
    FX.init(Pb, Givaro::Degree(R.coldim()-1));
    FX.init(Delta, Givaro::Degree(P.rowdim()-1));


	// row-major vectorization
    for(size_t i=0; i<ua.size(); ++i) Pa[i]=ua[i];
    for(size_t i=0; i<ub.size(); ++i) Pb[i]=ub[i];
    for(size_t i=0; i<wc.size(); ++i) Delta[i]=wc[i];

	// Direct multiplication
    FX.mul(Pc,Pa,Pb);

	// Diffrence with the value returned by the program
    FX.subin(Delta, Pc);

    typename FPol::Element irred;
    if (QP.size() >0) {
	irred.resize(QP.size());
	for(size_t i=0; i<QP.size(); ++i) {
	    FF.init(irred[i]);
	    typename Field::Element tmp; FF.init(tmp);
	    FF.init(tmp, QP[i].deno());
	    FF.init(irred[i],QP[i].nume());
	    FF.divin(irred[i],tmp);
	}

	    // Computation modulo the irreducible
	FX.modin(Delta, irred);
    }

	// =============================================
	// Both computations should agree
    if (FX.isZero (Delta)) {
	FF.write(std::clog <<"# \033[1;32mSUCCESS: correct "
		 << (L.coldim()-1) << 'o' << (R.coldim()-1) << 'o' << (P.rowdim()-1)
		 << " Polynomial-Multiplication \033[0m");
	if (irred.size()>0) FX.write(std::clog << " mod ", irred);
	std::clog << std::endl;
    } else {
	FF.write(std::cerr << "# \033[1;31m****** ERROR, not a "
		 << (L.coldim()-1) << 'o' << (R.coldim()-1) << 'o' << (P.rowdim()-1)
		 << " MM algorithm******\033[0m");
	if (irred.size()>0) FX.write(std::clog << " mod ", irred);
	std::clog << std::endl;

	ua.write(std::clog << "Ua:=", FileFormat::Maple ) << ';' << std::endl;
	va.write(std::clog << "Va:=", FileFormat::Maple ) << ';' << std::endl;
	ub.write(std::clog << "Ub:=", FileFormat::Maple ) << ';' << std::endl;
	vb.write(std::clog << "Vb:=", FileFormat::Maple ) << ';' << std::endl;
	vc.write(std::clog << "Vc:=", FileFormat::Maple ) << ';' << std::endl;
	wc.write(std::clog << "wc:=", FileFormat::Maple ) << ';' << std::endl;

	FX.write(std::clog << "Pa:=", Pa) << ';' << std::endl;
	FX.write(std::clog << "Pb:=", Pb) << ';' << std::endl;
	    // Correct value
	FX.write(std::clog << "Pc:=", Pc) << ';' << std::endl;
	    // Difference with computed value
	FX.write(std::clog << "Df:=", Delta) << ';' << std::endl;
	FX.addin(Delta, Pc);
	    // Computed value
	FX.write(std::clog << "Rc:=", Delta) << ';' << std::endl;

	return 1;
    }
    return 0;
}

// ===============================================================
int main(int argc, char ** argv) {

    size_t bitsize(32u);
    Givaro::Integer modulus(0u);
    PLinOpt::QRat QQ;
    QPol QQX(QQ, 'X');
    QPol::Element QP;
    std::vector<std::string> filenames;

	// Parse arguments

    for(int i=1; i<argc; ++i) {
	std::string args(argv[i]);
	if (args[0] == '-') {
	    if (args[1] == 'h') { usage(argv[0]); }
	    else if (args[1] == 'b') { bitsize = atoi(argv[++i]); }
	    else if ((args[1] == 'm') || (args[1] == 'q')) { modulus = Givaro::Integer(argv[++i]); }
	    else if (args[1] == 'I') {
		QQX.setIndeter(Givaro::Indeter(argv[++i][0]));
	    }
	    else if (args[1] == 'P') {
		std::stringstream sargs( argv[++i] );
		QQX.read(sargs, QP);
	    }
	} else { filenames.push_back(args); }
    }
    if (filenames.size() < 3) { usage(argv[0]); }


	// Select verification

    if (modulus>0) {
	    // Compute modular verification of rational polynomials
	Givaro::Modular<Givaro::Integer> F(modulus);
	return PMchecker(QQ, F, bitsize, QQX, QP, filenames);
    } else
	    // Compute rational verification of rational polynomials
	return PMchecker(QQ, QQ, bitsize, QQX, QP, filenames);
}
