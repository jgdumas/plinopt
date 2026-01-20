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
#include <givaro/givquotientdomain.h>
#include "plinopt_library.h"
#include "plinopt_polynomial.h"


// ===============================================================
void usage(const char* prgname) {
    std::clog << "Usage:" << prgname
              << "  [-h|-b #|-m #|-r # # #|-I x|-P P(x)] L.sms R.sms P.sms\n"
              << "  [-b b]: random check with values of size 'bitsize'\n"
              << "  [-m m]: check is modulo (mod) or /2^k (default no)\n"
              << "  [-r r e s]: check is modulo (r^e-s) or /2^k (default no)\n"
              << "  [-I x]: indeterminate (default is X)\n"
              << "  [-P P(x)]: modular polynomial(default is none)\n";

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


// ===============================================================
// Missing LinBox specialization of quotient homomorphism
template <class Field>
class LinBox::Hom<Field, Givaro::QuotientDom<Field> > {
public:
    typedef Field Source;
    typedef typename Givaro::QuotientDom<Field> Target;
    typedef typename Source::Element SrcElt;
    typedef typename Target::Element Elt;

    Hom(const Source& S, const Target& T) : _source(S), _target(T) {}

    Elt& image(Elt& t, const SrcElt& s) {
        _target.assign (t, s); // this will compute the modular image
        return t;
    }

    SrcElt& preimage(SrcElt& s, const Elt& t) {
        _source.assign (s, t); // this will just copy
        return s;
    }

    const Source& source() { return _source;}
    const Target& target() { return _target;}

private:
    const Source& _source;
    const Target& _target;
}; // end Hom
// ===============================================================


// ===============================================================
// Generic verification, within Field, of polynomials over Base
template<typename Base, typename Field>
int PMchecker(const Base& BB, const Field& FF, const size_t bitsize,
              const std::vector<std::string>& filenames) {
    typedef LinBox::MatrixStream<Base> BMstream;
    typedef LinBox::SparseMatrix<Base,
        LinBox::SparseMatrixFormat::SparseSeq > BMatrix;
    using FMatrix=typename BMatrix::template rebind<Field>::other;
    typedef LinBox::DenseVector<Field> FVector;

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

//     if ( (ua.size()!=(m*k)) || (ub.size()!=(k*n)) || (wc.size()!=(m*n)) ) {
//         std::cerr << "# \033[1;31m****** ERROR, outer dimension mismatch: "
//                   << L.coldim() << ':' << m << 'x' << k << ' '
//                   << R.coldim() << ':' << k << 'x' << n << ' '
//                   << P.rowdim() << ':' << m << 'x' << n
//                   << " ******\033[0m"
//                   << std::endl;
//         return 3;
//     }

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

    Givaro::Poly1Dom<Field, Givaro::Dense> PX(FF,Givaro::Indeter("X"));
    typename Givaro::Poly1Dom<Field, Givaro::Dense>::Element Pa, Pb, Delta;
    PX.init(Pa, Givaro::Degree(L.coldim()-1));
    PX.init(Pb, Givaro::Degree(R.coldim()-1));
    PX.init(Delta, Givaro::Degree(P.rowdim()-1));
    

        // row-major vectorization
    for(size_t i=0; i<ua.size(); ++i) Pa[i]=ua[i];
    for(size_t i=0; i<ub.size(); ++i) Pb[i]=ub[i];
    for(size_t i=0; i<wc.size(); ++i) Delta[i]=wc[i];

    typename Givaro::Poly1Dom<Field, Givaro::Dense>::Element Pc;
    PX.mul(Pc,Pa,Pb);
    PX.subin(Delta, Pc);

        // =============================================
        // Both computations should agree
    if (PX.isZero (Delta))
        FF.write(std::clog <<"# \033[1;32mSUCCESS: correct "
                 << L.coldim() << 'o' << R.coldim() << 'o' << P.rowdim()
                 << " Polynomial-Multiplication \033[0m") << std::endl;
    else{
        FF.write(std::cerr << "# \033[1;31m****** ERROR, not a "
                 << L.coldim() << 'o' << R.coldim() << 'o' << P.rowdim()
                 << " MM algorithm******\033[0m") << std::endl;

        ua.write(std::clog << "Ua:=", FileFormat::Maple ) << ';' << std::endl;
        va.write(std::clog << "Va:=", FileFormat::Maple ) << ';' << std::endl;
        ub.write(std::clog << "Ub:=", FileFormat::Maple ) << ';' << std::endl;
        vb.write(std::clog << "Vb:=", FileFormat::Maple ) << ';' << std::endl;
        vc.write(std::clog << "Vc:=", FileFormat::Maple ) << ';' << std::endl;
        wc.write(std::clog << "wc:=", FileFormat::Maple ) << ';' << std::endl;

        PX.write(std::clog << "Pa:=", Pa) << ';' << std::endl;
        PX.write(std::clog << "Pb:=", Pb) << ';' << std::endl;
            // Correct value
        PX.write(std::clog << "Pc:=", Pc) << ';' << std::endl;
            // Difference with computed value
        PX.write(std::clog << "Df:=", Delta) << ';' << std::endl;
        PX.addin(Delta, Pc);
            // Computed value
        PX.write(std::clog << "Rc:=", Delta) << ';' << std::endl;

        return 1;
    }
    return 0;
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
            else if (args[1] == 'm') { modulus = Givaro::Integer(argv[++i]); }
            else if (args[1] == 'r') {
                Givaro::Integer r(argv[++i]);
                size_t e(atoi(argv[++i]));
                Givaro::Integer s(argv[++i]);
                modulus = pow(r,e)-s;
            }
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
            // Remove powers of 2 for better probabilities
        while( (modulus % 2) == 0 ) { modulus >>=1; };
        if (modulus == 1) modulus = 2;
            // Compute modular verification of rational matrices
        Givaro::Modular<Givaro::Integer> F(modulus);
        return PMchecker(QQ, F, bitsize, filenames);
    }
    else if (QP.size()>0) {
            // Compute modular verification of polynomial matrices
        QQuo QQXm(QQX, QP);
        return PMchecker(QQX, QQXm, bitsize, filenames);
    } else
            // Compute rational verification of rational matrices
        return PMchecker(QQ, QQ, bitsize, filenames);
}
