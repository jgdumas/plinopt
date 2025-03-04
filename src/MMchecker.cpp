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
 * Usage: L.sms R.sms P.sms [bitsize [sq srep]]
 *          L.sms/R.sms/P.sms the 3 HM matrices
 *          bitsize: if present bitsize of random matrices
 *          sq srep: if present srep represents sqrt(sq)
 *                  then results is correct up to srep^2=sq,
 *                  thus checked modulo: (srep^2-sq) or /2 or /4
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

// ===============================================================
// argv[1-3]: L.sms R.sms P.sms
// argv[4]: bitsize
// argv[5-6]: modulus = (argv[6]^2-argv[5]), and /2 if argv[5]!=2
int main(int argc, char ** argv) {

    if ((argc <=3) || (std::string(argv[1]) == "-h")) {
        std::clog << "Usage:" << argv[0]
                  << "  L.sms R.sms P.sms [bitsize [sq srep]]\n"
                  << "  [bitsize]: random check with values of size 'bitsize'\n"
                  << "  [sq srep]: check is modulo (srep^2-sq) or /2 or /4\n";

        exit(-1);
    }

        // =============================================
        // if present, result is checked modulo
        // sq=2, check is mod (srep^2-sq) or /2 or /4
        // Indeed (2k+1)^2 - (2v+1) is even, so at least remove 2
        // Probability of correctness is better if modulus is prime
        // Examples:
        //    For sq=2, then srep=1013, gives mod: 1013^2-2     = 1026167
        //    For sq=3, then srep=1013, gives mod: (1013^2-3)/2 =  512083
        //    For sq=5, then srep=1013, gives mod: (1013^2-5)/4 =  256541
        //    For sq=7, then srep=1011, gives mod: (1011^2-7)/2 =  511057
    Givaro::Integer sq(argc>6?argv[5]:""), srep(argc>6?argv[6]:"");
    Givaro::Integer modulus(srep*srep-sq);
    if ( (sq%2) != 0) modulus >>=1; // (2k+1)^2 - (2v+1) is 0 mod 2
    if ( (sq%4) == 1) modulus >>=1; // (2k+1)^2 - (4v+1) is 0 mod 4
#ifdef VERBATIM_PARSING
    std::clog << std::string(30,'#') << std::endl;
    std::clog << "# Check is modulo: " << modulus << std::endl;
#endif
        // =============================================
        // Reading matrices
	std::ifstream left (argv[1]), right (argv[2]), product(argv[3]);
    size_t bitsize(argc>4?atoi(argv[4]):32u);

    QRat QQ;
    QMstream ls(QQ, left), rs(QQ, right), ss(QQ, product);
    Matrix L(ls), R(rs), P(ss);

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

    Tricounter mkn(LRP2MM(L,R,P));
    const size_t& m(std::get<0>(mkn)), k(std::get<1>(mkn)), n(std::get<2>(mkn));

        // =============================================
        // Random inputs
    QVector va(QQ,L.rowdim()), ua(QQ,L.coldim());
    QVector vb(QQ,R.rowdim()), ub(QQ,R.coldim());
    QVector wc(QQ,P.rowdim()), vc(QQ,P.coldim());

    if ( (ua.size()!=(m*k)) || (ub.size()!=(k*n)) || (wc.size()!=(m*n)) ) {
        std::cerr << "# \033[1;31m****** ERROR, outer dimension mismatch: "
                  << L.coldim() << ':' << m << 'x' << k << ' '
                  << R.coldim() << ':' << k << 'x' << n << ' '
                  << P.rowdim() << ':' << m << 'x' << n
                  << " ******\033[0m"
                  << std::endl;
        return 3;
    }

    Givaro::GivRandom generator;
    Givaro::Integer::seeding(generator.seed());
    for(auto &iter:ua) QQ.random(generator, iter, bitsize);
    for(auto &iter:ub) QQ.random(generator, iter, bitsize);

        // =============================================
        // Compute matrix product via the HM algorithm
    L.apply(va,ua);
    R.apply(vb,ub);

    for(size_t i=0; i<vc.size(); ++i) QQ.mul(vc[i],va[i],vb[i]);

    P.apply(wc,vc);

        // =============================================
        // Compute the matrix product directly
    LinBox::DenseMatrix<QRat> Ma(QQ,m,k), Mb(QQ,k,n), Delta(QQ,m,n);

        // row-major vectorization
    for(size_t i=0; i<ua.size(); ++i) Ma.setEntry(i/k,i%k,ua[i]);
    for(size_t i=0; i<ub.size(); ++i) Mb.setEntry(i/n,i%n,ub[i]);
    for(size_t i=0; i<wc.size(); ++i) Delta.setEntry(i/n,i%n,wc[i]);

    LinBox::MatrixDomain<QRat> BMD(QQ);
    LinBox::DenseMatrix<QRat> Mc(QQ,m,n);
    BMD.mul(Mc,Ma,Mb); // Direct matrix multiplication

    BMD.subin(Delta, Mc);

        // =============================================
        // Computations should agree modulo (srep^2-sq)
        //              as 'srep' represents sqrt('sq')
    if(argc>6) {
        for(size_t i=0; i<m; ++i) {
            for(size_t j=0; j<n; ++j) {
                Delta.setEntry(i,j, Givaro::Rational(
                    Delta.getEntry(i,j).nume() % modulus,
                    Delta.getEntry(i,j).deno() ) );
            }
        }
    }

        // =============================================
        // Both computations should agree
    if (BMD.isZero (Delta))
        std::clog <<"# \033[1;32mSUCCESS: correct "
                  << m << 'x' << k << 'x' << n
                  << " Matrix-Multiplication!\033[0m" << std::endl;
    else{
        std::cerr << "# \033[1;31m****** ERROR, not a "
                  << m << 'x' << k << 'x' << n
                  << " MM algorithm******\033[0m"
                  << std::endl;

        ua.write(std::clog << "Ua:=", FileFormat::Maple ) << ';' << std::endl;
        va.write(std::clog << "Va:=", FileFormat::Maple ) << ';' << std::endl;
        ub.write(std::clog << "Ub:=", FileFormat::Maple ) << ';' << std::endl;
        vb.write(std::clog << "Vb:=", FileFormat::Maple ) << ';' << std::endl;
        vc.write(std::clog << "Vc:=", FileFormat::Maple ) << ';' << std::endl;
        wc.write(std::clog << "wc:=", FileFormat::Maple ) << ';' << std::endl;

        Ma.write(std::clog << "Ma:=", FileFormat::Maple) << ';' << std::endl;
        Mb.write(std::clog << "Mb:=", FileFormat::Maple) << ';' << std::endl;
            // Correct value
        Mc.write(std::clog << "Mc:=", FileFormat::Maple) << ';' << std::endl;
            // Difference with computed value
        Delta.write(std::clog << "Df:=", FileFormat::Maple) << ';' << std::endl;
        BMD.addin(Delta, Mc);
            // Computed value
        Delta.write(std::clog << "Rc:=", FileFormat::Maple) << ';' << std::endl;

        return 1;
    }

    return 0;
}
