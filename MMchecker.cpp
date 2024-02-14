// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
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
 *                  that is modulo remainder=(srep^2-sq)
 * References:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *     Strassen's algorithm is not optimally accurate
 *     (https://hal.science/hal-04441653) ]
 *   [ J-G. Dumas, B. Grenet; Jul. 2023
 *     In-place accumulation of fast multiplication formulae
 *     (https://hal.science/hal-04167499) ]
 ****************************************************************/

#include "plinopt_library.h"

// ===============================================================
// argv[1-3]: L.sms R.sms P.sms
// argv[4]: bitsize
// argv[5-6]: modulus = argv[6]^2-argv[5], in that case 2*argv[6].bitsize is added to bitsize
int main(int argc, char ** argv) {

    if ((argc <=3) || (std::string(argv[1]) == "-h")) {
        std::clog << "Usage: " << argv[0] << " L.sms R.sms P.sms [bitsize [sq srep]]\n";
        exit(-1);
    }

        // =============================================
        // if present, result is checked modulo (srep^2-sq)
    Givaro::Integer sq(argc>6?argv[5]:""), srep(argc>6?argv[6]:"");
    Givaro::Integer modulus(srep*srep-sq);

        // =============================================
        // Reading matrices
	std::ifstream left (argv[1]), right (argv[2]), product(argv[3]);
    size_t bitsize(argc>4?atoi(argv[4]):32u);

        // To not interfere too much with probabilistic modular check ...
    if(argc>6) bitsize += srep.bitsize()<<1;


    QRat QQ;
    QMstream ls(QQ, left), rs(QQ, right), ss(QQ, product);
    Matrix L(ls), R(rs), P(ss);

	assert(L.rowdim() == R.rowdim());
	assert(L.rowdim() == P.coldim());


#ifdef VERBATIM_PARSING
    L.write(std::clog << "L:=",LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
    R.write(std::clog << "R:=",LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
    P.write(std::clog << "P:=",LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

        // =============================================
        // Random inputs
    Givaro::GivRandom generator;
    Givaro::Integer::seeding(generator.seed());

    QVector va(QQ,L.rowdim()), ua(QQ,L.coldim());
    QVector vb(QQ,R.rowdim()), ub(QQ,R.coldim());
    for(auto iter=ua.begin(); iter!=ua.end(); ++iter) QQ.random(generator, *iter, bitsize);
    for(auto iter=ub.begin(); iter!=ub.end(); ++iter) QQ.random(generator, *iter, bitsize);

        // =============================================
        // Compute matrix product via the HM algorithm
    L.apply(va,ua);
    R.apply(vb,ub);

    QVector wc(QQ,P.rowdim()), vc(QQ,P.coldim());
    for(size_t i=0; i<vc.size(); ++i) QQ.mul(vc[i],va[i],vb[i]);

    P.apply(wc,vc);

        // =============================================
        // Compute the matrix product directly
    const size_t n(std::sqrt(R.coldim()*P.rowdim()/L.coldim()));
    const size_t m(P.rowdim()/n), k(R.coldim()/n);
    LinBox::DenseMatrix<QRat> Ma(QQ,m,k), Mb(QQ,k,n), Mc(QQ,m,n);

        // row-major vectorization
    for(size_t i=0; i<ua.size(); ++i) Ma.setEntry(i/k,i%k,ua[i]);
    for(size_t i=0; i<ub.size(); ++i) Mb.setEntry(i/n,i%n,ub[i]);
    for(size_t i=0; i<wc.size(); ++i) Mc.setEntry(i/n,i%n,wc[i]);

    LinBox::MatrixDomain<QRat> BMD(QQ);
    LinBox::DenseMatrix<QRat> Rc(QQ,m,n);
    BMD.mul(Rc,Ma,Mb); // Direct matrix multiplication

    BMD.subin(Mc, Rc);

        // =============================================
        // Computations should agree modulo (srep^2-sq)
        //              as 'srep' represents sqrt('sq')
    if(argc>6) {
        for(size_t i=0; i<m; ++i) {
            for(size_t j=0; j<n; ++j) {
                Mc.refEntry(i,j) = (Mc.refEntry(i,j).nume() % modulus) / Mc.refEntry(i,j).deno();
            }
        }
    }
        // =============================================
        // Both computations should agree
    if (BMD.isZero (Mc))
        std::clog <<"# \033[1;32mSUCCESS: correct "
                  << m << 'x' << k << 'x' << n
                  << " Matrix-Multiplication!\033[0m" << std::endl;
    else{
        std::cerr << "# \033[1;31m****** ERROR, not a "
                  << m << 'x' << k << 'x' << n
                  << " MM algorithm******\033[0m"
                  << std::endl;

        ua.write(std::clog << "Ua:=", LinBox::Tag::FileFormat::Maple ) << ';' << std::endl;
        va.write(std::clog << "Va:=", LinBox::Tag::FileFormat::Maple ) << ';' << std::endl;
        ub.write(std::clog << "Ub:=", LinBox::Tag::FileFormat::Maple ) << ';' << std::endl;
        vb.write(std::clog << "Vb:=", LinBox::Tag::FileFormat::Maple ) << ';' << std::endl;
        vc.write(std::clog << "Vc:=", LinBox::Tag::FileFormat::Maple ) << ';' << std::endl;
        wc.write(std::clog << "wc:=", LinBox::Tag::FileFormat::Maple ) << ';' << std::endl;

        Ma.write(std::clog << "Ma:=", LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
        Mb.write(std::clog << "Mb:=", LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
        Mc.write(std::clog << "Mc:=", LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
        Rc.write(std::clog << "Rc:=", LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
    }

    return 0;
}
