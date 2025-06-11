// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * Randomly checks a matrix multiplication algorithm modulo a polynomial
 *
 * Usage: L.sms R.sms P.sms [bitsize]
 *          L.sms/R.sms/P.sms the 3 HM matrices
 *          bitsize: if present bitsize of random matrices
 *        The modular polynomial is then given in std::cin
 ****************************************************************/

#include <iostream>

#include <givaro/givrational.h>
#include <givaro/givquotientdomain.h>
#include <linbox/ring/polynomial-ring.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/util/matrix-stream.h>

using LinBox::Tag::FileFormat;

typedef Givaro::QField<Givaro::Rational> QRat;

typedef LinBox::PolynomialRing<QRat> QPol;
typedef Givaro::QuotientDom<QPol> PMods;

typedef LinBox::MatrixStream<PMods> QMstream;
typedef LinBox::SparseMatrix<PMods,
                             LinBox::SparseMatrixFormat::SparseSeq > Matrix;
typedef LinBox::DenseVector<PMods> QVector;


// ===============================================================
// argv[1-3]: L.sms R.sms P.sms
// argv[4]: bitsize
int main(int argc, char ** argv) {

    if ((argc <=3) || (std::string(argv[1]) == "-h")) {
        std::clog << "Usage: " << argv[0] << " L.sms R.sms P.sms [bitsize]\n";
        exit(-1);
    }

        // =============================================
        // Reading matrices
	std::ifstream left (argv[1]), right (argv[2]), product(argv[3]);
    size_t bitsize(argc>4?atoi(argv[4]):32u);

    QRat Rats;
    QPol QQX(Rats,'X');
    QPol::Element QP;
    std::clog << "# Enter quotient polynomial\n#\t(degree, then list of decreasing degree coefficients):" << std::endl;
    QQX.read(std::cin, QP);
    QQX.write(std::clog << "# Quotient polynomial: ", QP) << " = 0" << std::endl;
    
            

    PMods QQXm(QQX, QP);
    


    QMstream ls(QQXm, left), rs(QQXm, right), ss(QQXm, product);
    Matrix L(ls), R(rs), P(ss);

	assert(L.coldim() == R.coldim());
	assert(L.coldim() == P.rowdim());


#ifdef VERBATIM_PARSING
    L.write(std::clog << "L:=",FileFormat::Maple) << ';' << std::endl;
    R.write(std::clog << "R:=",FileFormat::Maple) << ';' << std::endl;
    P.write(std::clog << "P:=",FileFormat::Maple) << ';' << std::endl;
    std::clog << std::string(20,'#') << std::endl;
#endif

        // =============================================
        // Random inputs
    Givaro::GivRandom generator;
    Givaro::Integer::seeding(generator.seed());

    QVector va(QQXm,L.rowdim()), ua(QQXm,L.coldim());
    QVector vb(QQXm,R.rowdim()), ub(QQXm,R.coldim());
    for(auto iter=ua.begin(); iter!=ua.end(); ++iter) QQX.random(generator, *iter, bitsize);
    for(auto iter=ub.begin(); iter!=ub.end(); ++iter) QQX.random(generator, *iter, bitsize);

        // =============================================
        // Compute matrix product via the HM algorithm
    L.apply(va,ua);
    R.apply(vb,ub);

    QVector wc(QQXm,P.rowdim()), vc(QQXm,P.coldim());
    for(size_t i=0; i<vc.size(); ++i) QQXm.mul(vc[i],va[i],vb[i]);

    P.apply(wc,vc);


        // =============================================
        // Compute the matrix product directly
    size_t n(L.coldim()>>1);
    LinBox::DenseMatrix<PMods> Ma(QQXm,n,n), Mb(QQXm,n,n), Mc(QQXm,n,n);

    for(size_t i=0; i<ua.size(); ++i) Ma.setEntry(i/n,i%n,ua[i]);
    for(size_t i=0; i<ub.size(); ++i) Mb.setEntry(i/n,i%n,ub[i]);
    for(size_t i=0; i<wc.size(); ++i) Mc.setEntry(i/n,i%n,wc[i]);

    LinBox::MatrixDomain<PMods> BMD(QQXm);
    LinBox::DenseMatrix<PMods> Rc(QQXm,n,n);
    BMD.mul(Rc,Ma,Mb); // Direct matrix multiplication


    if (BMD.areEqual (Rc,Mc))
        std::clog <<"# \033[1;32mOK : correct Matrix-Multiplication!\033[0m" << std::endl;
    else{
        std::cerr << "# \033[1;31m****** ERROR, not a MM algorithm******\033[0m"
                  << std::endl;

        ua.write(std::clog << "Ua:=", FileFormat::Maple ) << ';' << std::endl;
        va.write(std::clog << "Va:=", FileFormat::Maple ) << ';' << std::endl;
        ub.write(std::clog << "Ub:=", FileFormat::Maple ) << ';' << std::endl;
        vb.write(std::clog << "Vb:=", FileFormat::Maple ) << ';' << std::endl;
        vc.write(std::clog << "Vc:=", FileFormat::Maple ) << ';' << std::endl;
        wc.write(std::clog << "wc:=", FileFormat::Maple ) << ';' << std::endl;

        Ma.write(std::clog << "Ma:=", FileFormat::Maple) << ';' << std::endl;
        Mb.write(std::clog << "Mb:=", FileFormat::Maple) << ';' << std::endl;
        Mc.write(std::clog << "Mc:=", FileFormat::Maple) << ';' << std::endl;
        Rc.write(std::clog << "Rc:=", FileFormat::Maple) << ';' << std::endl;
    }

    return 0;
}
