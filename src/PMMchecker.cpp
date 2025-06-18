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
using Givaro::Integer;
using Givaro::Degree;
using Givaro::Indeter;


typedef Givaro::QField<Givaro::Rational> QRat;

#if 0
typedef LinBox::PolynomialRing<QRat> QPol;
#else


class QPol : public LinBox::PolynomialRing<QRat> {
    typedef LinBox::PolynomialRing<QRat> Parent_t;
public:
    using Parent_t::Element;
    using Parent_t::Element_ptr;
    using Parent_t::ConstElement_ptr;
    using Parent_t::Type_t;
    using Parent_t::Parent_t; // using inherited constructors

    std::istream& read ( std::istream& i, Element& P) const {
        this->init(P); // reset P
        using Givaro::Degree;
        std::string s;
        std::getline(i, s);
        s.erase(std::remove_if(s.begin(), s.end(), ::isspace),s.end());

        std::ostringstream Xstream; Xstream << this->getIndeter();
        std::string X(Xstream.str());


            // Break into units at every '+' or '-'; the limits will be [p,q)
        size_t p = 0, q = p;
        while ( q < s.size() ) {
            for ( q = p + 1; q < s.size() && s[q] != '+' && s[q] != '-'; q++ );
            std::string unit(s.substr( p, q - p ));
            // should be of form "cxn", meaning c times x^n
            // Pick out coefficient and exponent
            QRat::Element c; this->getDomain().assign(c,this->getDomain().one);
            Givaro::Degree n(0);
            const size_t pos(unit.find( X )); // position of char X
            if ( pos == std::string::npos ) { // X not found; pure number
                std::stringstream( unit ) >> c;
            } else { // identify coefficient (c) and exponent (n)
                if ( pos != 0 ) { // pos == 0 would mean default c = 1
                    std::string first = unit.substr( 0, pos );
                    if ( first == "+" ) this->getDomain().assign(c,this->getDomain().one); // just "+" means +1
                    else if ( first == "-" ) this->getDomain().assign(c,this->getDomain().mOne); // just "-" means -1
                    else std::stringstream( first ) >> c;
                }
                n = 1; // n>=1 now
                if ( pos != unit.size() - 1 ) {

                        // monomial c * X^n
                    const std::size_t expo(unit.find_first_of("0123456789"));
                    if (expo != std::string::npos) {
                        std::stringstream( unit.substr( expo ) ) >> n;
                    }
                }
            }

//             this->write(std::clog << "Current P" << P << ':', P) << std::endl;
//             std::clog << "Read deg: " << n << std::endl;
//             std::clog << "Read coef: " << c << std::endl;

            Element monomial; this->init(monomial, Degree(n), c);
            this->addin(P,monomial);
            p = q;
        }

//         this->write(std::clog << "Parent read" << P << ':', P)
//                               << std::endl;

        return i;
    }
};
#endif


typedef Givaro::QuotientDom<QPol> PMods;

typedef LinBox::MatrixStream<PMods> QMstream;
typedef LinBox::SparseMatrix<PMods,
                             LinBox::SparseMatrixFormat::SparseSeq > Matrix;
typedef LinBox::DenseVector<PMods> QVector;




// ===============================================================
// argv[1-3]: L.sms R.sms P.sms
// argv[4]: bitsize
int main(int argc, char ** argv) {
    QRat Rats;
    QPol QQX(Rats, 'X');
    QPol::Element QP;
    std::vector<std::string> files;
    size_t bitsize(32u);

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << " [-h|-I x|-P P(x)|-b #] L.sms R.sms P.sms \n";
            std::clog
                << "  -I x: indeterminate (default is X)\n"
                << "  -P P(x): modular polynomial(default is none)\n"
                << "  -b #: random check with #-large values (default is 32)\n";

            exit(-1);
        }
        else if (args == "-I") { QQX.setIndeter(Indeter(argv[++i][0])); }
        else if (args == "-P") {
            std::stringstream sargs( argv[++i] );
            QQX.read(sargs, QP);
        } else if (args == "-b") { std::stringstream(argv[++i]) >> bitsize; }
        else { files.push_back(args); }
    }

        // =============================================
        // Reading matrices
	std::ifstream left (files[0]), right (files[1]), product(files[2]);

        // =============================================
        // Reading modular polynomial
//     std::clog << "# Enter quotient polynomial in " << QQX.getIndeter()
//               << " :" << std::endl;
//     QQX.read(std::cin, QP);
    QQX.write(QQX.write(std::clog << "# Quotient polynomial over ")
              << " : ", QP) << " = 0" << std::endl;

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
        QQX.write(std::clog <<"# \033[1;32mOK : correct Matrix-Multiplication modulo ",QP) << ".\033[0m" << std::endl;
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
