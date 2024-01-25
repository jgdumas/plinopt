// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear programs
// Authors: J-G. Dumas, B. Grenet
// ==========================================================================

/****************************************************************
 * Returns an in-place program 
 *          computing the bilinear function in HM representation
 *
 * Usage: L.sms R.sms P.sms [expansion]
 *          L.sms/R.sms/P.sms the 3 HM matrices
 *          expansion: if present intermediate results grouped by 2
 *
 * References: 
 *      [ In-place accumulation of fast multiplication formulae
 *        J-G. Dumas, B. Grenet
 *        https://hal.science/hal-04167499 ]
 *
 * Matrix syntax: SMS format, see:
 *      [Sparse Integer Matrix Collection](https://hpac.imag.fr)
 *		- Starts with: `m n 'R'`
 *		- then: `i j value`
 *		- ends with: `0 0 0`
 ****************************************************************/

// ============================================================
// Define to switch between explicit/in-place operator notation
//#define __INPLOP__
// ============================================================

// ============================================================
// Define to print comments during parsing
//#define VERBATIM_PARSING
// ============================================================

// ============================================================
// Define to print verification codes
//#define INPLACE_CHECKER
// ============================================================

#include "plinopt_inplace.h"




// ===============================================================
// argv[1-3]: L.sms R.sms P.sms 
int main(int argc, char ** argv) {


    if ((argc <=3) || (std::string(argv[1]) == "-h")) {
        std::clog << "Usage: " << argv[0] << " L.sms R.sms P.sms [bitsize]\n";
        exit(-1);
    }
     
        // ==========================
        // Reading matrices
	std::ifstream left (argv[1]), right (argv[2]), product(argv[3]);
    size_t bitsize(argc>4?atoi(argv[4]):32u);

    QRat QQ;
    QMstream ls(QQ, left), rs(QQ, right), ss(QQ, product);
    Matrix A(ls), B(rs), C(ss);

#ifdef VERBATIM_PARSING
    A.write(std::clog << "L:=",LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
    B.write(std::clog << "R:=",LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
    C.write(std::clog << "P:=",LinBox::Tag::FileFormat::Maple) << ';' << std::endl;
    std::clog << std::string(20,'#') << std::endl;
#endif

    Givaro::GivRandom generator;
    QVector va(QQ,A.rowdim()), ua(QQ,A.coldim());
    QVector vb(QQ,B.rowdim()), ub(QQ,B.coldim());
    
    for(auto iter=ua.begin(); iter!=ua.end(); ++iter) {
        QQ.random(generator, *iter, bitsize);
    }
    for(auto iter=ub.begin(); iter!=ub.end(); ++iter) {
        QQ.random(generator, *iter, bitsize);
    }

    A.apply(va,ua);
    B.apply(vb,ub);

    QVector wc(QQ,C.rowdim()), vc(QQ,C.coldim());
    for(size_t i=0; i<vc.size(); ++i) {
        QQ.mul(vc[i],va[i],vb[i]);
    }
    
    C.apply(wc,vc);
    

    size_t n(A.coldim()>>1);
        
    LinBox::DenseMatrix<QRat> Ma(QQ,n,n), Mb(QQ,n,n), Mc(QQ,n,n);
    
    for(size_t i=0; i<ua.size(); ++i) {
        Ma.setEntry(i/n,i%n,ua[i]);
    }
    
    for(size_t i=0; i<ub.size(); ++i) {
        Mb.setEntry(i/n,i%n,ub[i]);
    }
    
    for(size_t i=0; i<wc.size(); ++i) {
        Mc.setEntry(i/n,i%n,wc[i]);
    }
    
    LinBox::DenseMatrix<QRat> Res(QQ,n,n);
    LinBox::MatrixDomain<QRat> BMD(QQ);
    BMD.mul(Res,Ma,Mb);
    
    
    if (BMD.areEqual (Res,Mc))
        std::clog <<"# \033[1;32mOK : correct Matrix-Multiplication!\033[0m" << std::endl;
    else{
        std::cerr << "# \033[1;31m****** ERROR, not a MM algorithm******\033[0m"
                  << std::endl;

        ua.write(std::clog << "Ua: ") << std::endl;
        va.write(std::clog << "Va: ") << std::endl;
        ub.write(std::clog << "Ub: ") << std::endl;
        vb.write(std::clog << "Vb: ") << std::endl;
        vc.write(std::clog << "Vc: ") << std::endl;
        wc.write(std::clog << "wc: ") << std::endl;

        Ma.write(std::clog << "Ma: ", LinBox::Tag::FileFormat::Pretty) << std::endl;
        Mb.write(std::clog << "Mb: ", LinBox::Tag::FileFormat::Pretty) << std::endl;
        Mc.write(std::clog << "Mc: ", LinBox::Tag::FileFormat::Pretty) << std::endl;
        Res.write(std::clog << "Re: ", LinBox::Tag::FileFormat::Pretty) << std::endl;
    }

    return 0;
}
