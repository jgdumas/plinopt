// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
// ==========================================================================

/**********************************************************************
 * Computes the growth factors fro different norms
 * Usage:   L.sms R.sms P.sms
 * References:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Towards automated generation of fast and accurate algorithms
 *     for recursive matrix multiplication.
 *     Journal of Symbolic Computation 134:102524, 2026.
 *     (https://hal.science/hal-04995684) ]
 **********************************************************************/

#include "plinopt_norms.h"

// ===============================================================
// argv[1-3]: L.sms R.sms P.sms
int main(int argc, char ** argv) {

    if ((argc <=3) || (std::string(argv[1]) == "-h")) {
        std::clog << "Usage:" << argv[0] << " L.sms R.sms P.sms\n";
        exit(-1);
    }

        // =============================================
        // Reading matrices
	std::ifstream left (argv[1]), right (argv[2]), product(argv[3]);

    using PLinOpt::FileFormat;
    PLinOpt::QRat QQ;
    PLinOpt::QMstream ls(QQ, left), rs(QQ, right), ss(QQ, product);
    PLinOpt::Matrix L(ls), R(rs), P(ss);

    if ( (L.rowdim() != R.rowdim()) || (L.rowdim() != P.coldim()) ) {
        std::cerr << "# \033[1;31m****** ERROR, inner dimension mismatch: "
                  << L.rowdim() << "(.)" << R.rowdim() << '|' << P.coldim()
                  << " ******\033[0m"
                  << std::endl;
        return 2;
    }

#if VERBATIM_PARSING >= 2
    L.write(std::clog << "L:=",FileFormat::Maple) << ';' << std::endl;
    R.write(std::clog << "R:=",FileFormat::Maple) << ';' << std::endl;
    P.write(std::clog << "P:=",FileFormat::Maple) << ';' << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

    PLinOpt::Tricounter mkn(PLinOpt::LRP2MM(L,R,P));
    const size_t& m(std::get<0>(mkn)), k(std::get<1>(mkn)), n(std::get<2>(mkn));

    std::clog <<"# Norms of "
              << m << 'x' << k << 'x' << n
              << " Matrix-Multiplication:" << std::endl;

        // =============================================
        // Different norms

    const double ginfinf(PLinOpt::Ginfinf(L,R,P));
    const double ginf2(PLinOpt::Ginf2(L,R,P));
    const double g2inf(PLinOpt::G2inf(L,R,P));
    const double g22(PLinOpt::G22(L,R,P));
    const double g2(PLinOpt::G2(L,R,P));
    const double q0(PLinOpt::Q0(L,R,P));

    const double qkinfinf(_PLO_Qk_(q0,ginfinf,k));
    const double q1inf2(_PLO_Qk_(q0,ginf2,1));
    const double q12inf(_PLO_Qk_(q0,g2inf,1));
    const double sqrtk(std::sqrt(double(k)));
    const double kth(sqrtk*sqrtk*sqrtk);
    const double qk2inf(_PLO_Qk_(q0,g2inf,kth));
    const double q122(_PLO_Qk_(q0,g22,1));


    std::clog << std::fixed << std::setw(8)
              << "#  \t\tGamma \t\tlog_" << k << std::endl;
    std::clog << "## Ginfinf:\t" << ginfinf << '\t'
              << std::log(ginfinf)/std::log(k) << std::endl;
    std::clog << "## Ginf2:\t" << ginf2 << '\t'
              << std::log(ginf2 )/std::log(k) << std::endl;
    std::clog << "## G2inf:\t" << g2inf << '\t'
              << std::log(g2inf)/std::log(k) << std::endl;
    std::clog << "## G22:\t\t" << g22 << '\t'
              << std::log(g22)/std::log(k) << std::endl;
    std::clog << "## G2:\t\t" << g2 << '\t'
              << std::log(g2)/std::log(k) << std::endl;
    std::clog << "## Q0:\t\t" << q0 << '\t'
              << std::log(q0)/std::log(k) << std::endl;
    std::clog << "## Qkinfinf:\t" << qkinfinf << '\t'
              << std::log(qkinfinf)/std::log(k) << std::endl;
    std::clog << "## Q1inf2:\t" << q1inf2 << '\t'
              << std::log(q1inf2)/std::log(k) << std::endl;
    std::clog << "## Qk12inf:\t" << q12inf << '\t'
              << std::log(q12inf)/std::log(k) << std::endl;
    std::clog << "## Qk2inf:\t" << qk2inf << '\t'
              << std::log(qk2inf)/std::log(k) << std::endl;
    std::clog << "## Q122:\t" << q122 << '\t'
              << std::log(q122)/std::log(k) << std::endl;

    return 0;
}
