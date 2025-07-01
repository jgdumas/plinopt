// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
// ==========================================================================

/**********************************************************************
 * Reduces negative values by swapping signs of rows, 2 among L, R, P^T
 * Usage:   L.sms R.sms P.sms
 * outputs:	L.neg.sms R.neg.sms P.neg.sms
 * References:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
 *   [ J-G. Dumas, B. Grenet;
 *     In-place accumulation of fast multiplication formulae
 *     ISSAC 2024, Raleigh, NC USA, pp. 16-25.
 *     (https://hal.science/hal-04167499) ]
 **********************************************************************/

#include "plinopt_library.h"
#include <filesystem>

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
        // Optimizing negative values
    std::clog <<"# Optimizing negative values for "
              << m << 'x' << k << 'x' << n
              << " Matrix-Multiplication..." << std::endl;

    Matrix Pt(QQ); Transpose(Pt,P);
    size_t Gn(0), Nn(0), Sn(0);

    for(size_t i=0; i<L.rowdim(); ++i) {
        Sn += L[i].size(); Sn += R[i].size(); Sn += Pt[i].size();

        size_t Lnegs(0), Rnegs(0), Pnegs(0);
        for(const auto& it: L[i]) if (it.second < QQ.zero) ++Lnegs;
        for(const auto& it: R[i]) if (it.second < QQ.zero) ++Rnegs;
        for(const auto& it: Pt[i]) if (it.second < QQ.zero) ++Pnegs;

        size_t None(Lnegs+Rnegs+Pnegs); Gn += None;
        size_t NLR(L[i].size()-Lnegs+R[i].size()-Rnegs+Pnegs);
        size_t NLP(L[i].size()-Lnegs+Rnegs+Pt[i].size()-Pnegs);
        size_t NRP(Lnegs+R[i].size()-Rnegs+Pt[i].size()-Pnegs);


        std::clog << "# LRP" << '<' << i << '>'
                  << ' ' << Lnegs << '/' << L[i].size()
                  << ' ' << Rnegs << '/' << R[i].size()
                  << ' ' << Pnegs << '/' << Pt[i].size()
                  << "  -->  "
                  << None << ':' << NLR << ',' << NLP << ',' << NRP;

        if ((NLR < None) && (NLR <= NLP) && (NLR <= NRP)) {
            std::clog << "  -->  swap signs L[" << i << "] & R[" << i << ']';
            for(auto& it: L[i]) QQ.negin(it.second);
            for(auto& it: R[i]) QQ.negin(it.second);
            Nn += NLR;
        } else
        if ((NLP < None) && (NLP < NLR) && (NLP <= NRP)) {
            std::clog << "  -->  swap signs L[" << i << "] & Pt[" << i << ']';
            for(auto& it: L[i]) QQ.negin(it.second);
            for(auto& it: Pt[i]) QQ.negin(it.second);
            Nn += NLP;
        } else
        if ((NRP < None) && (NRP < NLP) && (NRP < NLR)) {
            std::clog << "  -->  swap signs R[" << i << "] & Pt[" << i << ']';
            for(auto& it: R[i]) QQ.negin(it.second);
            for(auto& it: Pt[i]) QQ.negin(it.second);
            Nn += NRP;
        } else
            Nn += None;

        std::clog << std::endl;
    }


    Transpose(P, Pt);

        // =============================================
        // Writing matrices
	std::ofstream
        oleft(std::filesystem::path(argv[1]).replace_extension(".neg.sms")),
        oright(std::filesystem::path(argv[2]).replace_extension(".neg.sms")),
        opost(std::filesystem::path(argv[3]).replace_extension(".neg.sms"));

    L.write(oleft,FileFormat(5));
    R.write(oright,FileFormat(5));
    P.write(opost,FileFormat(5));


       // =============================================
    std::clog << "# \033[1;32m" << Nn << '/' << Sn << "\tnegative\tinstead of "
              << Gn << '/' << Sn << "\033[0m" << std::endl;

    return 0;
}
