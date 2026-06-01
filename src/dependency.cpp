// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library
 * Reference:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/

#include "plinopt_sparsify.h"

namespace PLinOpt {

bool isZero(const QRat& QQ, const QArray& v) {
    for(const auto& it: v) { 
        if (! QQ.isZero(it))
            return false;
    }
    return true;
}

std::ostream& showOut(std::ostream& out, const QRat& QQ, 
                      const size_t& i, const Givaro::Rational& r) {
    out << (Fsign(QQ,r)<0?'-':'+') << 'o' << i;
    if (notAbsOne(QQ,r)) {
        out << '*' << Fabs(QQ,r);
    }
    return out;
}
 

// ============================================================
// 
int Depender(std::istream& input, const FileFormat& matformat,
             const size_t maxnumcoeff) {

        // ============================================================
        // Read Matrix of Linear Transformation
    QRat QQ;
    QMstream ms(QQ,input);
    Matrix M(ms); M.resize(M.rowdim(),M.coldim());

#ifdef VERBATIM_PARSING
    M.write(std::clog,FileFormat::Pretty);
#endif

    Givaro::Timer elapsed; elapsed.start();
    

        // ========================================
        // Try 0, 1, -1 and some coefficients in M
    std::vector<Givaro::Rational> Coeffs{1, -1};
    for( auto indices = M.IndexedBegin();
         (indices != M.IndexedEnd()) ; ++indices ) {
        augment(Coeffs, indices.value().nume(), QQ);
        augment(Coeffs, indices.value().deno(), QQ);
    }

    for(size_t i=2; Coeffs.size() < maxnumcoeff; ++i) {
        augment(Coeffs, Givaro::Integer(i), QQ);
    }
        // ========================================
        // reduce to at most maxnumcoeff
    if (Coeffs.size()>maxnumcoeff) Coeffs.resize(maxnumcoeff);

#ifdef VERBATIM_PARSING
    std::clog << "# linear combination coefficients: " << Coeffs << std::endl;
#endif

    QArray V(M.rowdim()),W(M.coldim());
    for(size_t i(0); i<M.rowdim(); ++i) {
        QQ.assign(V[i],QQ.one);
        for(size_t j(i+1); j<M.rowdim(); ++j) {
            for(size_t k(0); k<Coeffs.size(); ++k) {
                QQ.assign(V[j], Coeffs[k]);
                for(size_t l(j+1); l<M.rowdim(); ++l) {
                    for(size_t t(0); t<Coeffs.size(); ++t) {
                        QQ.assign(V[l], Coeffs[t]);
                        M.applyTranspose(W,V);
                        if (isZero(QQ,W)) {
#ifdef VERBATIM_PARSING
                            std::clog << "# W3 " << i << ' ' 
                                      << j << ' ' << Coeffs[k] << ' ' 
                                      << l << ' ' << Coeffs[t]
                                      << ':' << isZero(QQ,W) << " --> ";
#endif
                            std::cout << 'o' << i;
                            showOut(std::cout, QQ, j, Coeffs[k]);
                            showOut(std::cout, QQ, l, Coeffs[t]);
                            std::cout << ';' << std::endl;
#ifdef VERBATIM_PARSING
                            std::clog << '\n';
#endif
                        }
                        for(size_t m(l+1); m<M.rowdim(); ++m) {
                            for(size_t u(0); u<Coeffs.size(); ++u) {
                                QQ.assign(V[m], Coeffs[u]);
                                M.applyTranspose(W,V);
                                if (isZero(QQ,W)) {
#ifdef VERBATIM_PARSING
                                    std::clog << "# W4 " << i << ' ' 
                                              << j << ' ' << Coeffs[k] << ' ' 
                                              << l << ' ' << Coeffs[t] << ' ' 
                                              << m << ' ' << Coeffs[u]
                                              << ':' << isZero(QQ,W) << " --> ";
#endif
                                    std::cout << 'o' << i;
                                    showOut(std::cout, QQ, j, Coeffs[k]);
                                    showOut(std::cout, QQ, l, Coeffs[t]);
                                    showOut(std::cout, QQ, m, Coeffs[u]);
                                    std::cout << ';' << std::endl;
#ifdef VERBATIM_PARSING
                                    std::clog << '\n';
#endif
                                }
                            }
                            QQ.assign(V[m],QQ.zero);
                        }
                    }
                    QQ.assign(V[l],QQ.zero);
                }
            }
            QQ.assign(V[j],QQ.zero);
        }
        QQ.assign(V[i],QQ.zero);
    }
    

//         // Final check that we computed a factorization M=Res.CoB
//     std::clog << std::string(30,'#') << std::endl;
//     consistency(std::clog, M, Res, CoB)
//         << " \033[1;36m"
//         << sr << " non-zeroes (" << sb << " alt.) instead of " << sc
//         << "\033[0m:" << ' ' << elapsed << std::endl;

    return 0;
}


} // End of namespace PLinOpt
// ============================================


// ============================================================
// Main: select between file / std::cin
//       -c #: sets the max number of coefficients per iteration
//       -M/-P/-S: selects the ouput format
//       -b #: states the blocking dimension
//       -U [1|0]: uses an initial LU factorization | or not
int main(int argc, char** argv) {
    using PLinOpt::FileFormat;

    FileFormat matformat = FileFormat::Pretty;
    std::string filename;
    size_t maxnumcoeff(COEFFICIENT_SEARCH); // default max coefficients

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << " [-h|-M|-P|-S|-c #] [stdin|matfile.sms]\n"
                      << "  -c #: max number of coefficients per iteration\n"
                      << "  -M/-P/-S: selects the ouput format\n";

            exit(-1);
        }
        else if (args == "-M") { matformat = FileFormat(1); } // Maple
        else if (args == "-S") { matformat = FileFormat(5); } // SMS
        else if (args == "-P") { matformat = FileFormat(8); } // Pretty
        else if (args == "-c") { maxnumcoeff = atoi(argv[++i]); }
        else { filename = args; }
    }

    if (filename == "") {
        return PLinOpt::Depender(std::cin, matformat, maxnumcoeff);
    } else {
        std::ifstream inputmatrix(filename);
        return PLinOpt::Depender(inputmatrix, matformat, maxnumcoeff);
    }

    return -1;
}
// ============================================================
