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
        if (isAbsOne(QQ,r.nume())) {
            out << '/' << Fabs(QQ,r.deno());
        } else {
            out << '*' << Fabs(QQ,r);
        }
    }
    return out;
}

std::ostream& showLC(std::ostream& out, const QRat& QQ, const Matrix::Row& LC) {
    for(const auto& it: LC) {
        showOut(out, QQ, it.first, it.second);
    }
    return out << ';' << std::endl;
}

bool Explore(Matrix::Row& LC, QArray& W, const Matrix& M, const size_t m,
             const std::vector<Givaro::Rational>& Coeffs, const size_t level) {
    const auto& QQ(M.field());
    if (level>0) {
        for(size_t q(m+1); q<M.rowdim(); ++q) {
            Givaro::Rational prevv(QQ.zero), currv(QQ.zero);
            for(size_t v(0); v<Coeffs.size(); ++v) {
                currv = Coeffs[v]-prevv; prevv = Coeffs[v];
                LC.emplace_back(q,prevv);
                for(const auto& el:M[q]) QQ.axpyin(W[el.first],currv,el.second);
                if (isZero(QQ,W)) showLC(std::cout, QQ, LC);
                Explore(LC,W,M,q,Coeffs,level-1);
                LC.pop_back();
            }
            for(const auto& el:M[q]) QQ.maxpyin(W[el.first],prevv,el.second);
        }
        return true;
    } else
        return false;
}


// ============================================================
//
int Depender(std::istream& input, const FileFormat& matformat,
             const size_t maxnumcoeff, const size_t level) {

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
        // Try 1, -1 and some coefficients in M
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
    std::clog << "# [DEPND] linear combination coefficients: " << Coeffs
              << std::endl;
#endif

    Matrix::Row LC;
    QArray W(M.coldim(),QQ.zero);
    for(size_t i(0); i<M.rowdim(); ++i) {
        std::clog << "# [DEPND] o" << i << std::endl;
        LC.emplace_back(i,QQ.one);
        for(const auto& el:M[i]) QQ.assign(W[el.first],el.second);
        Explore(LC, W, M, i, Coeffs, level-1);
        for(const auto& el:M[i]) QQ.assign(W[el.first],QQ.zero);
        LC.pop_back();
    }

    elapsed.stop();
    std::clog << "# [DEPND]: " << elapsed << std::endl;
    return 0;
}


} // End of namespace PLinOpt
// ============================================


// ============================================================
// Main: select between file / std::cin
//       -c #: sets the max number of coefficients per iteration
//       -M/-P/-S: selects the ouput format
//       -l #: states the number of monomials
int main(int argc, char** argv) {
    using PLinOpt::FileFormat;

    FileFormat matformat = FileFormat::Pretty;
    std::string filename;
    size_t maxnumcoeff(COEFFICIENT_SEARCH); // default max coefficients
    size_t level(4);

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << " [-h|-M|-P|-S|-c #|-l #] [stdin|matfile.sms]\n"
                      << "  -c #: max number of coefficients per iteration\n"
                      << "  -l #: maximal number of monomials in the combination\n"
                      << "  -M/-P/-S: selects the ouput format\n";
            exit(-1);
        }
        else if (args == "-M") { matformat = FileFormat(1); } // Maple
        else if (args == "-S") { matformat = FileFormat(5); } // SMS
        else if (args == "-P") { matformat = FileFormat(8); } // Pretty
        else if (args == "-c") { maxnumcoeff = atoi(argv[++i]); }
        else if (args == "-l") { level = atoi(argv[++i]); }
        else { filename = args; }
    }

    if (filename == "") {
        return PLinOpt::Depender(std::cin, matformat, maxnumcoeff, level);
    } else {
        std::ifstream inputmatrix(filename);
        return PLinOpt::Depender(inputmatrix, matformat, maxnumcoeff, level);
    }

    return -1;
}
// ============================================================
