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
#include <set>
#include "plinopt_sparsify.h"

namespace PLinOpt {

template<typename Field, typename _Arr>
bool isZero(const Field& F, const _Arr& v) {
    for(const auto& it: v) {
        if (! F.isZero(it))
            return false;
    }
    return true;
}

template<typename Field, typename _Arr>
int isCano(const Field& F, const _Arr& v) {
    int loc(-1);
    for(size_t i(0); i<v.size(); ++i) {
        if (! F.isZero(v[i])) {
            if (loc != -1) {
                return -1;
            } else {
                loc = i;
            }
        }
    }
    return loc;
}

template<typename Field>
std::ostream& showOut(std::ostream& out, const Field& QQ,
                      const size_t& i, const typename Field::Element & r) {
    out << (Fsign(QQ,r)<0?'-':'+') << 'o' << i;
    if (notAbsOne(QQ,r)) out << '*' << Fabs(QQ,r);
    return out;
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

template<typename Field, typename Row_t>
std::ostream& showLC(std::ostream& out, const Field& QQ, const Row_t& LC) {
    for(const auto& it: LC) {
        showOut(out, QQ, it.first, it.second);
    }
    return out << ';' << std::endl;
}

template<typename _Mat, typename _Arr>
bool Explore(typename _Mat::Row& LC, _Arr& W, const _Mat& M, const size_t m,
             const _Arr&Coeffs, const size_t level) {
    typedef typename _Mat::Element _Elt;
    const auto& F(M.field());
    if (level>0) {
        for(size_t q(m+1); q<M.rowdim(); ++q) {
            _Elt prevv,currv; F.assign(prevv, F.zero); F.assign(currv, F.zero);
            for(size_t v(0); v<Coeffs.size(); ++v) {
                currv = Coeffs[v]-prevv; prevv = Coeffs[v];
                LC.emplace_back(q,prevv);
                for(const auto& el:M[q]) F.axpyin(W[el.first],currv,el.second);
                if (isZero(F,W))  showLC(std::cout, F, LC);
                else {
                    const int i(isCano(F,W));
                    if (i != -1) showLC(std::cout << "+i" << i, F, LC);
                }
                Explore(LC,W,M,q,Coeffs,level-1);
                LC.pop_back();
            }
            for(const auto& el:M[q]) F.maxpyin(W[el.first],prevv,el.second);
        }
        return true;
    } else
        return false;
}


// ============================================================
//
template<typename Field>
int Depender(std::istream& input, const size_t maxnumcoeff,
             const size_t level, const Field& F) {

    using FMatrix=typename Matrix::template rebind<Field>::other;
    using _Elt=typename Field::Element;
    using FArray=std::vector<_Elt>;
        // ============================================================
        // Read Matrix of Linear Transformation
    QRat QQ;
    QMstream ms(QQ,input);
    Matrix B(ms); B.resize(B.rowdim(),B.coldim());
    FMatrix M(B,F);

#ifdef VERBATIM_PARSING
    M.write(std::clog,FileFormat::Pretty);
#endif

    Givaro::Timer elapsed; elapsed.start();


        // ========================================
        // Try 1, -1 and some coefficients in M
    QArray Coeffs{1, -1};
    for( auto indices = B.IndexedBegin();
         (indices != B.IndexedEnd()) ; ++indices ) {
        augment(Coeffs, indices.value().nume(), QQ);
        augment(Coeffs, indices.value().deno(), QQ);
    }

    for(size_t i=2; Coeffs.size() < maxnumcoeff; ++i) {
        augment(Coeffs, Givaro::Integer(i), QQ);
    }
        // ========================================
        // reduce to at most maxnumcoeff
    if (Coeffs.size()>maxnumcoeff) Coeffs.resize(maxnumcoeff);

    std::vector<_Elt> FCoeffs;
    for(const auto& e: Coeffs) {
        _Elt num,den; F.init(num,e.nume()); F.init(den,e.deno());
        F.divin(num,den);
        if (! F.isZero(num)) {
            if(std::find(FCoeffs.begin(), FCoeffs.end(), num) == FCoeffs.end()) {
                FCoeffs.push_back(num);
            }
        }
    }

    std::clog << "# [DEPND] linear combination coefficients: " << FCoeffs
              << std::endl;

    typename FMatrix::Row LC;
    FArray W(M.coldim(),F.zero);
    for(size_t i(0); i<M.rowdim(); ++i) {
        std::clog << "# [DEPND] o" << i << std::endl;
        LC.emplace_back(i,F.one);
        for(const auto& el:M[i]) F.assign(W[el.first],el.second);
        Explore(LC, W, M, i, FCoeffs, level-1);
        for(const auto& el:M[i]) F.assign(W[el.first],F.zero);
        LC.pop_back();
    }

    elapsed.stop();
    std::clog << "# [DEPND]: " << elapsed << std::endl;
    return 0;
}


} // End of namespace PLinOpt
// ============================================



// ============================================================
//
int DependerQ(std::istream& input, const size_t maxnumcoeff,
             const size_t level, const Givaro::Integer& q) {
    if (Givaro::isZero(q)) {
        PLinOpt::QRat QQ;
        return PLinOpt::Depender(input, maxnumcoeff, level, QQ);
    } else {
        Givaro::Modular<Givaro::Integer> FF(q);
        return PLinOpt::Depender(input, maxnumcoeff, level, FF);
    }
}



// ============================================================
// Main: select between file / std::cin
//       -c #: sets the max number of coefficients per iteration
//       -M/-P/-S: selects the ouput format
//       -l #: states the number of monomials
int main(int argc, char** argv) {
    using PLinOpt::FileFormat;

    std::string filename;
    Givaro::Integer q(0u);
    size_t maxnumcoeff(COEFFICIENT_SEARCH); // default max coefficients
    size_t level(4);

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << " [-h|[-c|-l|-q] #] [stdin|matfile.sms]\n"
                      << "  -c #: max number of coefficients per iteration\n"
                      << "  -l #: maximal number of monomials in the combination\n"
                      << "  -q #: modular generation/check (default is Rationals)\n";
            exit(-1);
        }
        else if (args == "-q") { q = Givaro::Integer(argv[++i]); }
        else if (args == "-c") { maxnumcoeff = atoi(argv[++i]); }
        else if (args == "-l") { level = atoi(argv[++i]); }
        else { filename = args; }
    }

    if (filename == "") {
        return DependerQ(std::cin, maxnumcoeff, level, q);
    } else {
        std::ifstream inputmatrix(filename);
        return DependerQ(inputmatrix, maxnumcoeff, level, q);
    }

    return -1;
}
// ============================================================
