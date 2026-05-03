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
 * Examples:
 * > ./bin/MMchecker data/2x2x2_7_Strassen_{L,R,P}.sms
 * > ./bin/MMchecker data/2x2x2_7_DPS-accurate_{L,R,P}.sms -m 513083
 * > ./bin/MMchecker data/2x2x2_7_DPS-accurate_{L,R,P}.sms -r 1013 2 3
 * > ./bin/MMchecker data/2x2x2_7_DPS-accurate-X_{L,R,P}.sms -P "X^2-3"
 *
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

#include <filesystem>
#include <givaro/modular.h>
#include <givaro/givquotientdomain.h>
#include "plinopt_polynomial.h"
#include "plinopt_library.h"
#include "plinopt_sparsify.h"
#include "plinopt_optimize.h"

// ===============================================================
void usage(const char* prgname, const size_t randomloops) {
    std::clog << "Usage:" << prgname
              << "  [-h|-b #|-m/-q #|-r # # #|-I x|-P P(x)] L.sms R.sms P.sms\n"
              << "  [-b b]: random check with values of size 'bitsize'\n"
              << "  [-m/-q m]: check is modulo (mod) or (mod/2^k) (default no)\n"
              << "  [-r r e s]: check is modulo (r^e-s) or ((r^e-s)/2^k) (default no)\n"
              << "  [-I x]: indeterminate (default is X)\n"
              << "  [-P P(x)]: modular polynomial (default is none)\n"
              << "  [-s|-z]: search sparser|faster (default is sparser)\n"
              << "  [-O #]: randomized search with that many loops (default "
              << randomloops << " loops)\n";


    exit(-1);
}
// ===============================================================


namespace PLinOpt {

// Random {-1,0,1}, det=1, matrices
template<typename _Mat>
inline _Mat& zoiRandomMatrix(_Mat& M) {
    static thread_local Givaro::GivRandom
        generator(Givaro::BaseTimer::seed());
    const auto& FF(M.field());
    using Field = typename _Mat::Field;
    using Element=typename Field::Element;
    typedef LinBox::DenseVector<Field> FVector;

    std::vector<size_t> P(M.rowdim()),Q(M.coldim());
    std::iota(P.begin(), P.end(), 0); // Select a random permutation
    std::iota(Q.begin(), Q.end(), 0); // Select a random permutation
    std::shuffle(P.begin(), P.end(),
                 std::default_random_engine(generator()));
    std::shuffle(Q.begin(), Q.end(),
                 std::default_random_engine(generator()));

    std::vector<bool> D(M.rowdim()); for(auto&& it: D) it = bool(generator()&1);

#if defined(ACTION_FULL_PLUQ)
        // Random det 1 matrix
    _Mat L(FF,M.rowdim(),M.coldim());
    for(size_t i=0; i<L.rowdim(); ++i) xL.setEntry(i,i, D[i] ? FF.one : FF.mOne);
    Element tmp; FF.init(tmp);
    for(size_t i=0;i<L.rowdim();++i) {
        for(size_t j=0;j<i;++j) {
            if (! FF.isZero(zoRandomElt(tmp, FF, generator)))
                L.setEntry(i,j, tmp);
        }
    }
    FVector U(FF,M.coldim(),FF.zero), Row(FF,M.rowdim());
    for(size_t i=0;i<L.rowdim();++i) {
        FF.assign(U[Q[i]],FF.one);
        for(size_t j=0;j<Q[i];++j) zoRandomElt(U[j], FF, generator);
        for(size_t j=Q[i]+1;j<U.size(); ++j) FF.assign(U[j], FF.zero);
        L.apply(Row,U);
        setRow(M, P[i], Row);
    }

#elif defined(ACTION_HOUSEHOLDER)
        // Random orthogonal rank-1 update
    std::vector<Element> u(M.rowdim()), v(M.rowdim());
    Element dutu; FF.init(dutu);
    zoRandomVect(dutu, u, FF, generator);
    if (FF.isZero(dutu)) {
        for(size_t i=0; i<M.rowdim(); ++i) {
            M.setEntry(P[i],Q[i], D[i]?FF.one:FF.mOne );
        }
    } else {
        Element coef, tmp; FF.init(coef); FF.init(tmp);

        FF.invin(dutu);
        FF.add(tmp,dutu,dutu);
        FF.neg(dutu,tmp); // -2/u^Tu

        for(size_t i=0; i<M.rowdim(); ++i) {
            for(size_t j=0; j<M.coldim(); ++j) {
                FF.mul(coef,u[i],u[j]);
                FF.mulin(coef,dutu);
                if (i==j) FF.addin(coef,FF.one);
                if (!FF.isZero(coef))  M.setEntry(P[i],Q[j],
                                                  D[i]? coef : FF.negin(coef));
            }
        }
    }

#else
        // Only random triangular
    for(size_t i=0; i<M.rowdim(); ++i)
        M.setEntry(P[i],Q[i], D[i] ? FF.one : FF.mOne );
    Element tmp; FF.init(tmp);
    for(size_t i=0;i<M.rowdim();++i) {
        for(size_t j=i+1;j<M.rowdim();++j) {
            if (! FF.isZero(zoRandomElt(tmp, FF, generator)))
                M.setEntry(P[i],Q[j],tmp);
        }
    }
#endif

    return M;
}



template<typename _Mat>
Pair<size_t> nonzeroes(const _Mat& L, const _Mat& R, const _Mat& P) {
    Pair<size_t> nnz(nonzeroes(L));
    const Pair<size_t> nnzr(nonzeroes(R)), nnzp(nonzeroes(P));
    nnz += nnzr; return nnz += nnzp;
}

} // End of namespace PLinOpt



// Counting non-zero elements
template<int Measure>
struct Operations {
    template<typename _Mat>
    size_t operator()(const _Mat& L, const _Mat& R, const _Mat& P,
                      const size_t silent) {
        const size_t nnzl(PLinOpt::density(L)), nnzr(PLinOpt::density(R)),
            nnzp(PLinOpt::density(P));
        if (! silent) std::clog << nnzl << '+' << nnzr << '+' << nnzp;
        return (nnzl+nnzr+nnzp);
    }
};


// Counting number of SLP operations
template<>
struct Operations<1> {
    template<typename _Mat>
    size_t operator()(const _Mat& L, const _Mat& R, const _Mat& P,
                      const size_t subloops) {
        return nbOperations(L,R,P,subloops);
    }

    protected:
template<typename _Mat>
size_t nbOperations(const _Mat& M, const size_t rl) {
    auto nbops(PLinOpt::naiveOps(M));
    std::ostringstream sout;
    Givaro::Timer global;
    _Mat T(M.field(),M.coldim(),M.rowdim()); PLinOpt::Transpose(T,M);
    PLinOpt::CSEOptimiser(nbops, sout, M.field(), M, T, global, rl);
    return nbops.first+nbops.second;
}

template<typename Field>
size_t nbOperations(const LinBox::DenseMatrix<Field>& A, const size_t rl) {
    LinBox::SparseMatrix<Field> M(A.field(),A.rowdim(),A.coldim());
    PLinOpt::dense2sparse(M,A);
    return nbOperations(M, rl);
}


template<typename _Mat>
size_t nbOperations(const _Mat& L, const _Mat& R, const _Mat& P,
                    const size_t rl) {
    auto nbL(nbOperations(L,rl)),
        nbR(nbOperations(R,rl)),
        nbP(nbOperations(P,rl));
    return nbL+nbR+nbP;
}

};



// ===============================================================
// Reading matrices over Base
template<int Measure> struct Orbiter {
    template<typename Base, typename Field>
    int operator()(const Base& BB, const Field& FF, const size_t bitsize,
                   const size_t randomloops,
                   const std::vector<std::string>& filenames) {
    typedef LinBox::MatrixStream<Base> BMstream;
    typedef LinBox::SparseMatrix<Base,
        LinBox::SparseMatrixFormat::SparseSeq > BMatrix;
    using FMatrix=typename BMatrix::template rebind<Field>::other;
    using Element=typename Field::Element;
    typedef LinBox::DenseVector<Field> FVector;
    typedef LinBox::DenseMatrix<Field> DMatrix;
    using PLinOpt::FileFormat;


        // =============================================
        // Reading matrices
	std::ifstream left (filenames[0]), right (filenames[1]), post(filenames[2]);

    BMstream ls(BB, left), rs(BB, right), ss(BB, post);
    BMatrix BL(ls), BR(rs), BP(ss);
    FMatrix L(BL,FF), R(BR,FF), P(BP,FF);

    if ( (L.rowdim() != R.rowdim()) || (L.rowdim() != P.coldim()) ) {
        std::cerr << "# \033[1;31m****** ERROR, inner dimension mismatch: "
                  << L.rowdim() << "(.)" << R.rowdim() << '|' << P.coldim()
                  << " ******\033[0m"
                  << std::endl;
//         return 2;
    }

#ifdef VERBATIM_PARSING
    L.write(std::clog << "L:=",FileFormat::Maple) << ';' << std::endl;
    R.write(std::clog << "R:=",FileFormat::Maple) << ';' << std::endl;
    P.write(std::clog << "P:=",FileFormat::Maple) << ';' << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

    PLinOpt::MMchecker(FF, bitsize, L, R, P);

    PLinOpt::Tricounter mkn(PLinOpt::LRP2MM(L,R,P));
    const size_t& m(std::get<0>(mkn)), k(std::get<1>(mkn)), n(std::get<2>(mkn));
    const size_t mk(m*k), kn(k*n), mn(m*n);

    const size_t subloops( randomloops>16 ? randomloops>>4: 1);
    const size_t initopt = Operations<Measure>()(L,R,P,subloops);
    size_t bestopt(initopt);
    const auto initnz(PLinOpt::nonzeroes(L,R,P));
    auto bestnz(initnz);
    std::clog << "# Init. ops: " << initopt << ", " << initnz << std::endl;

    FMatrix bestLj(FF,L.rowdim(), L.coldim()),
        bestRg(FF, R.rowdim(), R.coldim()),
        besthP(FF, P.rowdim(), P.coldim());
    PLinOpt::sparse2sparse(bestLj,L);
    PLinOpt::sparse2sparse(bestRg,R);
    PLinOpt::sparse2sparse(besthP,P);
    Givaro::Timer chrono; chrono.start();

#pragma omp parallel for shared(L,R,P,m,k,n,bestopt,subloops,bestLj,bestRg,besthP)
    for(size_t i=0; i<randomloops; ++i) {
        FMatrix U(FF,m,m), V(FF,k,k), W(FF,n,n);
        PLinOpt::zoiRandomMatrix(U);
        PLinOpt::zoiRandomMatrix(V);
        PLinOpt::zoiRandomMatrix(W);

        FMatrix J(FF,mk,mk), G(FF,kn,kn), H(FF,mn,mn);
        FMatrix iVT(FF,k,k); PLinOpt::inverseTranspose(iVT, V);
        FMatrix iU(FF,m,m); PLinOpt::inverse(iU, U);
        FMatrix iW(FF,n,n); PLinOpt::inverse(iW, W);

        PLinOpt::Tensor(J,iU,V);
        PLinOpt::Tensor(G,iVT,W);
        PLinOpt::Tensor(H,U,iW);

        DMatrix Lj(FF,L.rowdim(), L.coldim()),
            Rg(FF, R.rowdim(), R.coldim()),
            hP(FF, P.rowdim(), P.coldim());
        LinBox::MatrixDomain<Field> BMD(FF);
        BMD.mul(Lj, L, J);
        BMD.mul(Rg, R, G);
        BMD.mul(hP, H, P);

        const size_t lbopt = Operations<Measure>()(Lj,Rg,hP,subloops);

        if (lbopt<=bestopt) {
            const auto lbnz(PLinOpt::nonzeroes(Lj,Rg,hP));
            if ((lbopt<bestopt) || (lbnz.first<bestnz.first) || (lbnz.second<bestnz.second)) {
#pragma omp critical
                {
#ifdef VERBATIM_PARSING
                    U.write(std::clog << "U:=",FileFormat::Maple) << ';' << std::endl;
                    V.write(std::clog << "V:=",FileFormat::Maple) << ';' << std::endl;
                    W.write(std::clog << "W:=",FileFormat::Maple) << ';' << std::endl;
                    J.write(std::clog << "J:=",FileFormat::Maple) << ';' << std::endl;
                    G.write(std::clog << "G:=",FileFormat::Maple) << ';' << std::endl;
                    H.write(std::clog << "H:=",FileFormat::Maple) << ';' << std::endl;
                    std::clog << "# Curr. nnz: " << lbnnz << '=' << nnzlj << '+'
                              << nnzrg << '+' << nnzhp << std::endl;
#endif

                std::clog << "# Found opt: " << lbopt << (lbopt<bestopt?'<':'=') << bestopt
                          << '\t' << lbnz << '&' << bestnz
                          << "\t[" << i << '/' << omp_get_thread_num() << ']'
                          << std::endl;
                bestopt = lbopt;
                bestnz = lbnz;
                PLinOpt::dense2sparse(bestLj, Lj);
                PLinOpt::dense2sparse(bestRg, Rg);
                PLinOpt::dense2sparse(besthP, hP);
                }
            }
        }
    }

    chrono.stop();

    std::clog << "# Search(" << randomloops <<"): " << chrono << std::endl;

    if ((bestopt<=initopt) &&
        ((bestopt<initopt) || (bestnz.first<initnz.first) || (bestnz.second<initnz.second)) ) {
        std::clog << "# \033[1;36mRdcd. opt: " << bestopt << '=';
        Operations<0>()(bestLj,bestRg,besthP,0);
        std::clog << '<' << initopt << "\033[0m" << std::endl;

            // =============================================
            // Writing matrices
        std::ofstream
            oleft(std::filesystem::path(filenames[0])
                  .replace_extension(".nnz.sms")),
            oright(std::filesystem::path(filenames[1])
                   .replace_extension(".nnz.sms")),
            opost(std::filesystem::path(filenames[2])
              .replace_extension(".nnz.sms"));

        bestLj.write(oleft,FileFormat(5));
        bestRg.write(oright,FileFormat(5));
        besthP.write(opost,FileFormat(5));

#ifdef VERBATIM_PARSING
        bestLj.write(std::clog << "Lj:=",FileFormat::Maple) << ';' << std::endl;
        bestRg.write(std::clog << "Rg:=",FileFormat::Maple) << ';' << std::endl;
        besthP.write(std::clog << "hP:=",FileFormat::Maple) << ';' << std::endl;
#endif
        PLinOpt::MMchecker(FF, bitsize, bestLj, bestRg, besthP);
    }

    return 0;
    }
};



// ===============================================================
int main(int argc, char ** argv) {

    size_t bitsize(32u);
    size_t randomloops(DORANDOMSEARCH?DEFAULT_RANDOM_LOOPS:1);
    Givaro::Integer modulus(0u);
    PLinOpt::QRat QQ;
    using QPol = PLinOpt::PRing<PLinOpt::QRat>;
    using QQuo = Givaro::QuotientDom<QPol>;
    QPol QQX(QQ, 'X');
    QPol::Element QP;
    std::vector<std::string> filenames;

    size_t optim(0); // 0 for sparsity, 1 for nbops

        // Parse arguments

    for(int i=1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args[0] == '-') {
            if (args[1] == 'h') { usage(argv[0],randomloops); }
            else if (args[1] == 'b') { bitsize = atoi(argv[++i]); }
            else if ((args[1] == 'm') || (args[1] == 'q')) { modulus = Givaro::Integer(argv[++i]); }
            else if (args[1] == 'r') {
                Givaro::Integer r(argv[++i]);
                size_t e(atoi(argv[++i]));
                Givaro::Integer s(argv[++i]);
                modulus = pow(r,e)-s;
            }
            else if (args[1] == 'I') {
                QQX.setIndeter(Givaro::Indeter(argv[++i][0]));
            }
            else if (args[1] == 'P') {
                std::stringstream sargs( argv[++i] );
                QQX.read(sargs, QP);
            }
            else if (args[1] == 'O') {
                randomloops = atoi(argv[++i]);
            }
            else if (args[1] == 's') {
                optim = 0;
            }
            else if (args[1] == 'z') {
                optim = 1;
            }
        } else { filenames.push_back(args); }
    }
    if (filenames.size() < 3) { usage(argv[0],randomloops); }


        // Select verification

    if (modulus>0) {
//		// Remove powers of 2 for better probabilities
        while( (modulus % 2) == 0 ) { modulus >>=1; };
        if (modulus == 1) modulus = 2;
//		// Compute modular verification of rational matrices
        Givaro::Modular<Givaro::Integer> F(modulus);
        return Orbiter<0>()(QQ, F, bitsize, randomloops, filenames);
    } else if (QP.size()>0) {
            // Compute modular verification of polynomial matrices
        QQuo QQXm(QQX, QP);
        return Orbiter<0>()(QQX, QQXm, bitsize, randomloops, filenames);
    } else {
		// Compute rational verification of rational matrices
        if (optim == 1) {
            return Orbiter<1>()(QQ, QQ, bitsize, randomloops, filenames);
        } else {
            return Orbiter<0>()(QQ, QQ, bitsize, randomloops, filenames);
        }
    }
}
