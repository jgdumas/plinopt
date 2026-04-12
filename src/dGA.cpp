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
#include "plinopt_library.h"
#include "plinopt_polynomial.h"
#include "plinopt_sparsify.h"

// ===============================================================
void usage(const char* prgname, const size_t randomloops) {
    std::clog << "Usage:" << prgname
              << "  [-h|-b #|-m/-q #|-r # # #|-I x|-P P(x)] L.sms R.sms P.sms\n"
              << "  [-b b]: random check with values of size 'bitsize'\n"
              << "  [-m/-q m]: check is modulo (mod) or (mod/2^k) (default no)\n"
              << "  [-r r e s]: check is modulo (r^e-s) or ((r^e-s)/2^k) (default no)\n"
              << "  [-I x]: indeterminate (default is X)\n"
              << "  [-P P(x)]: modular polynomial (default is none)\n"
              << "  [-O #]: randomized search with that many loops (default "
              << randomloops << " loops)\n";


    exit(-1);
}

template<typename Domain>
typename Domain::Element& zoRandomElt(typename Domain::Element& e,
                                      const Domain& D,
                                      Givaro::GivRandom& generator,
                                      const size_t bitsize) {
    D.init(e, (generator() % 3));
    return D.subin(e, D.one);
}
// ===============================================================


namespace PLinOpt {

// ===============================================================
// Generic verification, within Field, of matrices over Base
template<typename _Matrix>
_Matrix& Tensor(_Matrix& T, const _Matrix& A, const _Matrix& B) {
    T.resize(A.rowdim()*B.rowdim(), A.coldim()*B.coldim());
    using Element=typename _Matrix::Element;
    Element tmp; A.field().init(tmp);
    for(auto itA = A.IndexedBegin(); itA != A.IndexedEnd(); ++itA) {
    for(auto itB = B.IndexedBegin(); itB != B.IndexedEnd(); ++itB) {
        T.setEntry(itA.rowIndex()*B.rowdim()+itB.rowIndex(),
                   itA.colIndex()*B.coldim()+itB.colIndex(),
                   A.field().mul(tmp,itA.value(),itB.value()));
    }
    }
    return T;
}

template<typename _Matrix>
size_t nNonZero(const _Matrix& A) {
    size_t nnz(0);
    for(auto row=A.rowBegin(); row != A.rowEnd(); ++row) {
        nnz += row->size(); // # of non-zero elements
    }
    return nnz;
}

template<typename Field>
size_t nNonZero(const LinBox::DenseMatrix<Field>& A) {
    size_t nnz(0);
    for(auto itA = A.IndexedBegin(); itA != A.IndexedEnd(); ++itA) {
        if (! A.field().isZero(itA.value())) ++nnz;
    }
    return nnz;

}

	// Update col j of A, by v
template<typename _Mat, typename _Vector>
inline _Mat& updCol(_Mat& A, size_t j, const _Vector& v) {
    using Field = typename _Mat::Field;
    const Field& F(A.field());
    for(size_t i=0; i<v.size(); ++i)
        if (! F.isZero(v[i])) A.setEntry(i,j,v[i]);
    return A;
}

	// Computes the inverse of A
template<typename _Mat1, typename _Mat2>
inline _Mat1& inverse(_Mat1& T, const _Mat2& A) {
    assert(A.rowdim() == A.coldim());
    const size_t n(A.rowdim());

    using Field = typename _Mat1::Field;
    using FVector = LinBox::DenseVector<Field>;
    const Field& FF(A.field());

    LinBox::GaussDomain<Field> GD(FF);

    T.resize(n,n);
    _Mat1 U(FF); matrixCopy(U, A);

    typename Field::Element Det;
    size_t Rank;
    _Mat1 L(FF, n, n);
    LinBox::Permutation<Field> Q(FF,n);
    LinBox::Permutation<Field> P(FF,n);

    GD.QLUPin(Rank, Det, Q, L, U, P, n, n );

    FVector x(FF,n), w(FF,n), Iv(FF, n);
    for(size_t j=0; j<Iv.size(); ++j) FF.assign(Iv[j],FF.zero);

    for(size_t j=0; j<n; ++j) {
        FF.assign(Iv[j],FF.one);
        GD.solve(x, w, Rank, Q, L, U, P, Iv);
        FF.assign(Iv[j],FF.zero);
        updCol(T, j, x);
    }

    return T;
}

// Random {-1,0,1}, det=1, matrices
template<typename _Mat>
inline _Mat& zoiRandomMatrix(_Mat& M, const size_t bitsize) {
    static thread_local Givaro::GivRandom
        generator(Givaro::BaseTimer::seed());
    const auto& FF(M.field());
    using Field = typename _Mat::Field;
    using Element=typename Field::Element;
    typedef LinBox::DenseVector<Field> FVector;
#ifdef ACTION_FULL_PLUQ
    _Mat L(FF,M.rowdim(),M.coldim());
    for(size_t i=0; i<L.rowdim(); ++i) L.setEntry(i,i,FF.one);
    Element tmp; FF.init(tmp);
    for(size_t i=0;i<L.rowdim();++i) {
        for(size_t j=0;j<i;++j) {
            if (! FF.isZero(zoRandomElt(tmp, FF, generator, bitsize)))
                L.setEntry(i,j, tmp);
        }
    }
    FVector U(FF,M.coldim(),FF.zero), Row(FF,M.rowdim());
    std::vector<size_t> P(M.rowdim()),Q(M.coldim());
    std::iota(P.begin(), P.end(), 0); // Select a random permutation
    std::iota(Q.begin(), Q.end(), 0); // Select a random permutation
    std::shuffle(P.begin(), P.end(),
                 std::default_random_engine(generator()));
    std::shuffle(Q.begin(), Q.end(),
                 std::default_random_engine(generator()));
    for(size_t i=0;i<L.rowdim();++i) {
        FF.assign(U[Q[i]],FF.one);
        for(size_t j=0;j<Q[i];++j) zoRandomElt(U[j], FF, generator, bitsize);
        for(size_t j=Q[i]+1;j<U.size(); ++j) FF.assign(U[j], FF.zero);
        L.apply(Row,U);
        setRow(M, P[i], Row);
    }
#else
        // Only random triangular
    std::vector<size_t> P(M.rowdim()),Q(M.coldim());
    std::iota(P.begin(), P.end(), 0); // Select a random permutation
    std::iota(Q.begin(), Q.end(), 0); // Select a random permutation
    std::shuffle(P.begin(), P.end(),
                 std::default_random_engine(generator()));
    std::shuffle(Q.begin(), Q.end(),
                 std::default_random_engine(generator()));

    for(size_t i=0; i<M.rowdim(); ++i) M.setEntry(P[i],Q[i],FF.one);
    Element tmp; FF.init(tmp);
    for(size_t i=0;i<M.rowdim();++i) {
        for(size_t j=i+1;j<M.rowdim();++j) {
            if (! FF.isZero(zoRandomElt(tmp, FF, generator, bitsize)))
                M.setEntry(P[i],Q[j], tmp);
        }
    }
#endif
    return M;
}

} // End of namespace PLinOpt


// ===============================================================
// Reading matrices over Base
template<typename Base, typename Field>
int deGrooteAction(const Base& BB, const Field& FF, const size_t bitsize,
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
    size_t nnzl(PLinOpt::nNonZero(L)), nnzr(PLinOpt::nNonZero(R)),
        nnzp(PLinOpt::nNonZero(P));
    PLinOpt::Tricounter mkn(PLinOpt::LRP2MM(L,R,P));
    const size_t& m(std::get<0>(mkn)), k(std::get<1>(mkn)), n(std::get<2>(mkn));
    const size_t mk(m*k), kn(k*n), mn(m*n);


    size_t bestnnz(nnzl+nnzr+nnzp);
    const size_t initnnz(bestnnz);
    std::clog << "# Init. nnz: " << bestnnz << '=' << nnzl << '+'
              << nnzr << '+' << nnzp << std::endl;

    FMatrix bestLj(FF,L.rowdim(), L.coldim()),
        bestRg(FF, R.rowdim(), R.coldim()),
        besthP(FF, P.rowdim(), P.coldim());
    PLinOpt::sparse2sparse(bestLj,L);
    PLinOpt::sparse2sparse(bestRg,R);
    PLinOpt::sparse2sparse(besthP,P);
    Givaro::Timer chrono; chrono.start();

#pragma omp parallel for shared(L,R,P,m,k,n,bestnnz,nnzl,nnzr,nnzp)
    for(size_t i=0; i<randomloops; ++i) {
        FMatrix U(FF,m,m), V(FF,k,k), W(FF,n,n);
        PLinOpt::zoiRandomMatrix(U, bitsize);
        PLinOpt::zoiRandomMatrix(V, bitsize);
        PLinOpt::zoiRandomMatrix(W, bitsize);

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
        const size_t nnzlj(PLinOpt::nNonZero(Lj)), nnzrg(PLinOpt::nNonZero(Rg)),
            nnzhp(PLinOpt::nNonZero(hP));
        const size_t lbnnz(nnzlj+nnzrg+nnzhp);





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
        bool show(false);
        if (nnzlj<=nnzl) {
            show = true;
            nnzl = nnzlj;
        }
        if (nnzrg<=nnzr) {
            show = true;
            nnzr = nnzrg;
        }
        if (nnzhp<=nnzp) {
            show = true;
            nnzp = nnzhp;
        }
        if (lbnnz<bestnnz) {
            bestnnz = lbnnz;
            PLinOpt::dense2sparse(bestLj, Lj);
            PLinOpt::dense2sparse(bestRg, Rg);
            PLinOpt::dense2sparse(besthP, hP);
        }
        if(show)
            std::clog << "# Found nnz: " << lbnnz << '=' << nnzlj << '+'
                      << nnzrg << '+' << nnzhp << '|' << bestnnz << std::endl;
        }
    }

    chrono.stop();

    std::clog << "# Search(" << randomloops <<"): " << chrono << std::endl;

    if (bestnnz<initnnz) {
        std::clog << "# \033[1;36mRdcd. nnz: " << bestnnz <<
            '<' << initnnz << "\033[0m" << std::endl;

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
            else if (args == "-O") {
                randomloops = atoi(argv[++i]);
            }
        } else { filenames.push_back(args); }
    }
    if (filenames.size() < 3) { usage(argv[0],randomloops); }


        // Select verification

    if (modulus>0) {
            // Remove powers of 2 for better probabilities
        while( (modulus % 2) == 0 ) { modulus >>=1; };
        if (modulus == 1) modulus = 2;
            // Compute modular verification of rational matrices
        Givaro::Modular<Givaro::Integer> F(modulus);
        return deGrooteAction(QQ, F, bitsize, randomloops, filenames);
    }
    else if (QP.size()>0) {
            // Compute modular verification of polynomial matrices
        QQuo QQXm(QQX, QP);
        return deGrooteAction(QQX, QQXm, bitsize, randomloops, filenames);
    } else
            // Compute rational verification of rational matrices
        return deGrooteAction(QQ, QQ, bitsize, randomloops, filenames);
}
