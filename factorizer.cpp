// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library
 * Reference:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *     Strassen's algorithm is not optimally accurate
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/

#include "plinopt_sparsify.h"

	// Replace row i of A, by row j of B
inline Matrix& setRow(Matrix& A, size_t i, const Matrix& B, size_t j, const QRat& QQ) {
    A[i].resize(0);
    for(const auto& iter: B[j]) {
        if (! QQ.isZero(iter.second)) A.setEntry(i,iter.first,iter.second);
    }
    return A;
}


int backSolver(Matrix& CoB, Matrix& Res, const Matrix& iM, const QRat& QQ) {

    Matrix M(QQ, iM.rowdim(), iM.coldim()); matrixCopy(M,iM,QQ);
    const size_t r(M.rowdim()), n(M.coldim()), k(CoB.rowdim());
    const size_t s(r-k);
    assert( Res.rowdim() == r );  // Res is (r x k)
    assert( Res.coldim() == k );
    assert( CoB.coldim() == n );  // CoB is (k x n)

    LinBox::Permutation<QRat> T(QQ,r);

    Matrix A2(QQ,s,n);
    size_t t(0);
    for(size_t i=0; i<n; ++i) {
        for(size_t j=i+t;j<r;++j) {
            setRow(CoB, i, M, j, QQ);
            size_t r;
            if (rank(r,CoB) == (i+1)) {
                if (i != j) {
                    T.permute(i,j);
                    std::swap(M[i],M[j]);
                }
                break;
            }
        }
    }
    for(size_t j=n; j<k; ++j) {
        setRow(CoB,j, M, j, QQ);
    }
    for(size_t j=k+t; j<r; ++j)
        setRow(A2,t++, M, j, QQ);

#ifdef VERBATIM_PARSING
    T.write(std::clog << "# Initial permutation: ") << std::endl;
    CoB.write(std::clog << "# Full row rank: ", FileFormat::Pretty) << std::endl;
    A2.write(std::clog << "# Free profile : ", FileFormat::Pretty) << std::endl;
#endif

    Matrix U(QQ,n,k); Transpose(U, CoB); // U is (n x k)
    Matrix B(QQ,n,s); Transpose(B, A2); // B is (n x s)

    Givaro::Rational Det;
    size_t Rank;
    Matrix L(QQ, n, n);
    LinBox::Permutation<QRat> P(QQ,n);
    LinBox::Permutation<QRat> Q(QQ,k);

    static LinBox::GaussDomain<QRat> GD(QQ);
    GD.QLUPin(Rank, Det, P, L, U, Q, n, k );


#ifdef VERBATIM_PARSING
    P.write(std::clog << "P: ") << std::endl;
    L.write(std::clog << "L: ", FileFormat::Pretty) << std::endl;
    U.write(std::clog << "U: ", FileFormat::Pretty) << std::endl;
    Q.write(std::clog << "Q: ") << std::endl;
#endif

    for(size_t i=0; i<k; ++i) Res.setEntry(i,i,QQ.one);

    DenseMatrix TY(QQ,s,n); Transpose(TY,B);
    QVector x(QQ,k), w(QQ,k), bi(QQ,n);
    for(size_t i=0; i<s; ++i) {
        for(size_t j=0; j<n; ++j)
            bi[j] = TY[i][j];
        GD.solve(x, w, r, P, L, U, Q, bi);
        setRow(Res, k+i, x, QQ);
    }

    DenseMatrix R(QQ, r, k);
    T.solveRight(R, Res);
    dense2sparse(Res, R, QQ);
    return 0;
}

// ============================================================
int Factorizer(std::istream& input, const FileFormat& matformat,
             const size_t selectinnerdim) {

        // ============================================================
        // Read Matrix of Linear Transformation
    QRat QQ;
    QMstream ms(QQ,input);
    Matrix M(ms); M.resize(M.rowdim(),M.coldim());
    size_t sc,sb,sr;
    densityProfile(std::clog << "# Initial profile: ", sc, M)
                             << std::endl;

#ifdef VERBATIM_PARSING
    M.write(std::clog,FileFormat::Pretty) << std::endl;
#endif
    size_t innerdim(selectinnerdim == 0 ? M.coldim() : selectinnerdim);
    if (innerdim >= M.rowdim()) {
        std::cerr << "# \033[1;36mFail: inner dimension too large.\033[0m\n";
        return -1;
    }
    if (innerdim < M.coldim()) {
        std::cerr << "# \033[1;36mFail: inner dimension too small.\033[0m\n";
        return -1;
    }

    Givaro::Timer elapsed;
    Matrix CoB(QQ, innerdim, M.coldim());
    Matrix Res(QQ, M.rowdim(), innerdim);

    elapsed.start();
    backSolver(CoB, Res, M, QQ);
    elapsed.stop();

        // ============================================================
        // Print resulting matrices

        // change of basis to stdout
    densityProfile(std::clog << "# Alternate basis profile: \033[1;36m",
                   sb, CoB) << "\033[0m" << std::endl;
    CoB.write(std::cout, matformat) << std::endl;


        // residuum sparse matrix to stdlog
    densityProfile(std::clog << "# Sparse residuum profile: \033[1;36m",
                   sr, Res) << "\033[0m" << std::endl;
    Res.write(std::clog, matformat)<< std::endl;


        // Final check that we computed a factorization M=Res.CoB
    std::clog << std::string(30,'#') << std::endl;
    consistency(std::clog, M, Res, CoB)
        << " \033[1;36m"
        << sr << " non-zeroes (" << sb << " alt.) instead of " << sc
        << "\033[0m:" << ' ' << elapsed << std::endl;

    return 0;
}


// ============================================================
// Main: select between file / std::cin
//       -k #: sets the inner dimension (default is column dimension)
//       -M/-P/-S/-L: selects the ouput format
int main(int argc, char** argv) {

    FileFormat matformat = FileFormat::Pretty;
    std::string filename;
    size_t innerdim(0);					// will default to columndimension

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0] << " [-h|-M|-P|-S|-c #] [stdin|matrixfile.sms]\n";
            exit(-1);
        }
        else if (args == "-M") { matformat = FileFormat(1); } // Maple
        else if (args == "-S") { matformat = FileFormat(5); } // SMS
        else if (args == "-P") { matformat = FileFormat(8); } // Pretty
        else if (args == "-L") { matformat = FileFormat(12); }// Linalg
        else if (args == "-k") { innerdim = atoi(argv[++i]); }
        else { filename = args; }
    }

    if (filename == "") {
        return Factorizer(std::cin, matformat, innerdim);
    } else {
        std::ifstream inputmatrix(filename);
        return Factorizer(inputmatrix, matformat, innerdim);
    }

    return -1;
}
// ============================================================
