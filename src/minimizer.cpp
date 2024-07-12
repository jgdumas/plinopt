// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// ==========================================================================

#include "plinopt_sparsify.h"

using Field = Givaro::Modular<double>;
using FMatrix=typename Matrix::template rebind<Field>::other;

void usage(char ** argv) {
    std::clog << "Usage: " << argv[0] << " L.sms R.sms P.sms p\n";
    exit(-1);
}


// ===============================================================
// argv[1-3]: L.sms R.sms P.sms
int main(int argc, char ** argv) {

    if (argc<5) usage(argv);
        // =================================
        // Reading matrices
	std::ifstream left (argv[1]), right (argv[2]), product(argv[3]);
    size_t p(atoi(argv[4]));
    

    const QRat QQ;
    QMstream ls(QQ, left), rs(QQ, right), ps(QQ, product);
    Matrix A(ls), B(rs), C(ps);
    assert( A.rowdim() == B.rowdim() ); assert( A.rowdim() == C.coldim() );

    const Field FF(p);
    FMatrix L(A,FF), R(B,FF), P(C,FF);
    const size_t dimension(L.rowdim());
    
    L.write(std::clog << "L:=Matrix(",FileFormat::Maple) << ");" << std::endl;
    R.write(std::clog << "R:=Matrix(",FileFormat::Maple) << ");" << std::endl;
    P.write(std::clog << "P:=Matrix(",FileFormat::Maple) << ");" << std::endl;


    std::vector<FMatrix> RankOnes;
    
    Field::Element tmp; FF.init(tmp);
    
    for(size_t k=0; k<dimension; ++k) {
        FMatrix M(FF, L.coldim(), R.coldim());
        const auto& Lk(L[k]);
        const auto& Rk(R[k]);
        for(auto& row: Lk) {
            for(auto& col: Rk) {
                FF.mul(tmp, row.second, col.second);
                M.setEntry(row.first,col.first, tmp);
            }
        }
        RankOnes.push_back(M);
    }
    
    for(const auto& it: RankOnes)
        it.write(std::clog << "M:",FileFormat::Pretty) << std::endl;

    LinBox::MatrixDomain<Field> BMD(FF);
    std::vector<FMatrix> Tensor;
    for(size_t i=0; i<P.rowdim(); ++i) {
        FMatrix M(FF, L.coldim(), R.coldim());
        const auto& Pi(P[i]);
        for(auto& row: Pi) {
            FMatrix Ri(FF); matrixCopy(Ri,RankOnes[row.first]);
            BMD.mulin(Ri, row.second);
            BMD.addin(M, Ri);
        }
        Tensor.push_back(M);
    }
    
    for(const auto& it: Tensor)
        it.write(std::clog << "T:",FileFormat::Pretty) << std::endl;

    return 0;
}
