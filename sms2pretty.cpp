// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * Pretty print SMS matrix
 * Matrix syntax: SMS format, see:
 *                [Sparse Integer Matrix Collection](https://hpac.imag.fr)
 *				- Starts with: `m n 'R'`
 *				- then: `i j value`
 *				- ends with: `0 0 0`
 ****************************************************************/

#include <linbox/linbox-config.h>
#include <iostream>
#include <givaro/qfield.h>
#include <givaro/givrational.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/matrix/sparse-matrix.h>

using namespace LinBox;


int PrettyPrint(std::istream& input) 
{
    Givaro::QField<Givaro::Rational> QQ;
    MatrixStream<Givaro::QField<Givaro::Rational>> ms( QQ, input );
    
    Givaro::Timer chrono; chrono.start();
    SparseMatrix<Givaro::QField<Givaro::Rational>, SparseMatrixFormat::SparseSeq> A ( ms );
    chrono.stop();
    
    size_t nnz(0);
    for(auto iter=A.rowBegin(); iter != A.rowEnd(); ++iter) nnz += iter->size();
    
    std::clog << "[READ]: " << A.rowdim() << 'x' << A.coldim() << ' ' << nnz << ' ' << chrono << std::endl;
    
    A.write(std::cout, Tag::FileFormat::Pretty);
    
    return 0;
}

    
    
int main(int argc, char ** argv) {
    if (argc <=1) {
        PrettyPrint(std::cin);    
    } else {
        for(size_t i=1; i<argc; ++i) {
            std::ifstream input (argv[i]);
            PrettyPrint(input);
        }
    }
    return 0;
}

    
    
