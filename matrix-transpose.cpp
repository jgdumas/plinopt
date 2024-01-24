// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * Computes the transposed matrix
 * Matrix syntax: SMS format, see:
 *                [Sparse Integer Matrix Collection](https://hpac.imag.fr)
 *				- Starts with: `m n 'R'`
 *				- then: `i j value`
 *				- ends with: `0 0 0`
 ****************************************************************/

#include <iostream>
#include <givaro/qfield.h>
#include <givaro/givrational.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/matrix/sparse-matrix.h>

// ============================================================
// Define to print comments during parsing
//#define VERBATIM_PARSING
// ============================================================

int MatrixTranspose(std::istream& input) {

    Givaro::QField<Givaro::Rational> QQ;
    LinBox::MatrixStream<Givaro::QField<Givaro::Rational>> ms( QQ, input );

    Givaro::Timer chrono; chrono.start();
    LinBox::SparseMatrix<Givaro::QField<Givaro::Rational>,
        LinBox::SparseMatrixFormat::SparseSeq> A ( ms );
    A.resize(A.rowdim(), A.coldim());
    chrono.stop();

#ifdef VERBATIM_PARSING
    std::clog << "# [READ]: " << A.rowdim() << 'x' << A.coldim() << ' ' << chrono << std::endl;
#endif

    chrono.start();
    LinBox::SparseMatrix<Givaro::QField<Givaro::Rational>,
        LinBox::SparseMatrixFormat::SparseSeq> AT(QQ, A.coldim(), A.rowdim());


    for(auto it = A.IndexedBegin(); it != A.IndexedEnd(); ++it)
        AT.setEntry(it.colIndex(),it.rowIndex(), it.value());

    chrono.stop();
#ifdef VERBATIM_PARSING
    std::clog << "# [TRSP]: " << AT.rowdim() << 'x' << AT.coldim() << ' ' << chrono << std::endl;
#endif

    chrono.start();
    AT.write(std::cout, LinBox::Tag::FileFormat::SMS) << std::flush;

    chrono.stop();

#ifdef VERBATIM_PARSING
    std::clog << "# [WRIT]: " << AT.rowdim() << 'x' << AT.coldim() << ' ' << chrono << std::endl;
#endif
    return 0;
}

int main(int argc, char ** argv) {

    if (argc <= 1) {
        return MatrixTranspose(std::cin);
    } else {
        std::ifstream inputmatrix(argv[1]);
        if ( inputmatrix ) {
            int rt=MatrixTranspose(inputmatrix);
            inputmatrix.close();
            return rt;
        }
    }
    return -1;
}
