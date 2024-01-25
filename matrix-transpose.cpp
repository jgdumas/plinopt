// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
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

int MatrixTranspose(std::istream& input, const LinBox::Tag::FileFormat& matformat) {

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
    AT.write(std::cout, matformat) << std::flush;

    chrono.stop();

#ifdef VERBATIM_PARSING
    std::clog << "# [WRIT]: " << AT.rowdim() << 'x' << AT.coldim() << ' ' << chrono << std::endl;
#endif
    return 0;
}


/* Matrix formats:
  Maple = 1,
  SMS = 5,
  Matlab = 7,
  Pretty = 8,
*/
int main(int argc, char ** argv) {
    LinBox::Tag::FileFormat matformat = LinBox::Tag::FileFormat(5); // SMS output by default

    if (argc <= 1) {
        return MatrixTranspose(std::cin, matformat);
    } else {
        std::string args(argv[1]);
        if (args[0] == '-') {
            if (args[1] == 'h') {
                std::clog << "Usage: " << argv[0] << " [-#] [stdin|matrixfile.sms]\n";
                exit(-1);
            } else {
                matformat = LinBox::Tag::FileFormat(-atoi(argv[1]));
                args = std::string(argv[2]);
            }
        }

        std::ifstream inputmatrix(args);
        if ( inputmatrix ) {
            int rt=MatrixTranspose(inputmatrix, matformat);
            inputmatrix.close();
            return rt;
        }
    }
    return -1;
}
