// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************************
 * Swapping columns of the matrix M (first argument):
 *   each row is interpreted as a small matrix ri
 *   each ri is transpozed and then re-vectorized as the new i-th row of M
 * The result is a permutation of the columns of M
 * If present, the second argument is the row dimension of the ri
 *   default is the square root of the row dimension of M
 *   in both cases it should be a divisor of the row dimension of M
 ****************************************************************************/

#include <iostream>
#include <givaro/qfield.h>
#include <givaro/givrational.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/matrix/sparse-matrix.h>

int InputTranspose(std::istream& input,
                   const LinBox::Tag::FileFormat& matformat, const size_t im) {

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

    const size_t m(im==0?std::sqrt(A.coldim()) : im);
    const size_t n(A.coldim()/m);

    LinBox::SparseMatrix<Givaro::QField<Givaro::Rational>,
        LinBox::SparseMatrixFormat::SparseSeq> As(QQ, A.rowdim(), A.coldim());


    for(auto it = A.IndexedBegin(); it != A.IndexedEnd(); ++it) {
        const size_t i(it.colIndex() % m);
        const size_t j( (it.colIndex()-i)/m );
        As.setEntry(it.rowIndex(), i*n+j, it.value());
    }

    chrono.stop();
    std::clog << "# [SWAP]: \033[1;32m" << As.rowdim() << 'x' << As.coldim()
              << ' ' << As.rowdim() << '(' << m << 'x' << n << ')'
              << "^T: " << As.rowdim() << '(' << n << 'x' << m << ')'
              << "\033[0m " << chrono << std::endl;

    chrono.start();
    As.write(std::cout, matformat) << std::flush;

    chrono.stop();

    return 0;
}


int main(int argc, char ** argv) {
	// SMS output
    LinBox::Tag::FileFormat matformat = LinBox::Tag::FileFormat(5);
    size_t m(0); // ri row dimension

    std::string filename;
    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0] << " [stdin|file.sms] [-m #]\n"
                      << "  -m #: row dim. of Matricized rows (default sqrt)\n";
            exit(-1);
        }
        else if (args == "-m") { m = atoi(argv[++i]); }
        else { filename = args; }
    }
    
    if (filename.size() > 0) {
        std::ifstream inputmatrix(filename);
        if ( inputmatrix ) {
            int rt=InputTranspose(inputmatrix, matformat, m);
            inputmatrix.close();
            return rt;
        }
    }
    
    return InputTranspose(std::cin, matformat, m);
}
