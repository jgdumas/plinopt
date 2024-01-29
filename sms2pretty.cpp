// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
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


// ============================================================
// Special latex formatting for matrices of rationals
// ============================================================
template<typename Matrix>
std::ostream& latex_print(std::ostream& out, const Matrix& A) {
    out << "\\begin{bmatrix}\n";
    Givaro::Rational tmp;
    for(size_t i=0; i<A.rowdim();++i) {
        for(size_t j=0; j<A.coldim(); ++j) {
            A.getEntry(tmp,i,j);
            if (tmp.nume()<0) {
                out << '-'; tmp=-tmp;
            }
            if (! isOne(tmp.deno())) {
                out << "\\frac{" << tmp.nume() << "}{" << tmp.deno() << '}';
            } else {
                out << tmp;
            }
            out << '&';
        }
        out << "\\\\\n";
    }
    return out << "\\end{bmatrix}";
}
// ============================================================



// ============================================================
// Pretty print SMS matrix, number of non-zero & Frobenius norm
// ============================================================
int PrettyPrint(std::istream& input, const Tag::FileFormat& matformat)
{
    Givaro::QField<Givaro::Rational> QQ;
    MatrixStream<Givaro::QField<Givaro::Rational>> ms( QQ, input );

        // ====================================================
        // Read Matrix
    Givaro::Timer chrono; chrono.start();
    SparseMatrix<Givaro::QField<Givaro::Rational>,
        SparseMatrixFormat::SparseSeq> A ( ms );
    chrono.stop();

        // ====================================================
        // number of non-zero & Frobenius norm
    size_t nnz(0); Givaro::Rational FNormSq(0);
    for(auto row=A.rowBegin(); row != A.rowEnd(); ++row) {
        nnz += row->size(); // # of non-zero elements
        for(const auto& iter:*row)
            FNormSq += (iter.second)*(iter.second); // Square of Frobenius norm
    }

    std::clog << "[READ]: " << A.rowdim() << 'x' << A.coldim() << ' ' << nnz << ' ' << FNormSq << ' ' << chrono << std::endl;

        // ====================================================
        // Special latex for rationals or one of LinBox formats
    if (matformat == Tag::FileFormat::LaTeX) {
        latex_print(std::cout, A) << std::endl;
    } else {
        A.write(std::cout, matformat);
    }

    return 0;
}
// ============================================================


// ============================================================
// LinBox Matrix formats:
/* Maple = 1,
 * LaTeX = 3,
 * SMS = 5,
 * Matlab = 7,
 * Pretty = 8,
 */
int main(int argc, char ** argv) {
        // ====================================================
        // argv[i]: '-i' selects format 'i', otherwise filename

        // ====================================================
        // Pretty output by default
    Tag::FileFormat matformat = Tag::FileFormat::Pretty;

    if (argc <=1) {
        PrettyPrint(std::cin, matformat);
    } else {
        for(int i=1; i<argc; ++i) {
            std::string args(argv[i]);
            if (args[0] == '-') {
                if (args[1] == 'h') {
                    std::clog << "Usage: " << argv[0] << " [-#] [stdin|matrixfile.sms]\n";
                    exit(-1);
                } else {
                    matformat = Tag::FileFormat(-atoi(args.c_str()));
                    args = std::string(argv[++i]);
                }
            }
            std::ifstream input (args);
            PrettyPrint(input, matformat);
        }
    }

    return 0;
}
