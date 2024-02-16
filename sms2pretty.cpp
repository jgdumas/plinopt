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
 * References:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *     Strassen's algorithm is not optimally accurate
 *     (https://hal.science/hal-04441653) ]
 *   [ J-G. Dumas, B. Grenet; Jul. 2023
 *     In-place accumulation of fast multiplication formulae
 *     (https://hal.science/hal-04167499) ]
 ****************************************************************/

#include <linbox/linbox-config.h>
#include <iostream>
#include <givaro/qfield.h>
#include <givaro/givrational.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/matrix/sparse-matrix.h>

#include "plinopt_library.h"
using LinBox::Tag::FileFormat;

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
            if ( (j+1) < A.coldim()) out << '&';
        }
        out << "\\\\\n";
    }
    return out << "\\end{bmatrix}";
}
// ============================================================



// ============================================================
// Pretty print SMS matrix, number of non-zero & Frobenius norm
// ============================================================
int PrettyPrint(std::istream& input, const FileFormat& matformat)
{
    Givaro::QField<Givaro::Rational> QQ;
    LinBox::MatrixStream<Givaro::QField<Givaro::Rational>> ms( QQ, input );

        // ====================================================
        // Read Matrix
    Givaro::Timer chrono; chrono.start();
    Matrix A ( ms );
    chrono.stop();

        // ====================================================
        // number of non-zero & Frobenius norm
    size_t nnz(0); Givaro::Rational FNormSq(0);
    for(auto row=A.rowBegin(); row != A.rowEnd(); ++row) {
        nnz += row->size(); // # of non-zero elements
        for(const auto& iter:*row)
            FNormSq += (iter.second)*(iter.second); // Square of Frobenius norm
    }

    std::clog << "# [READ]: \033[1;32m" << A.rowdim() << 'x' << A.coldim() << ' ' << nnz << ' ' << FNormSq << "\033[0m " << chrono << std::endl;

        // ====================================================
        // Special latex for rationals or one of LinBox formats
    if (matformat == FileFormat::LaTeX) {
        latex_print(std::cout, A) << std::endl;
    } else {

        if ( (matformat == FileFormat::Plain) ||
             (matformat == FileFormat::HTML)  ||
             (matformat == FileFormat(6))||
             (matformat == FileFormat::linalg)  ) {
            DenseMatrix B(QQ,A.rowdim(),A.coldim()); matrixCopy(B,A,QQ);
            if (matformat == FileFormat(6))
                std::cout << B << std::endl;
            else
                B.write(std::cout, matformat);
        } else {
            A.write(std::cout, matformat);
        }
    }

    return 0;
}
// ============================================================


// ============================================================
// LinBox Matrix formats:
/* Plain = 0,
 * Maple = 1,
 * HTML = 2,
 * LaTeX = 3,
 * SMS = 5,
 * Dense (column-major MatrixMarket) = 6,
 * Matlab = 7,
 * Pretty = 8,
 * OneBased = 10,
 * MatrixMarket = 11,
 * linalg = 12,
 */
int main(int argc, char ** argv) {
        // ====================================================
        // argv[i]: '-i' selects format 'i', otherwise filename
        // ====================================================
    std::vector<FileFormat> matformats;
    std::vector<std::string> filenames;

    for(int i=1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args[0] == '-') {
            if (args[1] == 'h') {
                std::clog << "Usage: " << argv[0] << " [stdin| ([-#] matrixfile.sms)*]\n";
                exit(-1);
            } else {
                matformats.push_back(FileFormat(-atoi(args.c_str())));
            }
        } else {
            filenames.push_back(args);
        }
    }

        // Pretty output by default
    if (matformats.size() == 0)
        matformats.emplace_back(FileFormat::Pretty);

    if (filenames.size() > 0) {
        for(size_t i(0); i<filenames.size(); ++i) {
            std::ifstream input (filenames[i]);
            PrettyPrint(input,
                        i<matformats.size() ?
                        matformats[i] :
                        matformats.back());
        }
    } else {
        return PrettyPrint(std::cin, matformats.front());
    }

    return 0;
}
