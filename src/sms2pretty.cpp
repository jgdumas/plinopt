// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
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
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
 *   [ J-G. Dumas, B. Grenet;
 *     In-place accumulation of fast multiplication formulae
 *     ISSAC 2024, Raleigh, NC USA, pp. 16-25.
 *     (https://hal.science/hal-04167499) ]
 ****************************************************************/

#include <linbox/linbox-config.h>
#include <iostream>
#include <givaro/qfield.h>
#include <givaro/givrational.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/matrix/sparse-matrix.h>

#include "plinopt_library.h"
using PLinOpt::FileFormat;

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
std::ostream& PrettyPrint(std::ostream& out, std::istream& input,
                const FileFormat& matformat)
{
    Givaro::QField<Givaro::Rational> QQ;
    LinBox::MatrixStream<Givaro::QField<Givaro::Rational>> ms( QQ, input );

        // ====================================================
        // Read Matrix
    Givaro::Timer chrono;

    do {
        chrono.start();
        PLinOpt::Matrix A ( ms );
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

#ifdef VERBATIM_PARSING
        std::clog << "# row weights: ";
        for(size_t i(0); i<A.rowdim(); ++i) std::clog << A[i].size() << ' ';
        std::clog << std::endl;
#endif
            // ====================================================
            // Special latex for rationals or one of LinBox formats
        if (matformat == FileFormat::LaTeX) {
            latex_print(out, A) << std::endl;
        } else {

            if ( (matformat == FileFormat::Plain) ||
                 (matformat == FileFormat::HTML)  ||
                 (matformat == FileFormat(6))     ||	// dense
                 (matformat == FileFormat(12))  ) {		// linalg
                PLinOpt::DenseMatrix B(QQ,A.rowdim(),A.coldim());
                PLinOpt::any2dense(B,A);
                if (matformat == FileFormat(6))
                    out << B << std::endl;
                else
                    B.write(out, matformat);
            } else {
                A.write(out, matformat);
            }
        }
            // Stop when no new matrix or empty file
        try { ms.newmatrix(); } catch(LinBox::MatrixStreamError& e) { break; };
    } while(! input.eof());

    return out;
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

                std::clog << "  # is the matrix format among:";
                std::clog << "  Plain = 0, Maple = 1, HTML = 2, LaTeX = 3, SMS = 5, Dense = 6, Matlab = 7, Pretty = 8, OneBased = 10, MatrixMarket = 11, linalg = 12, etc.\n";

                exit(-1);
            }
                // Maple
            else if (args[1] == 'M') { matformats.push_back(FileFormat(1)); }
                // HTML
            else if (args[1] == 'H') { matformats.push_back(FileFormat(2)); }
                // TeX
            else if (args[1] == 'T') { matformats.push_back(FileFormat(3)); }
                // SMS
            else if (args[1] == 'S') { matformats.push_back(FileFormat(5)); }
                // Dense
            else if (args[1] == 'D') { matformats.push_back(FileFormat(6)); }
                // Matlab
            else if (args[1] == 'B') { matformats.push_back(FileFormat(7)); }
                // Pretty
            else if (args[1] == 'P') { matformats.push_back(FileFormat(8)); }
                // OneBased
            else if (args[1] == 'O') { matformats.push_back(FileFormat(10)); }
                // MatrixMarket
            else if (args[1] == 'X') { matformats.push_back(FileFormat(11)); }
                // Linalg
            else if (args[1] == 'L') { matformats.push_back(FileFormat(12)); }
                // Directly by number
            else { matformats.push_back(FileFormat(-atoi(args.c_str()))); }
        }
        else { filenames.push_back(args); }
    }

        // Pretty output by default
    if (matformats.size() == 0)
        matformats.emplace_back(FileFormat::Pretty);

    if (filenames.size() > 0) {
        for(size_t i(0); i<filenames.size(); ++i) {
            std::ifstream input (filenames[i]);
            PrettyPrint(std::cout, input,
                        i<matformats.size() ?
                        matformats[i] :
                        matformats.back()) << std::flush;
        }
    } else {
        PrettyPrint(std::cout, std::cin, matformats.front()) << std::flush;
    }

    return 0;
}
