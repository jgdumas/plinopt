// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, inline implementations
 ****************************************************************/

#include "plinopt_library.h"

Matrix& Transpose(Matrix& T, const Matrix& A) {
    T.resize(0,0);
    T.resize(A.coldim(), A.rowdim());
    for(auto it = A.IndexedBegin(); it != A.IndexedEnd(); ++it)
        T.setEntry(it.colIndex(),it.rowIndex(), it.value());
    return T;
}

Matrix& NegTranspose(Matrix& T, const Matrix& A) {
    T.resize(0,0);
    T.resize(A.coldim(), A.rowdim());
    typename Matrix::Element tmp; T.field().init(tmp);
    for(auto it = A.IndexedBegin(); it != A.IndexedEnd(); ++it)
        T.setEntry(it.colIndex(),it.rowIndex(), T.field().neg(tmp,it.value()));
    return T;
}
