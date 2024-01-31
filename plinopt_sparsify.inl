// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Sparsification  inline implementations
 ****************************************************************/

#include "plinopt_sparsify.h"


	// Find most frequent Rational between begin and end
    // element is suppoded to be the second in a pair
template <typename Fwd>
Givaro::Rational most_frequent_element(const Fwd& begin, const Fwd& end)
{
    std::map<Givaro::Rational, int> count;

    for (auto it = begin; it != end; ++it) ++count[it->second];

    return std::max_element(
        count.begin(), count.end(),
        [](const auto& a, const auto& b) { return a.second < b.second;})->first;
}

	// Replace row i of A, by v
template<typename Vector>
Matrix& setRow(Matrix& A, size_t i, const Vector& v, const QRat& QQ) {
    A[i].resize(0);
    for(size_t j=0; j<v.size(); ++j) {
        if (! QQ.isZero(v[j])) A.setEntry(i,j,v[j]);
    }
    return A;
}

	// copy dense matrix M into sparse matrix A
Matrix& dense2sparse(Matrix& A, const DenseMatrix& M, const QRat& QQ) {
    A.resize(M.rowdim(), M.coldim());
    for(size_t i=0; i<A.rowdim(); ++i) {
        setRow(A,i,M[i],QQ);
    }
    return A;
}

	// copy sparse matrix B into sparse matrix A
Matrix& sparse2sparse(Matrix& A, const Matrix& B) {
    A.resize(B.rowdim(), B.coldim());
    auto rowB=B.rowBegin();
    for(auto rowA=A.rowBegin(); rowA!=A.rowEnd();++rowA, ++rowB) {
        *rowA = *rowB;
    }

    return A;
}

	// v is augmented by i, -i, 1/i and -1/i
template<typename Vector>
void augment(Vector& v, const Givaro::Integer& i, const QRat& QQ) {
    const Givaro::Rational r(i);
    if (std::find(v.begin(), v.end(), r) == v.end()) {
        v.push_back(r);
        v.push_back(-r);
        Givaro::Rational tmp;
        QQ.inv(tmp, r);
        v.push_back(tmp);
        QQ.negin(tmp);
        v.push_back(tmp);
    }
}

	// Computes the rank of A
size_t& rank(size_t& r, const Matrix& A) {
    static QRat QQ;
    static LinBox::GaussDomain<QRat> GD(QQ);
    Matrix copyA(QQ); sparse2sparse(copyA, A);
    return GD.rankInPlace(r, copyA);
}


    // Prints out density profile of M
    // returns total density
std::ostream& densityProfile(std::ostream& out, size_t& ss, const Matrix& M) {
    ss = 0;
    out << "# Density profile: ";
    for(auto it=M.rowBegin();it!=M.rowEnd();++it) {
        ss += it->size();
        out << it->size() << ' ';
    }
    return out << ':' << ss;
}

// ============================================================
// Sparsifying a Matrix M
//   uses a limited number of coefficients for the linear comb.
//   (that is at most: maxnumcoeff)
Matrix& Sparsifier(Matrix& CoB, Matrix& M, const size_t maxnumcoeff) {

        // ========================================
        // Read Matrix of Linear Transformation
    static QRat QQ;
    Matrix TCoB(QQ,4,4);

        // ========================================
        // Try 0, 1, -1 and some coefficients in M
    std::vector<Givaro::Rational> Coeffs{0, 1, -1};
    for(auto row=M.rowBegin(); row != M.rowEnd(); ++row) {
        for(auto it=row->begin(); it != row->end(); ++it) {
            augment(Coeffs, it->second.nume(), QQ);
            augment(Coeffs, it->second.deno(), QQ);
         }
    }

        // ========================================
        // reduce to at most maxnumcoeff
    if (Coeffs.size()>maxnumcoeff) Coeffs.resize(maxnumcoeff);


        // ========================================
        // Try each column of CoB, one at a time
    for(size_t num=0; num<TCoB.rowdim(); ++num) {

        QArray v(M.rowdim());
        Matrix A(QQ); sparse2sparse(A, TCoB);
        int hw(-1); // Best Hamming weight so far


            // ========================================
            // Ech column has 4 coefficients
        for(size_t i=0; i<Coeffs.size(); ++i) {
        for(size_t j=0; j<Coeffs.size(); ++j) {
        for(size_t k=0; k<Coeffs.size(); ++k) {
        for(size_t l=0; l<Coeffs.size(); ++l) {
                // Try linear combination
            QArray w{Coeffs[i],Coeffs[j],Coeffs[k],Coeffs[l]};
                // Check independency
            setRow(A,num,w,QQ);
            size_t r; rank(r, A);

            if (r>num) {
                M.apply(v,w);
                    // Count number of produced zeroes
                int hwl = std::count_if(v.begin(), v.end(),
                                        [](const auto& e) { return isZero(e);} );
                    // Find the best one so far
               if (hwl > hw) {
#ifdef VERBATIM_PARSING
                    std::clog << "Found: " << w << ' ' << hwl << std::endl;
#endif
                    hw = hwl;
                        // Register best linear combination so far
                    setRow(TCoB,num,w,QQ);
                }
            }
        }}}}
    }

        // ========================================
        // Supplementary change of basis found
    Matrix C(QQ,CoB.rowdim(), CoB.coldim());
    Transpose(C,TCoB);

        // ========================================
        // Apply it:
        //   - to the sparsified matrix
        //   - to previous change of Basis
    static LinBox::MatrixDomain<QRat> BMD(QQ);
    DenseMatrix R(M), S(CoB);
    BMD.mul(R,M,C); // Direct matrix multiplication
    BMD.mul(S,CoB,C);


    dense2sparse(M,R,QQ);
    dense2sparse(CoB,S,QQ);

    return CoB;
}


// ============================================================
// Factoring out some coefficients:
//    from most frequent in a column of M, towards a row of CoB
Matrix& FactorDiagonals(Matrix& CoB, Matrix& M) {
    static QRat QQ;
    Matrix T(QQ), C(QQ);
    Transpose(T, M);
    Transpose(C, CoB);

    for(size_t i=0; i<T.rowdim(); ++i) {
        Givaro::Rational r(most_frequent_element(T[i].begin(), T[i].end() ) );
        for(auto& iter : T[i]) iter.second /= r; // scale T
        for(auto& iter : C[i]) iter.second /= r; // scale C
    }

    Transpose(M,T);
    Transpose(CoB, C);

    return CoB;
}



// ============================================================
// Consistency check of M == R.C
std::ostream& consistency(std::ostream& out, const Matrix& M,
                          const DenseMatrix& R, const Matrix& C) {
    static QRat QQ;
    static LinBox::MatrixDomain<QRat> BMD(QQ);
    DenseMatrix A(R);

    BMD.mul(A,R,C);
    BMD.subin(A,M);

    if (BMD.isZero (A))
        out <<"# \033[1;32mOK: consistent factorization!\033[0m";
    else{
        std::cerr << "# \033[1;31m****** ERROR inconsistency ******\033[0m"
                  << std::endl;
        M.write(out,FileFormat::Pretty) << std::endl;
        out << " != " << std::endl;
        R.write(out,FileFormat::Pretty) << std::endl;
        out << " * " << std::endl;
        C.write(out,FileFormat::Pretty) << std::endl;
        out << std::string(20,'#') << std::endl;
    }

    return out;
}
