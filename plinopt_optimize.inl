// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Optimization inline implementations
 ****************************************************************/

#include "plinopt_optimize.h"






template<typename Ring>
std::ostream& printmulorjustdiv(std::ostream& out,
                                const char c, const size_t i,
                                const typename Ring::Element& e,
                                size_t& nbmul, const Ring& F) {
    out << c << i;
    if (notAbsOne(F,e)) {
        ++nbmul;
        out << '*' << e;
    }
    return out;
}

template<>
std::ostream& printmulorjustdiv(std::ostream& out,
                                const char c, const size_t i,
                                const Givaro::Rational& r,
                                size_t& nbmul, const QRat& QQ) {
    out << c << i;
    if (!QQ.isOne(r)) {
        ++nbmul;
        if (Givaro::isOne(r.nume()))
            out << '/' << r.deno();
        else
            out << '*' << r;
    }
    return out;
}



template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::map<T1,T2>& v) {
    out << '{';
    for (const auto& [element, frequency] : v)
        out << element << '=' << frequency << ' ';
    return out << '}';
}

// Build pairs of indices with normalized coefficient (ratio of the two values)
template<typename Container, typename Ring>
std::vector<Etriple<Ring>> listpairs (const Container& c, const Ring& F) {
    std::vector<Etriple<Ring>> v;
    typename Ring::Element tmp;
    for(auto iter=c.begin(); iter != c.end(); ++iter) {
        auto next(iter);
        for(++next; next!= c.end(); ++next) {
            v.emplace_back(iter->first, next->first,
                           F.div(tmp,next->second,iter->second));
        }
    }
    return v;
}

// If cse is present in row, add square of row density to score
template<typename triple>
inline size_t score(const std::vector<std::vector<triple>>& AllPairs,
                    const std::vector<size_t>& Density,
                    const triple& cse) {
    size_t score(0);
    for(size_t k=0; k<Density.size(); ++k) {
        if (std::find(AllPairs[k].begin(), AllPairs[k].end(), cse) != AllPairs[k].end()) {
            score += Density[k]*Density[k];
        }
    }
    return score;
}


template<typename triple, typename _Mat>
Pair<size_t> RemOneCSE(std::ostream& ssout, _Mat& lM, size_t& nbmul,
                       std::vector<triple>& lmultiples, const triple& cse,
                       const std::vector<std::vector<triple>>& AllPairs,
                       const char tev, const char rav) {
    const auto& FF(lM.field());
    size_t savedadds(0), savedmuls(0), lm(lM.coldim()); // The new pair

        // Factor out cse, in all rows
        // adds a column for the new factor
    for(size_t i=0; i<AllPairs.size(); ++i) {
        const auto& rows(AllPairs[i]);
        if (std::find(rows.begin(), rows.end(), cse) != rows.end()) {
            typename _Mat::Element coeff;
            for(auto iter=lM[i].begin(); iter!= lM[i].end(); ++iter) {
                if (iter->first==std::get<0>(cse)) {
                    coeff = iter->second;
                    lM[i].erase(iter); ++savedadds;
                    break;
                }
            }
            for(auto iter=lM[i].begin(); iter!= lM[i].end(); ++iter) {
                if (iter->first==std::get<1>(cse)) {
                    if (notAbsOne(FF,iter->second)) ++savedmuls;
                    lM[i].erase(iter);
                    lM[i].emplace_back(lm,coeff);
                    break;
                }
            }
        }
    }

        // If coefficient was already applied
        //    then reuse the multiplication
        //    else put it in a temporary for future reuse
    auto asgs(abs(std::get<2>(cse)));
    size_t rindex(lm), moremul(0);
    if (!FF.isOne(asgs)) {
        for(const auto& iter: lmultiples) {
            if ((std::get<1>(iter) == std::get<1>(cse)) &&
                (std::get<2>(iter) == asgs)) {
                rindex = std::get<0>(iter);
                break;
            }
        }
        if (rindex == lm) {
            ssout << rav << lm << ":=";
            if (FF.isMOne(asgs)) ssout << '-';
            printmulorjustdiv(ssout, tev, std::get<1>(cse),
                              asgs, moremul, FF) << ';' << std::endl;
            lmultiples.emplace_back(lm,std::get<1>(cse),asgs);
        }
    }

    nbmul += moremul;
    savedmuls -= moremul;

        // Outputs the factor into a temporary variable
    ssout << tev << lm << ":="
          << tev << std::get<0>(cse)
          << (sign(std::get<2>(cse)) >= 0?'+':'-');
    --savedadds;

    if (FF.isOne(asgs))
        ssout << tev << std::get<1>(cse);
    else
        ssout << rav << rindex;
    ssout << ';' << std::endl;

    ++lm;
    lM.resize(lM.rowdim(), lm);

    return Pair<size_t>{savedadds,savedmuls};
}





// Removing one pair
template<typename triple,typename _Mat>
bool OneSub(std::ostream& sout, _Mat& M, std::vector<triple>& multiples,
            size_t& nbmul, const char tev, const char rav) {
// M.write(std::clog << "# BEG OS\n",FileFormat::Pretty) << ';' << std::endl;
    size_t m(M.coldim());
    const auto& FF(M.field());


    std::vector<std::vector<triple>> AllPairs;
    std::vector<size_t> Density;

        // Compute initial density, and all pairs, in a row
    for(auto iter=M.rowBegin(); iter != M.rowEnd(); ++iter) {
        Density.emplace_back(iter->size());
        AllPairs.push_back(listpairs(*iter, FF));
    }

        // Count occurences of each pair in whole matrix
    std::map<triple,size_t> PairMap;
    for(const auto& rows: AllPairs) {
        for (const auto& iter: rows) {
            PairMap[iter]++;
        }
    }

        // Found some pairs
    if (PairMap.size()) {
        size_t maxfrq(0);
        std::vector<triple> MaxCSE;
        triple cse;
            // Find all pairs with maximal frequency
        for (const auto& [element, frequency] : PairMap) {
            if (frequency == maxfrq) {
                MaxCSE.push_back(element);
            }
            if (frequency > maxfrq) {
                maxfrq = frequency;
                cse = element;
                MaxCSE.resize(0); MaxCSE.push_back(element);
            }
        }

            // Factoring will gain something if maximal frequency > 1
        if (maxfrq > 1) {
                // More than one pair with maximal frequency
            if (MaxCSE.size()>1) {
                    // Tie breaking heuristics
#ifdef RANDOM_TIES
                std::shuffle ( MaxCSE.begin(), MaxCSE.end(),
                               std::default_random_engine(Givaro::BaseTimer::seed()) );
                cse = MaxCSE.front();
#else
                size_t maxscore(0);
                for(const auto& element: MaxCSE) {
                    const size_t newscore(score(AllPairs,Density,element));
                    if (newscore == maxscore) {
                            // Tie breaking by multiplier
                        if (FF.isOne(std::get<2>(element)) ||
                            FF.isMOne(std::get<2>(element)) ) {
                            cse = element;
                        }
                    }
                        // Tie breaking by score
                    if (newscore > maxscore) {
                        maxscore = newscore;
                        cse = element;
                    }
                }
#endif
            }

#ifdef VERBATIM_PARSING
            std::clog << "# Found: " << cse << '=' << maxfrq
                      << ',' << score(AllPairs,Density,cse) << std::endl;
            for (const auto& [element, frequency] : PairMap) {
                if ( (frequency == maxfrq) && (element != cse)) {
                    std::clog << "# tied : " << element << '=' << maxfrq
                              << ',' << score(AllPairs,Density,element)
                              << std::endl;
                }
            }
#endif

                // Now factoring out that CSE from the matrix
            const Pair<size_t> savings = RemOneCSE(sout, M, nbmul, multiples,
                                                   cse, AllPairs, tev, rav);

// M.write(std::clog << "# END OS\n",FileFormat::Pretty) << ';' << std::endl;
            return true;
        }
    }
    return false;
}

// Factors out same coefficient in a column
template<typename Iter, typename triple, typename _Mat>
void FactorOutColumns(std::ostream& sout, _Mat& T,
                      std::vector<triple>& multiples, size_t& nbmul,
                      const char tev, const char rav,
                      const size_t j, const Iter& start, const Iter& end) {
    if (start == end) return;
    std::map<typename _Mat::Element, size_t> MapVals;
    for(Iter iter = start; iter != end; ++iter) {
        MapVals[abs(iter->second)]++;
    }

    const auto& FF(T.field());

// T.write(std::clog << "BEG FOC\n",FileFormat::Pretty) << ';'
//                   << "\nMap: " << MapVals << std::endl;

    for (const auto& [element, frequency] : MapVals) {
        size_t m(T.rowdim());

            // Found repeated coefficient
        if ((frequency>1) && (!FF.isOne(element)) ) {
            size_t rindex(m);
                // If coefficient was already applied
                //    then reuse the multiplication
                //    else put it in a temporary for future reuse
             for(const auto& iter: multiples) {
                if ((std::get<1>(iter) == j) &&
                    (std::get<2>(iter) == element)) {
                    rindex = std::get<0>(iter);
                    break;
                }
            }
            if (rindex == m) {
                sout << rav << m << ":=";
                if (FF.isMOne(element)) sout << '-';
                printmulorjustdiv(sout, tev, j,
                                  element, nbmul, FF) << ';' << std::endl;
                multiples.emplace_back(m,j,element);
            }

                // Outputs the coefficient into a temporary variable
            sout << tev << m << ":=";
            sout << rav << rindex << ';' << std::endl;
            ++m;
            T.resize(m, T.coldim());

                // Replace the coefficient in all rows by 1 or -1
            for(size_t k=0; k<frequency; ++k) {
                auto& row(T[j]);
                for(auto iter=row.begin(); iter != row.end(); ++iter) {
                    if (abs(iter->second) == element) {
                        T[m-1].emplace_back(iter->first, (sign(iter->second) >= 0 ? 1 : -1));
                        row.erase(iter);
                        break;
                    }
                }
            }
// T.write(std::clog << "END FOC\n",FileFormat::Pretty) << ';' << std::endl;
        }
    }
}

// Factors out same coefficient in a row
template<typename Iter, typename _Mat>
void FactorOutRows(std::ostream& sout, _Mat& M, size_t& nbadd, const char tev,
                   const size_t i, const Iter& start, const Iter& end) {
    if (start == end) return;
    std::map<typename _Mat::Element, size_t> MapVals;
    for(Iter iter = start; iter != end; ++iter) {
        MapVals[abs(iter->second)]++;
    }

    const auto& FF(M.field());

// M.write(std::clog << "BEG FOR\n",FileFormat::Pretty) << ';'
//                   << "\nMap: " << MapVals << std::endl;

    size_t m(M.coldim());
    for (const auto& [element, frequency] : MapVals) {
            // Found repeated coefficient
        if ((frequency>1) && (!FF.isOne(element)) ) {
            sout << tev << m << ":=";
            ++m;
                // Add a column with coefficient multiplying a new sum
            M.resize(M.rowdim(), m);
            M[i].emplace_back(m-1, element);

                // Remove elements that will be in the new sum
            auto& row(M[i]);
            for(auto iter=row.begin(); iter != row.end(); ++iter) {
                if (abs(iter->second) == element) {
                    if (sign(iter->second) < 0) sout << '-';
                    sout << tev << iter->first;
                    row.erase(iter);
                    break;
                }
            }
                // Precompute the new sum (to be multiplied afterwards)
            for(size_t k=1; k<frequency; ++k) {
                for(auto iter=row.begin(); iter != row.end(); ++iter) {
                    if (abs(iter->second) == element) {
                        ++nbadd;
                        sout << (sign(iter->second) < 0 ? '-' : '+')
                                  << tev << iter->first;
                        row.erase(iter);
                        break;
                    }
                }
            }
            sout << ';' << std::endl;
        }
    }
// M.write(std::clog << "END FOR\n",FileFormat::Pretty) << ';' << std::endl;
}

// Sets new temporaries with the input values
void input2Temps(std::ostream& sout, const size_t N,
                 const char inv, const char tev) {
    // Inputs to temporaries
    for(size_t i=0; i<N; ++i) {
        sout << tev << i << ":="
                  << inv << i << ';' << std::endl;
    }
}
// Sets new temporaries with the input values
template<typename _Mat>
void input2Temps(std::ostream& sout, const size_t N,
                 const char inv, const char tev,
                 const _Mat& trsp) {
    // Inputs to temporaries
    for(size_t i=0; i<N; ++i) {
        if (trsp[i].size()) {
            sout << tev << i << ":="
                      << inv << i << ';' << std::endl;
        } // otherwise variable is not used
    }
}


// Direct program generateur from a matrix
template<typename _Mat, typename triple>
std::ostream& ProgramGen(std::ostream& sout, _Mat& M,
                         std::vector<triple>& multiples,
                         size_t& addcount, size_t& nbmul,
                         const char inv, const char ouv,
                         const char tev, const char rav) {

    const auto& FF(M.field());

        // Factoring multiplier by colums
    _Mat T(FF, M.coldim(), M.rowdim()); Transpose(T, M);
    for(size_t j=0; j<M.coldim(); ++j) {
        FactorOutColumns(sout, T, multiples, nbmul, tev, rav,
                         j, T[j].begin(), T[j].end());
    }
    Transpose(M,T);

        // Factoring multiplier by rows
    for(size_t i=0; i<M.rowdim(); ++i) {
        FactorOutRows(sout, M, addcount, tev, i, M[i].begin(), M[i].end());
    }

        // Computing remaining (simple) linear combinations
    for(size_t i=0; i<M.rowdim(); ++i) {
        const auto& row(M[i]);
        sout << ouv << i << ":=";
        if (row.size()>0) {

            if ( (sign(row.begin()->second) < 0) || FF.isMOne(row.begin()->second) ) sout << '-';
            auto arbs(abs(row.begin()->second));

                // If already multiplied, reuse it
            size_t rindex(M.coldim());
            for(const auto& miter: multiples) {
                if ((std::get<1>(miter) == row.begin()->first) &&
                    (std::get<2>(miter) == arbs)) {
                    rindex = std::get<0>(miter);
                    break;
                }
            }
            if (rindex != M.coldim()) {
                sout << rav << rindex;
            } else {
                printmulorjustdiv(sout, tev, row.begin()->first,
                                  arbs, nbmul, FF);
            }

                // For all the monomials of the linear combination
            auto iter(row.begin());
            for(++iter; iter!= row.end(); ++iter) {
                ++addcount;
                sout << ( ( (sign(iter->second) <0) || FF.isMOne(iter->second) ) ? '-' : '+');
                auto ais(abs(iter->second));

                    // If already multiplied, reuse it
                size_t rindex(M.coldim());
                for(const auto& miter: multiples) {
                    if ((std::get<1>(miter) == iter->first) &&
                        (std::get<2>(miter) == ais)) {
                        rindex = std::get<0>(miter);
                        break;
                    }
                }
                if (rindex != M.coldim()) {
                    sout << rav << rindex;
                } else {
                        // otherwise next function will sout << '-'
                    printmulorjustdiv(sout, tev, iter->first,
                                      ais, nbmul, FF);
                }
            }
        } else {
            sout << '0';
        }
        sout << ';' << std::endl;
    }


// std::clog << std::string(30,'#') << std::endl;
// M.write(std::clog << "M:=",FileFormat::Maple) << ';' << std::endl;
    return sout;
}


// Global random optimization function (pairs and factors)
template<typename _Mat>
Pair<size_t> Optimizer(std::ostream& sout, _Mat& M,
                       const char inv, const char ouv,
                       const char tev, const char rav) {

    using triple=std::tuple<size_t, size_t, typename _Mat::Element>;
    size_t nbadd(0), nbmul(0);


        // Factoring sums
    std::vector<triple> multiples;
    for( ; OneSub(sout, M, multiples, nbmul, tev, rav) ; ++nbadd) { }

        // No more useful CSE, thus
        // generate the rest of the program from the new M
    ProgramGen(sout, M, multiples, nbadd, nbmul, inv, ouv, tev, rav);

    return Pair<size_t>(nbadd,nbmul);
}



	// Postcondition _Matrix A is nullified
template<typename _Mat>
Pair<size_t> nullspacedecomp(std::ostream& sout, _Mat& x, _Mat& A) {
    const auto& FF(A.field());
    typename _Mat::Element Det;
    size_t Rank;
    size_t Ni(A.rowdim()),Nj(A.coldim());

    _Mat FreePart(FF, A.coldim(), A.rowdim()); Transpose(FreePart,A);

        // ============================================
        // Find the rows to start the kernel elimination
    std::vector<size_t> Q(Nj);
    std::iota(Q.begin(), Q.end(), 0); // Q will be this permutation


#ifdef RANDOM_TIES
        // Randomly swap initial rows of FreePart
    std::vector<size_t> l(Nj);
    std::iota(l.begin(), l.end(), 0); // Q will be this permutation
    std::shuffle ( l.begin(), l.end(),
                   std::default_random_engine(Givaro::BaseTimer::seed()));

    for(size_t i=1; i<Nj; ++i) {
        if (i != l[i]) {
            std::swap(FreePart[i],FreePart[l[i]]);
            std::swap(Q[i], Q[l[i]]);
        }
    }
#else
        // Find sparsest initial rows of FreePart
    for(size_t i=1; i<Nj; ++i) {
        for(size_t j=0; j<i; ++j) {
            if (FreePart[i].size()<FreePart[j].size()) {
                std::swap(FreePart[i],FreePart[j]);
                std::swap(Q[i], Q[j]);
                --i;
                break;
            }
        }
    }
#endif

    Transpose(A, FreePart); // Apply the Q permutation on A

        // ============================================
        // Now find subset of independent rows
    LinBox::Permutation<typename _Mat::Field> P(FF,(int)Nj);
    LinBox::GaussDomain<typename _Mat::Field> GD(FF);

        // LUP decomposition
    GD.InPlaceLinearPivoting(Rank, Det, A, P, Ni, Nj );

        // Removing zero rows
    for(size_t i=0; i< Ni; ++i) {
        if (A[i].size() == 0) {
            size_t j(i);
            if (nextnonzero(j,Ni,A)) {
                A[i] = A[j];
                A[j].resize(0);
            }
            else {
                break;
            }
        }
    }

        // Compute NullSpace from U
    size_t nullity = A.coldim()-Rank;
    x.resize(x.rowdim(),nullity);
#ifdef VERBATIM_PARSING
    std::clog << "# NullSpace dimensions:" << x.rowdim() << 'x' << x.coldim() << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif
    if ( (Rank != 0) && (nullity != 0) ) {
            // compute U2T s.t. U = [ U1 | -U2T^T ]
        _Mat U2T(FF,nullity,Rank);

        for(auto uit=A.IndexedBegin(); uit != A.IndexedEnd(); ++uit) {
            if (uit.colIndex() >= Rank)
                U2T.setEntry(uit.colIndex()-Rank,uit.rowIndex(),uit.value());
        }
        for(auto u2it=U2T.Begin(); u2it != U2T.End(); ++u2it)
            FF.negin(*u2it);

            // Compute the basis vector by vector
        typedef LinBox::Sparse_Vector<typename _Mat::Element> SparseVect;
        for(size_t i=0; i<nullity; ++i) {
            SparseVect W1Ti;
                // Solve for upper part of basis
            upperTriangularSparseSolve(W1Ti, Rank, A, U2T[i]);

            for(size_t j=0; j<W1Ti.size(); ++j) {
                    // P.applyTranspose(x[i],W1T[i]);
                    // Transposein(x)
                x.setEntry( (size_t)P.getStorage()[ (size_t)W1Ti[(size_t)j].first ], i, W1Ti[j].second );
            }
        }


            // ============================================
            // Remove dependent rows from FreePart
        for(size_t i=0; i<nullity; ++i) {
            FreePart[P.getStorage()[ Rank+i ]].resize(0);
        }

            // ============================================
            // Unapply Q
        std::vector<size_t> cQ(Q);
        for(size_t i=0; i<Nj; ++i) {
            while(cQ[i] != i) {
                std::swap(FreePart[i], FreePart[ cQ[i] ]);
                std::swap(x[i], x[ cQ[i] ]);
                std::swap(cQ[i], cQ[ cQ[i] ]);
           }
        }



#ifdef VERBATIM_PARSING
        FreePart.write(std::clog << "# FreePart:", FileFormat::Maple)
                                 << std::endl;
        std::clog << std::string(30,'#') << std::endl;
#endif


            // ============================================
            // Optimize the set of independent rows
        input2Temps(sout, FreePart.coldim(), 'i', 't');
            // o <-- Free . i
        auto Fops( Optimizer(sout, FreePart, 'i', 'o', 't', 'r') );


            // ============================================
            // Optimize the dependent rows (transposed nullspace)
            // [ Free^T | Dep^T ] . [ x^T | I ]^T = 0
            // So Dep . i = (- x^T) Free . i = (- x^T) o
         _Mat Tx(x.field(), x.coldim(), x.rowdim()); NegTranspose(Tx, x);

#ifdef VERBATIM_PARSING
        std::clog << std::string(30,'#') << std::endl;
        Tx.write(std::clog << "# Dependent:", FileFormat::Maple) << std::endl;
#endif

        input2Temps(sout, Tx.coldim(), 'o', 'v', x);
            // x <-- Tx . o = (- x^T) o
        auto Kops( Optimizer(sout, Tx, 'o', 'x', 'v', 'g') );

            // ============================================
            // Recover final output of NullSpace
            // Applying both permutations, P then Q
            // back into 'o' variables
        for(size_t i=0; i<nullity; ++i) {
            sout << 'o' << Q[ P.getStorage()[ Rank+i ] ] << ":="
                 << 'x' << i << ';' << std::endl;
        }

            // ============================================
            // Total number of operations
        Fops.first += Kops.first;
        Fops.second += Kops.second;

        return Fops;
    }
    return Pair<size_t>{-1,-1}; // +infty,+infty
}



#include <deque>

// Recusive search for the best cse
template<typename triple,typename _Mat>
bool RecSub(std::deque<std::string>& out, _Mat& Mat,
            std::vector<triple>& multiples,
            size_t& nbadd, size_t& nbmul, const size_t lvl,
            const char tev, const char rav) {

    size_t m(Mat.coldim());
    const auto& FF(Mat.field());

    std::vector<std::vector<triple>> AllPairs;
    std::vector<size_t> Density;
        // Compute initial density, and all pairs, in a row
    for(auto iter=Mat.rowBegin(); iter != Mat.rowEnd(); ++iter) {
        AllPairs.push_back(listpairs(*iter, FF));
    }


        // Count occurences of each pair in whole matrix
    std::map<triple,size_t> PairMap;
    for(const auto& rows: AllPairs) {
        for (const auto& iter: rows) {
            PairMap[iter]++;
        }
    }
        // Found some pairs
    if (PairMap.size()) {
        size_t maxfrq(0);
            // Find all pairs with maximal frequency
        for (const auto& [element, frequency] : PairMap) {
            if (frequency > maxfrq) {
                maxfrq = frequency;
            }
        }
            // Factoring will gain something if maximal frequency > 1
        if (maxfrq <= 1) {
            return false;
        }
    }




    triple cse{0,0,0};
    for(const auto& rows: AllPairs) {
        if (rows.size() > 0) cse = rows.front();
    }

    size_t bestadds(nbadd), bestmuls(nbmul);
    _Mat bestM(FF,Mat.rowdim(),Mat.coldim()); sparse2sparse(bestM, Mat);
    std::vector<triple> bestmultiples(multiples);
    std::deque<std::string> bestdout;
//     std::clog << std::string(lvl,'#') << " lvl(" << lvl << "), " << nbadd << '|' << nbmul << ", allpairs: ";
//     size_t subexps(0); for(const auto& rows: AllPairs) subexps+=rows.size();
//    std::clog << subexps << std::endl;


    for(const auto& rows: AllPairs) {
        for (const auto& cse: rows) {
            if (PairMap[cse] > 1) {
                _Mat lM(FF,Mat.rowdim(),Mat.coldim()); sparse2sparse(lM, Mat);
                std::vector<triple> lmultiples(multiples);
                std::ostringstream ssout;
                size_t moremul(0);

                    // Factoring out that CSE from the matrix
                const Pair<size_t> savings = RemOneCSE(ssout, lM, moremul,
                                                       lmultiples, cse,
                                                       AllPairs, tev, rav);

                size_t ladditions(nbadd-savings.first);
                size_t lmuls(nbmul-savings.second);

                std::deque<std::string> sdout;
                RecSub(sdout, lM, lmultiples, ladditions, lmuls, lvl+1, tev, rav);

                if ( (ladditions < bestadds) ||
                     ( (ladditions == bestadds) && (lmuls < bestmuls) ) ) {
                    bestadds = ladditions;
                    bestmuls = lmuls;
                    sparse2sparse(bestM, lM);
                    bestmultiples = lmultiples;
                    bestdout = sdout;
                    bestdout.push_front(ssout.str());
                }
            }
        }
    }

    nbadd = bestadds;
    nbmul = bestmuls;
    multiples = bestmultiples;
    sparse2sparse(Mat, bestM);
    out.insert(out.end(),bestdout.begin(), bestdout.end());

#ifdef VERBATIM_PARSING
    std::clog << std::string(lvl,'#') << "# lvl(" << lvl << "), best: " << bestadds << '|' << bestmuls << std::endl;
#endif

    return true;
}



// Global exhaustive optimization function (pairs and factors)
template<typename _Mat>
Pair<size_t> RecOptimizer(std::ostream& sout, _Mat& M,
                          const char inv, const char ouv,
                          const char tev, const char rav) {
    using triple=std::tuple<size_t, size_t, typename _Mat::Element>;
    const auto& FF(M.field());


    size_t nbadd(0), nbmul(0);
    for(auto iter=M.rowBegin(); iter != M.rowEnd(); ++iter)
        nbadd += std::max((int)iter->size()-1,0);
    for(auto it = M.IndexedBegin(); it != M.IndexedEnd(); ++it)
        if (notAbsOne(FF,it.value())) ++nbmul;

        // Factoring sums
    std::deque<std::string> dout;
    std::vector<triple> multiples;
    size_t lmul(0);
    RecSub(dout, M, multiples, nbadd, nbmul, 0, tev, rav);

    size_t alreadymuls(nbmul);
    for(auto it = M.IndexedBegin(); it != M.IndexedEnd(); ++it)
        if (notAbsOne(FF,it.value())) --alreadymuls;

    for(auto iter = dout.rbegin(); iter != dout.rend(); ++iter) {
        sout << *iter;
    }

    size_t addcount(0);


    ProgramGen(sout, M, multiples, addcount, alreadymuls, inv, ouv, tev, rav);

    return Pair<size_t>(nbadd,alreadymuls);
}
