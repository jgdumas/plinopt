// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Optimization inline implementations
 ****************************************************************/

#include "plinopt_optimize.h"
#include "plinopt_sparsify.h"



namespace PLinOpt {
// ===============================================================

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
    typename Ring::Element tmp; F.init(tmp);
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
        if (std::find(AllPairs[k].begin(), AllPairs[k].end(), cse)
            != AllPairs[k].end()) {
            score += Density[k]*Density[k];
        }
    }
    return score;
}

// Removes cse, updating lM (and AllPairs/PairMap if updateAPM)
template<typename triple, typename _Mat>
Pair<long> RemOneCSE(std::ostream& ssout, _Mat& lM, size_t& nbmul,
                     std::vector<triple>& lmultiples, const triple& cse,
                     std::vector<std::vector<triple>>& AllPairs,
                     std::map<triple,size_t>& PairMap,
                     const char tev, const char rav, const bool updateAPM) {
    const auto& FF(lM.field());
    size_t lm(lM.coldim());
    long savedadds(0), savedmuls(0);

        // Looking for the most number of ones in either columns of cse
    _Mat lT(FF, lM.coldim(), lM.rowdim()); Transpose(lT, lM);
    size_t count0(0), count1(0);
    for(const auto& iter: lT[std::get<0>(cse)]) {
        if (isAbsOne(FF,iter.second)) ++count0;
    }
    for(const auto& iter: lT[std::get<1>(cse)]) {
        if (isAbsOne(FF,iter.second)) ++count1;
    }

    triple lcse;
    if (count0<count1) { // More ones in std::get<1>(cse);
        std::get<0>(lcse) = std::get<1>(cse);
        std::get<1>(lcse) = std::get<0>(cse);
        FF.inv( std::get<2>(lcse), std::get<2>(cse) );
    } else {             // More ones in std::get<0>(cse);
        std::get<0>(lcse) = std::get<0>(cse);
        std::get<1>(lcse) = std::get<1>(cse);
        FF.assign( std::get<2>(lcse), std::get<2>(cse) );
    }

        // Factor out cse, in all rows
        // adds a column for the new factor
    for(size_t i=0; i<AllPairs.size(); ++i) {
        const auto& rows(AllPairs[i]);
        if (std::find(rows.begin(), rows.end(), cse) != rows.end()) {
            typename _Mat::Element coeff; FF.init(coeff);
            for(auto iter=lM[i].begin(); iter!= lM[i].end(); ++iter) {
                if (iter->first==std::get<0>(lcse)) {
                    coeff = iter->second;
                    lM[i].erase(iter); ++savedadds;
                    break;
                }
            }
            for(auto iter=lM[i].begin(); iter!= lM[i].end(); ++iter) {
                if (iter->first==std::get<1>(lcse)) {
                    if (notAbsOne(FF,iter->second)) ++savedmuls;
                    lM[i].erase(iter);
                    lM[i].emplace_back(lm,coeff);
                    break;
                }
            }

            if (updateAPM) {
                    // Update AllPairs[i] and PairMap
                    // Decrease PairMap by old AP[i]
                for (const auto& iter: AllPairs[i]) {
                    PairMap[iter]--;
                    if (PairMap[iter] == 0) PairMap.erase(iter);
                }

                    // Preserve pairs not involving the change
                std::vector<triple> newrow;
                for (const auto& iter: AllPairs[i]) {
                    if ( (std::get<0>(iter) != std::get<0>(lcse)) &&
                         (std::get<1>(iter) != std::get<0>(lcse)) &&
                         (std::get<0>(iter) != std::get<1>(lcse)) &&
                         (std::get<1>(iter) != std::get<1>(lcse))
                         ) {
                        newrow.push_back(iter);
                    }
                }
                     // Add all pairs with new (last) element
                typename _Mat::Element tmp; FF.init(tmp);
                auto endlMi(lM[i].end()); --endlMi;
                for(auto iter=lM[i].begin(); iter != endlMi; ++iter) {
                    newrow.emplace_back(iter->first, endlMi->first,
                                        FF.div(tmp, endlMi->second, iter->second));
                }

                     // Increase PairMap by new AP[i]
                for (const auto& iter: newrow) {
                     PairMap[iter]++;
                }

                    // Update AllPairs[i]
                AllPairs[i].assign(newrow.begin(),newrow.end());
            }
        }
    }

        // If coefficient was already applied
        //    then reuse the multiplication
        //    else put it in a temporary for future reuse
    auto asgs(Fabs(FF,std::get<2>(lcse)));
    size_t rindex(lm), moremul(0);
    if (notAbsOne(FF,asgs)) {
        for(const auto& iter: lmultiples) {
            if ((std::get<1>(iter) == std::get<1>(lcse)) &&
                (std::get<2>(iter) == asgs)) {
                rindex = std::get<0>(iter);
                break;
            }
        }
        if (rindex == lm) {
            ssout << rav << lm << ":=";
            printmulorjustdiv(ssout, tev, std::get<1>(lcse),
                              asgs, moremul, FF) << ';' << std::endl;
            lmultiples.emplace_back(lm,std::get<1>(lcse),asgs);
        }
    }

    nbmul += moremul;
    savedmuls -= moremul;

        // Outputs the factor into a temporary variable
    ssout << tev << lm << ":=" << tev << std::get<0>(lcse);
    if ( (FF.isMOne(asgs)) || (Fsign(FF,std::get<2>(lcse))<0) ) {
        ssout << '-';
    } else {
        ssout << '+';
    }
    if (isAbsOne(FF,asgs)) {
        ssout << tev << std::get<1>(lcse);
    } else {
        ssout << rav << rindex;
    }
    ssout << ';' << std::endl;

    --savedadds;

    ++lm;
    lM.resize(lM.rowdim(), lm);

    return Pair<long>{savedadds,savedmuls};
}

template<typename triple, typename _Mat>
Pair<long> RemOneCSE(std::ostream& ssout, _Mat& lM, size_t& nbmul,
                     std::vector<triple>& lmultiples, const triple& cse,
                     std::vector<std::vector<triple>>& AllPairs,
                     const char tev, const char rav) {
    std::map<triple,size_t> PairMap;
    return RemOneCSE(ssout, lM, nbmul, lmultiples, cse, AllPairs,
                     PairMap, tev, rav, false);
}


// Removing one pair (CSE)
template<typename triple,typename _Mat>
bool OneSub(std::ostream& sout, _Mat& M, std::vector<triple>& multiples,
            size_t& nbadd, size_t& nbmul, const char tev, const char rav) {
    const auto& FF(M.field());

        // Compute all pairs, in a row
    std::vector<std::vector<triple>> AllPairs;
    for(auto iter=M.rowBegin(); iter != M.rowEnd(); ++iter) {
        AllPairs.push_back(listpairs(*iter, FF));
    }

        // Count occurences of each pair in whole matrix
    std::map<triple,size_t> PairMap;
    for(const auto& rows: AllPairs) {
        for (const auto& iter: rows) {
            PairMap[iter]++;
        }
    }

#ifdef VERBATIM_PARSING
    size_t apcount(0); for(const auto& rows: AllPairs) apcount += rows.size();
    std::clog << "# Pairs: " << PairMap.size() << '|' << apcount << '|'
              << double(apcount)/double(AllPairs.size()) << std::endl;
#endif

    if (PairMap.size()==0) return false;

    bool goodfreq(false);

    while(PairMap.size() > 0) {
            // Found some pairs
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

        if (maxfrq <=1) return goodfreq;
            // Factoring will gain something if maximal frequency > 1
        goodfreq = true;

            // More than one pair with maximal frequency
        if (MaxCSE.size()>1) {
                // Tie breaking heuristics
#if defined(RANDOM_TIES) && !defined(DENSITY_OPTIMIZATION)
            static thread_local Givaro::GivRandom
                generator(Givaro::BaseTimer::seed());
            cse = MaxCSE[ generator() % MaxCSE.size() ];
#else
                // Compute density in a row
            std::vector<size_t> Density;
            for(auto iter=M.rowBegin(); iter != M.rowEnd(); ++iter) {
                Density.emplace_back(iter->size());
            }

            size_t maxscore(0);
            for(const auto& element: MaxCSE) {
                const size_t newscore(score(AllPairs,Density,element));
                if (newscore == maxscore) {
                        // Tie breaking by multiplier
                    if (isAbsOne(FF,std::get<2>(element))) {
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

            // Now factoring out that CSE from the matrix
        ++nbadd;
        auto savings(RemOneCSE(sout, M, nbmul, multiples, cse,
                               AllPairs, PairMap, tev, rav, true));

#ifdef VERBATIM_PARSING
#  if VERBATIM_PARSING >= 2
        printEtriple(std::clog << "# Found: ", FF, cse) << '=' << maxfrq
                               << ", Saved: " << savings.first << "+|"
                               << savings.second << 'x' << std::endl;
#  endif
#  if VERBATIM_PARSING >= 3
        for (const auto& [element, frequency] : PairMap) {
            if ( (frequency == maxfrq) && (element != cse)) {
                printEtriple(std::clog << "# tied : ", FF, element)
                                       << '=' << maxfrq << ','
                                       << score(AllPairs,Density,element)
                                       << std::endl;
            }
        }
#  endif
#endif
    }
    return true;
}

// Factors out same coefficient in a column
template<typename Iter, typename triple, typename _Mat>
void FactorOutColumns(std::ostream& sout, _Mat& T,
                      std::vector<triple>& multiples, size_t& nbmul,
                      const char tev, const char rav,
                      const size_t j, const Iter& start, const Iter& end) {
    if (start == end) return;
    const auto& FF(T.field());

    std::map<typename _Mat::Element, size_t> MapVals;
    for(Iter iter = start; iter != end; ++iter) {
        MapVals[Fabs(FF,iter->second)]++;
    }


    for (const auto& [element, frequency] : MapVals) {
        size_t m(T.rowdim());

            // Found repeated coefficient
        if ( (frequency>1) && (notAbsOne(FF,element)) ) {
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
                printmulorjustdiv(sout, tev, j,
                                  element, nbmul, FF) << ';' << std::endl;
                multiples.emplace_back(m,j,element);
            }
                // Outputs the coefficient into a temporary variable
            sout << tev << m << ":=";
            sout << rav << rindex << ';' << std::endl;
            T.resize(++m, T.coldim());

                // Replace the coefficient in all rows by 1 or -1
            for(size_t k=0; k<frequency; ++k) {
                auto& row(T[j]);
                for(auto iter=row.begin(); iter != row.end(); ++iter) {
                    if (Fabs(FF,iter->second) == element) {
                        T[m-1].emplace_back(iter->first, (Fsign(FF,iter->second) >= 0 ? FF.one : FF.mOne));
                        row.erase(iter);
                        break;
                    }
                }
            }
        }
    }
}

// Factors out same coefficient in a row
template<typename Iter, typename _Mat>
void FactorOutRows(std::ostream& sout, _Mat& M, size_t& nbadd, const char tev,
                   const size_t i, const Iter& start, const Iter& end) {
    if (start == end) return;
    const auto& FF(M.field());

    std::map<typename _Mat::Element, size_t> MapVals;
    for(Iter iter = start; iter != end; ++iter) {
        MapVals[Fabs(FF,iter->second)]++;
    }

    size_t m(M.coldim());
    for (const auto& [element, frequency] : MapVals) {
            // Found repeated coefficient
        if ((frequency>1) && (notAbsOne(FF,element)) ) {
            sout << tev << m << ":=";
            ++m;
                // Add a column with coefficient multiplying a new sum
            M.resize(M.rowdim(), m);
            M[i].emplace_back(m-1, element);

                // Remove elements that will be in the new sum
            auto& row(M[i]);
            for(auto iter=row.begin(); iter != row.end(); ++iter) {
                if (Fabs(FF,iter->second) == element) {
                    if (Fsign(FF,iter->second) < 0) sout << '-';
                    sout << tev << iter->first;
                    row.erase(iter);
                    break;
                }
            }
                // Precompute the new sum (to be multiplied afterwards)
            for(size_t k=1; k<frequency; ++k) {
                for(auto iter=row.begin(); iter != row.end(); ++iter) {
                    if (Fabs(FF,iter->second) == element) {
                        ++nbadd;
                        sout << (Fsign(FF,iter->second) < 0 ? '-' : '+')
                                  << tev << iter->first;
                        row.erase(iter);
                        break;
                    }
                }
            }
            sout << ';' << std::endl;
        }
    }
}

// Factors out triangles:
// < ab | b > replaced by < 0 | b | b > then < 0 | 0 | 0 | 1>
// < a  | . >             < 0 | . | 1 >      < 0 | . | 1 | 0>
// With only 2 multiplications, by a, then by b, instead of 3
template<typename triple, typename _Mat>
bool Triangle(std::ostream& sout, _Mat& M, _Mat& T,
              std::vector<triple>& multiples, size_t& nbadd, size_t& nbmul,
              const char tev, const char rav, const size_t j) {

    if (T[j].size() == 0) return false;
    const auto& FF(T.field());

    bool found = false;
    bool over(true);

    do { over = true;

    for(auto iter = T[j].begin(); iter != T[j].end();
        ++iter) { if (notAbsOne(FF, iter->second)) {
        for(auto next = T[j].begin(); next != T[j].end();
            ++next) { if ((next != iter) && notAbsOne(FF,next->second)) {
            const size_t i(next->first);
            typename _Mat::Element quot; FF.init(quot);
            const auto& rowi(M[i]);
            FF.div(quot, next->second, iter->second); // nnz should be inv.
            for(auto third=rowi.begin(); third != rowi.end();
                ++third) { if ( (third->first != j)
                                && (notAbsOne(FF,third->second)) ) {
                typename _Mat::Element coeff; FF.init(coeff);
                FF.div(coeff, quot, third->second);
                if (isAbsOne(FF,coeff)) {
                    size_t m(T.rowdim());
                    found = true; // Triangle found !!!
                    over = false; // will have to loop again

                        // First, record one multiplication by a
					sout << tev << m << ":=";
                    auto ais(Fabs(FF,iter->second));
                    if ( (Fsign(FF,iter->second) <0)
                         || FF.isMOne(iter->second) ) sout << '-';
                    printmulorjustdiv(sout, tev, j, ais,
                                      nbmul, FF) << ';' << std::endl;
                    multiples.emplace_back(m,j,iter->second);

                        // Second, in j-th column, divide both elements by a
#ifdef VERBATIM_PARSING
#  if VERBATIM_PARSING >= 2
                    size_t dummy(0);
                    printmulorjustdiv(std::clog << "# Found Triangle [",
                                      tev, j, ais, dummy, FF);
                    std::clog << "]: (" << iter->first << ',' << j << ')' << '('
                              << i << ',' << j << ')' << '('
                              << i << ',' << third->first << ')' << std::endl;
#  endif
#endif

                        //   second.i) add column with 1 and quot
                    T.resize(++m, T.coldim());
                    T[m-1].emplace_back(iter->first, FF.one);
                    T[m-1].emplace_back(next->first, quot);

                        //   second.ii) remove iter & next from T[j]
                    const auto& diter(*iter), & dnext(*next);
                    auto isIorN {[diter,dnext](
                        const typename _Mat::Row::value_type& p)
                        { return ((p == diter) || (p == dnext));}};
                    T[j].erase(std::remove_if(T[j].begin(), T[j].end(), isIorN),
                               T[j].end());

                    Transpose(M,T);

                        // Third, record multiplication by b
                        // Fourth, divide both row elements by b
                    FactorOutRows(sout, M, nbadd, tev,
                                  i, M[i].begin(), M[i].end());

                    Transpose(T,M);
                    break;
            } } }
            if (found) break;
        } }
        if (found) break;
    } }
    } while(!over);
    return found;
}



// Direct program generateur from a matrix
template<typename outstream, typename _Mat, typename triple>
std::ostream& ProgramGen(outstream& sout, _Mat& M,
                         std::vector<triple>& multiples,
                         size_t& addcount, size_t& nbmul,
                         const char ouv, const char tev, const char rav) {

#ifdef VERBATIM_PARSING
    std::clog << "# Program Generation:" << ouv << ' ' << tev << ' ' << rav
              << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

    const auto& FF(M.field());

        // Factoring multiplier by columns
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

    Transpose(T,M);

        // Factoring triangle multiplier
    for(size_t j=0; j<M.coldim(); ++j) {
        Triangle(sout, M, T, multiples, addcount, nbmul, tev, rav, j);
    }

        // Computing remaining (simple) linear combinations
    for(size_t i=0; i<M.rowdim(); ++i) {
        const auto& row(M[i]);
        if (row.size()>0) {
            sout << ouv << i << ":=";

            auto arbs(Fabs(FF,row.begin()->second));

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
                if (! FF.areEqual(arbs,row.begin()->second)) sout << '-';
                sout << rav << rindex;
            } else {
                if ( (Fsign(FF,row.begin()->second) < 0)
                     || FF.isMOne(row.begin()->second) ) sout << '-';
                printmulorjustdiv(sout, tev, row.begin()->first,
                                  arbs, nbmul, FF);
            }

                // For all the monomials of the linear combination
            auto iter(row.begin());
            for(++iter; iter!= row.end(); ++iter) {
                ++addcount;
                auto ais(Fabs(FF,iter->second));

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
                    sout << (FF.areEqual(ais,iter->second) ? '+' : '-')
                         << rav << rindex;
                } else {
                        // otherwise next function will sout << '-'
                    sout << ( ( (Fsign(FF,iter->second) <0)
                                || FF.isMOne(iter->second) ) ? '-' : '+');
                    printmulorjustdiv(sout, tev, iter->first,
                                      ais, nbmul, FF);
                }
            }
            sout << ';' << std::endl;
        }
        else {
            sout << ouv << i << ":=0;" << std::endl;
        }
    }
#ifdef VERBATIM_PARSING
    std::clog << "# Program Generation done." << std::endl;
    std::clog << std::string(30,'#') << std::endl;
#endif

    return sout;
}


// Global random optimization function (pairs and factors)
template<typename outstream, typename _Mat>
Pair<size_t> Optimizer(outstream& sout, _Mat& M,
                       const char ouv, const char tev, const char rav) {

    using triple=std::tuple<size_t, size_t, typename _Mat::Element>;
    size_t nbadd(0), nbmul(0);

        // Factoring sums
    std::vector<triple> multiples;
    for( ; OneSub(sout, M, multiples, nbadd, nbmul, tev, rav) ; ) { }

        // No more useful CSE, thus
        // generate the rest of the program from the new M
    ProgramGen(sout, M, multiples, nbadd, nbmul, ouv, tev, rav);

    return Pair<size_t>(nbadd,nbmul);
}



	// Postcondition _Matrix A is nullified
template<typename outstream, typename _Mat>
Pair<size_t> nullspacedecomp(outstream& sout, _Mat& x, _Mat& A,
                             const bool mostCSE) {
	std::vector<size_t> l;
    return nullspacedecomp(sout, x, A, l, mostCSE);

    const auto& FF(A.field());
    const size_t Nj(A.coldim());
        // Add identity

    x.resize(x.rowdim()+A.rowdim(),x.coldim());
    A.resize(A.rowdim(),Nj+A.rowdim());
    for(size_t i=0; i<A.rowdim(); ++i)
        A.setEntry(i,Nj+i,FF.one);



//     std::clog << std::string(20,'#') << std::endl;
//     A.write(std::clog << "# Aug:\n", FileFormat::Pretty) << std::endl;
//     std::clog << std::string(20,'#') << std::endl;

    const auto P = nullspacedecomp(sout, x, A, l, mostCSE);

//     std::clog << sout.str() << std::flush;
//     std::clog << std::string(20,'#') << std::endl;


    return P;
}

	// Postcondition _Matrix A is nullified
template<typename outstream, typename _Mat>
Pair<size_t> nullspacedecomp(outstream& sout, _Mat& x, _Mat& A,
                             std::vector<size_t>& l, const bool mostCSE) {
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
        // If permutation is not fixed,
        // Randomly swap initial rows of FreePart
    if (l.size() != Nj) {
        l.resize(Nj);
        std::iota(l.begin(), l.end(), 0); // Select a random permutation
        std::shuffle ( l.begin(), l.end(),
                       std::default_random_engine(Givaro::BaseTimer::seed()));
    } // Otherwise only use the prescribed permutation
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
    LinBox::Permutation<typename _Mat::Field> P(FF,(long)Nj);
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
    std::clog << "# NullSpace dimensions:\t" << x.rowdim() << 'x' << x.coldim() << std::endl;
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
        Givaro::GivRandom generator;
        const size_t NotIndep(generator() % nullity);

        for(size_t i=0; i+NotIndep<nullity; ++i) {
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
        std::clog << "# FreePart dimensions:\t(" << Rank
                  << '+' << NotIndep << ')'
                  << 'x' << FreePart.coldim() << std::endl;
        std::clog << std::string(30,'#') << std::endl;
#endif

            // ============================================
            // Optimize the set of independent rows
        input2Temps(sout, FreePart.coldim(), 'i', 't');
            // o <-- Free . i
        Pair<size_t> Fops;
        if (mostCSE) {
            Fops = RecOptimizer(sout, FreePart, 'o', 't', 'r');
        } else {
            Fops = Optimizer(sout, FreePart, 'o', 't', 'r');
        }

            // ============================================
            // Optimize the dependent rows (transposed nullspace)
            // [ Free^T | Dep^T ] . [ x^T | I ]^T = 0
            // So Dep . i = (- x^T) Free . i = (- x^T) o
         _Mat Tx(x.field(), x.coldim(), x.rowdim()); NegTranspose(Tx, x);
         Tx.resize(Tx.rowdim()-NotIndep,Tx.coldim());

#ifdef VERBATIM_PARSING
        std::clog << std::string(30,'#') << std::endl;
        std::clog << "# Dependent dimensions (" << x.coldim() << '-' << NotIndep
                  << "):\t" << Tx.rowdim() << 'x' << Tx.coldim() << std::endl;
#endif

        Transpose(x, Tx);  // Remove the NotIndep columns of x for i2T checks
        input2Temps(sout, Tx.coldim(), 'o', 'v', x);
            // x <-- Tx . o = (- x^T) o
        Pair<size_t> Kops;
        if (mostCSE) {
            Kops = RecOptimizer(sout, Tx, 'x', 'v', 'g');
        } else {
            Kops = Optimizer(sout, Tx, 'x', 'v', 'g');
        }

            // ============================================
            // Recover final output of NullSpace
            // Applying both permutations, P then Q
            // back into 'o' variables
        for(size_t i=0; i+NotIndep<nullity; ++i) {
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


// Recusive search for the best cse
template<typename triple,typename _Mat>
bool RecSub(std::vector<std::string>& out, _Mat& Mat,
            std::vector<triple>& multiples,
            size_t& nbadd, size_t& nbmul, const size_t lvl,
            const char tev, const char rav) {

    const auto& FF(Mat.field());

    std::vector<std::vector<triple>> AllPairs;
    std::vector<size_t> Density;
        // Compute all pairs, in a row
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
    std::vector<std::string> bestdout;

#pragma omp parallel for shared(AllPairs,Mat,bestadds,bestmuls,bestM,bestmultiples,bestdout,tev,rav,multiples)
    for(size_t i=0; i<AllPairs.size(); ++i) {
        const auto& rows(AllPairs[i]);
        for (const auto& cse: rows) {
            if (PairMap[cse] > 1) {
                _Mat lM(FF,Mat.rowdim(),Mat.coldim()); sparse2sparse(lM, Mat);
                std::vector<triple> lmultiples(multiples);
                std::ostringstream ssout;
                size_t moremul(0);

                    // Factoring out that CSE from the matrix
                const Pair<long> savings = RemOneCSE(ssout, lM, moremul,
                                                     lmultiples, cse,
                                                     AllPairs, tev, rav);

                size_t ladditions(nbadd-savings.first);
                size_t lmuls(nbmul-savings.second);

                std::vector<std::string> sdout(1,ssout.str());
                RecSub(sdout, lM, lmultiples, ladditions, lmuls, lvl+1, tev, rav);

#pragma omp critical
                {
                    if ( (ladditions < bestadds) ||
                         ( (ladditions == bestadds) && (lmuls < bestmuls) ) ) {
                        bestadds = ladditions;
                        bestmuls = lmuls;
                        sparse2sparse(bestM, lM);
                        bestmultiples = lmultiples;
                        bestdout = sdout;
                    }
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



// Global exhaustive greedy CSE optimization function (pairs and factors)
template<typename outstream, typename _Mat>
Pair<size_t> RecOptimizer(outstream& sout, _Mat& M,
                          const char ouv, const char tev, const char rav) {
    using triple=std::tuple<size_t, size_t, typename _Mat::Element>;
    const auto& FF(M.field());

    Pair<size_t> nbops(naiveOps(M));

        // Factoring sums
    std::vector<std::string> dout;
    std::vector<triple> multiples;
    RecSub(dout, M, multiples, nbops.first, nbops.second, 0, tev, rav);

    size_t alreadymuls(nbops.second);
    for(auto it = M.IndexedBegin(); it != M.IndexedEnd(); ++it)
        if (notAbsOne(FF,it.value())) --alreadymuls;

    for(auto iter = dout.begin(); iter != dout.end(); ++iter) {
        sout << *iter;
    }

    size_t addcount(0);

    ProgramGen(sout, M, multiples, addcount, alreadymuls, ouv, tev, rav);

    return Pair<size_t>(nbops.first,alreadymuls);
}




// ============================================================
// Optimizing a linear program (Gaussian elimination method)
template<typename Field>
Pair<size_t>& LUOptimiser(Pair<size_t>& nbops, std::ostringstream& sout,
                          const Field& F, const Matrix& M, const Matrix& T,
                          Givaro::Timer& global, const size_t randomloops) {

        // ============================================================
        // Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;
    Givaro::Timer chrono; chrono.start();

        // ============================================================

    FMatrix U(M, F);

#ifdef DEBUG
    U.write(std::clog << "## M:=Matrix(", FileFormat(8)) << ");" << std::endl;
#endif

    typename Field::Element Det;
    size_t Rank;
    FMatrix L(F, U.rowdim(), U.rowdim());
    LinBox::Permutation<Field> Q(F,U.rowdim());
    LinBox::Permutation<Field> P(F,U.coldim());
    LinBox::GaussDomain<Field> GD(F);
    GD.QLUPin(Rank, Det, Q, L, U, P, U.rowdim(), U.coldim() );

#ifdef DEBUG
    Q.write(std::clog << "## Q:=Matrix(", FileFormat(8)) << ");" << std::endl;
    L.write(std::clog << "## L:=Matrix(", FileFormat(8)) << ");" << std::endl;
    U.write(std::clog << "## U:=Matrix(", FileFormat(8)) << ");" << std::endl;
    P.write(std::clog << "## P:=Matrix(", FileFormat(8)) << ");" << std::endl;
#endif

    std::ostringstream gout;
    auto gops(nbops);

#pragma omp parallel for shared(Q,L,U,P,gout,gops)
    for(size_t i=0; i<randomloops; ++i) {
        std::ostringstream luout;
        FMatrix lU(U, F);
        FMatrix lL(L, F);

            // ============================================
            // Applying permutation P to 'i' variables
        input2Temps(luout, P.rowdim(), 'i', 't', P.getStorage());

            // ============================================
            // Applying Upper matrix to 't' variables
        auto Uops = Optimizer(luout, lU, 'v', 't', 'r');

            // ============================================
            // Applying Lower matrix to 'v' variables
        auto Lops = Optimizer(luout, lL, 'x', 'v', 'g');

            // ============================================
            // Applying permutation Q back into 'o' variables
        input2Temps(luout, Q.rowdim(), 'x', 'o', Q.getStorage());

        Uops.first += Lops.first;
        Uops.second += Lops.second;

#pragma omp critical
        {
            const bool better(cmpOpCount(Uops,gops));
            if ( (gout.tellp() == std::streampos(0) ) || better ) {
                gout.clear(); gout.str(std::string());
                gout.str(luout.str());
                if (better) {
                    std::clog << "# Found G: "
                              << Uops.first << '|' << Uops.second
                              << " instead of "
                              << gops.first << '|' << gops.second
                              << std::endl;
                    gops = Uops;
                }
            }
        }
    }

    chrono.stop(); global += chrono;
    if (cmpOpCount(gops,nbops)) {
        nbops = gops;
        sout.clear(); sout.str(std::string());
        sout << gout.str();
    }
    return nbops;
}

// ============================================================
// Optimizing a linear program (Alternative Factrization method)
template<typename Field>
Pair<size_t>& ABOptimiser(Pair<size_t>& nbops, std::ostringstream& sout,
                          const Field& F, const Matrix& M, const Matrix& T,
                          const size_t selectinnerdim,
                          Givaro::Timer& global, const size_t randomloops) {

        // ============================================================
        // Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;
    Givaro::Timer chrono; chrono.start();

    FMatrix FM(M, F);

    FMatrix Alt(F, FM.rowdim(), selectinnerdim);
    FMatrix CoB(F, selectinnerdim, FM.coldim());

    const size_t factoloops(1+randomloops>>3);
    Factorizer(Alt, CoB, FM, factoloops, selectinnerdim, false);

    const size_t sa(density(Alt)), sc(density(CoB));

    std::ostringstream gout;
    auto gops(nbops);

#pragma omp parallel for shared(Alt,CoB,sa,sc,gout,gops)
    for(size_t i=0; i<randomloops; ++i) {
        FMatrix lC(CoB, F);
        FMatrix lA(Alt, F);

        std::ostringstream about;
        input2Temps(about, lC.coldim(), 'i', 't');

            // ============================================
            // Applying CoB matrix from 't' to 'v' variables
        auto Bops = Optimizer(about, lC, 'v', 't', 'r');
            // ============================================
            // Applying Alt matrix from 'v' to 'o' variables
        auto Aops = Optimizer(about, lA, 'o', 'v', 'g');

        Pair<size_t> ABops(Bops.first + Aops.first,
			   Bops.second + Aops.second);

#pragma omp critical
        {
            const bool better(cmpOpCount(ABops,gops));
            if ( better || (gout.tellp() == std::streampos(0) ) ) {
                gout.clear(); gout.str(std::string());
                gout << about.str();
                if (better) {
                    std::clog << "# Found A: ("
                              << Alt.rowdim() << 'x' << Alt.coldim() << 'x'
                              << CoB.coldim() << ' ' << sa << '/' << sc
                              << ')' << '\t' << Bops.first << '+' << Aops.first
			      << '|' << Bops.second << '+' << Aops.second
                              << " instead of "
                              << gops.first << '|' << gops.second
                              << std::endl;
                    gops = ABops;
                }
            }
        }
    }

    chrono.stop(); global += chrono;
    if (cmpOpCount(gops,nbops)) {
        nbops = gops;
        sout.clear(); sout.str(std::string());
        sout << gout.str();
    }
    return nbops;
}



// ============================================================
// Optimizing a linear program (Common Subexpressions Elimination)
template<typename Field>
Pair<size_t>& CSEOptimiser(Pair<size_t>& nbops, std::ostringstream& sout,
                           const Field& F, const Matrix& M, const Matrix& T,
                           Givaro::Timer& global, const size_t randomloops) {
        // ============================================================
        // Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;

    Givaro::Timer chrono; chrono.start();
    std::ostringstream dout;
    auto dops(nbops);

#pragma omp parallel for shared(M,T,dout,dops)
    for(size_t i=0; i<randomloops; ++i) {
        FMatrix lM(M, F);
        FMatrix lT(T, F);

        std::ostringstream ldout;
            // Cancellation-free optimization
        input2Temps(ldout, lM.coldim(), 'i', 't', lT);
        auto lnbops( Optimizer(ldout, lM, 'o', 't', 'r') );


#pragma omp critical
        {
#ifdef VERBATIM_PARSING
            std::clog << "# Found, direct: "
                      << lnbops.first << "\tadditions, "
                      << lnbops.second << "\tmultiplications." << std::endl;
#endif
            const bool better(cmpOpCount(lnbops,dops));
            if ( (dout.tellp() == std::streampos(0)) || better ) {
                dout.clear(); dout.str(std::string());
                dout << ldout.str();
                    if (better) {
                        std::clog << "# Found D: "
                                  << lnbops.first << '|' << lnbops.second
                                  << " instead of "
                                  << dops.first << '|' << dops.second
                                  << std::endl;
                        dops = lnbops;
                    }
            }
        }
    }

    chrono.stop(); global += chrono;
    if (cmpOpCount(dops,nbops)) {
        nbops = dops;
        sout.clear(); sout.str(std::string());
        sout << dout.str();
    }
    return nbops;
}

// ============================================================
// Optimizing a linear program (whole CSE tree)
template<typename Field>
Pair<size_t>& AllCSEOpt(Pair<size_t>& nbops, std::ostringstream& sout,
                        const Field& F, const Matrix& M, const Matrix& T,
                        Givaro::Timer& global) {
        // ============================================================
        // Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;

    Givaro::Timer chrono; chrono.start();
    FMatrix lM(M, F);
    FMatrix lT(T, F);

    std::ostringstream iout;
        // Cancellation-free optimization
    input2Temps(iout, lM.coldim(), 'i', 't', lT);
    auto rnbops( RecOptimizer(iout, lM, 'o', 't', 'r') );

    chrono.stop(); global += chrono;

    if (cmpOpCount(rnbops, nbops)) {
        sout.clear(); sout.str(std::string());
        sout << iout.str();
        nbops = rnbops;
    } else {
            // Optimal was already found
        std::clog << "# \033[1;36m"
                  << "No greedy CSE schedule has less additions.\033[0m \t"
                  << chrono << std::endl;
    }
    return nbops;
}



// ============================================================
// Optimizing a linear program (kernel method)
template<typename Field>
Pair<size_t>& KernelOptimiser(Pair<size_t>& nbops, std::ostringstream& sout,
                              const Field& F, const Matrix& T,
                              Givaro::Timer& global, const size_t randomloops) {
        // ============================================================
        // Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;
    Givaro::Timer chrono; chrono.start();
    std::ostringstream kout;
    auto kops(nbops);
    bool nonzerokernel(true);

#pragma omp parallel for shared(T,kout,kops)
    for(size_t i=0; i<randomloops; ++i) {
        if (nonzerokernel) {
            FMatrix lT(T, F);
            std::ostringstream lkout;
            FMatrix NullSpace(F,lT.coldim(),lT.coldim());
            auto lnbops( nullspacedecomp(lkout, NullSpace, lT) );

#pragma omp critical
            {
#ifdef VERBATIM_PARSING
                std::clog << "# Found, kernel: "
                          << lnbops.first << "\tadditions, "
                          << lnbops.second << "\tmultiplications." << std::endl;
#endif
                if (lnbops == Pair<size_t>{-1,-1}) {
                    nonzerokernel = false;
                }
                const bool better(cmpOpCount(lnbops,kops));
                if ( (kout.tellp() == std::streampos(0) ) || better ) {
                    kout.clear(); kout.str(std::string());
                    kout << lkout.str();
                    if (better) {
                        std::clog << "# Found K: "
                                  << lnbops.first << '|' << lnbops.second
                                  << " instead of "
                                  << kops.first << '|' << kops.second
                                  << std::endl;
                        kops = lnbops;
                    }
                }
            }
        }
    }

    chrono.stop(); global += chrono;
    if (! nonzerokernel) {
            // Zero dimensional kernel
        std::clog << "# \033[1;36mZero dimensional kernel.\033[0m"
                  << std::endl;
    } else if (cmpOpCount(kops,nbops)) {
        nbops = kops;
        sout.clear(); sout.str(std::string());
        sout << kout.str();
    }
    return nbops;
}

// ============================================================
// Optimizing a linear program (kernel method)
template<typename Field>
Pair<size_t>& AllKernelOpt(Pair<size_t>& nbops, std::ostringstream& sout,
                           const Field& F, const Matrix& T, const bool mostCSE,
                           Givaro::Timer& global, const size_t randomloops) {
        // ============================================================
        // Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;
    Givaro::Timer chrono; chrono.start();
    const size_t m(T.coldim());
    std::vector<size_t> Fm { factorial(m) };
    std::ostringstream aout;
    auto aops(nbops);

#pragma omp parallel for shared(Fm,T,aout,aops)
    for(size_t i=0; i<Fm.back(); ++i) {
        std::vector<size_t> l{kthpermutation(i,m,Fm)};
        FMatrix lT(T, F);
        std::ostringstream laout;
        FMatrix NullSpace(F,lT.coldim(),T.coldim());
        auto lkops( nullspacedecomp(laout, NullSpace, lT, l, mostCSE) );

#pragma omp critical
        {
#ifdef VERBATIM_PARSING
            std::clog << "# Found, kernel: "
                      << lkops.first << "\tadditions, "
                      << lkops.second << "\tmultiplications." << std::endl;
#endif
            const bool better(cmpOpCount(lkops,aops));
            if ( (aout.tellp() == std::streampos(0)) || better) {
                aout.clear(); aout.str(std::string());
                aout << laout.str();
                if (better) {
                    std::clog << "# Found E: "
                              << lkops.first << '|' << lkops.second
                              << " instead of "
                              << aops.first << '|' << aops.second
                              << std::endl;
                    aops = lkops;
                }
            }
        }
    }

    chrono.stop(); global += chrono;

    if (cmpOpCount(aops,nbops)) {
        sout.clear(); sout.str(std::string());
        sout << aout.str();
        nbops = aops;
    } else {
            // Optimal was already found
        std::clog << "# \033[1;36m"
                  << "No kernel permutation has less additions.\033[0m \t"
                  << chrono << std::endl;
    }
    return nbops;
}


// ============================================================
// Optimizing a linear program (Direct CSE or Kernel methods)
template<typename Field>
Pair<size_t> OptMethods(const Pair<size_t> opsinit,
                        std::ostringstream& sout, const Field& F,
                        const Matrix& M, const Matrix& T,
                        Givaro::Timer& global, const size_t randomloops,
                        const bool printMaple, const bool printPretty,
                        const bool tryDirect, const bool tryKernel,
                        const bool tryLU, bool tryAB,
                        const bool mostCSE, const bool allkernels) {

    Pair<size_t> nbops(opsinit);

        // ============================================================
        // Optimize the Alternative factorization
    if (tryAB) {
//         for(int id(M.rowdim()-1); id >= M.coldim(); --id)
        for(size_t id(M.coldim()); id < M.rowdim(); ++id)
            ABOptimiser(nbops, sout, F, M, T, id, global, randomloops);
    }


        // ============================================================
        // Otimize the whole matrix
    if (tryDirect) {
        CSEOptimiser(nbops, sout, F, M, T, global, randomloops);
    }

        // ============================================================
        // Separate independent and dependent rows
    if (tryKernel) {
        KernelOptimiser(nbops, sout, F, T, global, randomloops);
    }

        // ============================================================
        // Optimize the factorized matrix
    if (tryLU) {
        LUOptimiser(nbops, sout, F, M, T, global, randomloops);
    }

        // ============================================================
        // Exhaustive nullspace permutation search (if # is <= 12!)
    if (allkernels && (M.rowdim() < 13)) {
        AllKernelOpt(nbops, sout, F, T, mostCSE, global, randomloops);
    }

        // ============================================================
        // Greedy CSE search
    if (mostCSE) {
        AllCSEOpt(nbops, sout, F, M, T, global);
    }

    if (cmpOpCount(opsinit,nbops) || (opsinit == nbops)) {
            // Found nothing better than direct matrix expression
        typedef typename Matrix::template rebind<Field>::other FMatrix;
        FMatrix lM(M, F), lT(T, F);
        size_t nbadd(0), nbmul(0);
        std::vector<Etriple<Field>> multiples;

        input2Temps(sout, lM.coldim(), 'i', 't', lT);
        ProgramGen(sout, lM, multiples, nbadd, nbmul, 'o', 't', 'r');

        nbops.first = nbadd;
        nbops.second= nbmul;
    }

    return nbops;
}


} // End of namespace PLinOpt
// ============================================
