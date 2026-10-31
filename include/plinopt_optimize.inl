// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * PLinOpt Library, Optimization inline implementations
 ****************************************************************/

#include <limits>
#include <unordered_map>
#include "plinopt_optimize.h"
#include "plinopt_sparsify.h"
#include "plinopt_programs.h"




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
		     const char tev, const char rav, const bool updateAPM,
		     triple* applied = nullptr) {
    const auto& FF(lM.field());
    size_t lm(lM.coldim());
    long savedadds(0), savedmuls(0);

	// Looking for the most number of ones in either columns of cse
	// (scanning the rows: a transposed copy per CSE was a fifth of the
	// optimizer's time)
    size_t count0(0), count1(0);
    for(auto rit=lM.rowBegin(); rit != lM.rowEnd(); ++rit) {
	for(const auto& iter: *rit) {
	    if (iter.first == std::get<0>(cse)) {
		if (isAbsOne(FF,iter.second)) ++count0;
	    } else if (iter.first == std::get<1>(cse)) {
		if (isAbsOne(FF,iter.second)) ++count1;
	    }
	}
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
    if (applied) *applied = lcse;      // the new temporary is t_a + c t_b

	// Factor out cse, in all rows
	// adds a column for the new factor
    typename _Mat::Element ratio; FF.init(ratio);
    for(size_t i=0; i<AllPairs.size(); ++i) {
	    // Is cse a pair of row i: its two columns, in this order, with
	    // this ratio (what listpairs would have listed)
	bool found(false);
	{
	    auto first(lM[i].end()), second(lM[i].end());
	    for(auto iter=lM[i].begin(); iter != lM[i].end(); ++iter) {
		if (iter->first == std::get<0>(cse)) first = iter;
		else if (iter->first == std::get<1>(cse)) second = iter;
	    }
	    if ( (first != lM[i].end()) && (second != lM[i].end())
		 && (first < second) ) {
		FF.div(ratio, second->second, first->second);
		found = FF.areEqual(ratio, std::get<2>(cse));
	    }
	}
	if (found) {
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
		    // Update AllPairs[i] and PairMap: only the pairs with one
		    // of the merged columns leave, only the pairs with the
		    // new column arrive
		std::vector<triple> newrow;
		for (const auto& iter: AllPairs[i]) {
		    if ( (std::get<0>(iter) != std::get<0>(lcse)) &&
			 (std::get<1>(iter) != std::get<0>(lcse)) &&
			 (std::get<0>(iter) != std::get<1>(lcse)) &&
			 (std::get<1>(iter) != std::get<1>(lcse))
			 ) {
			newrow.push_back(iter);
		    } else {
			auto pit(PairMap.find(iter));
			if (pit != PairMap.end()) {
			    if (--(pit->second) == 0) PairMap.erase(pit);
			}
		    }
		}
		     // Add all pairs with new (last) element
		typename _Mat::Element tmp; FF.init(tmp);
		auto endlMi(lM[i].end()); --endlMi;
		for(auto iter=lM[i].begin(); iter != endlMi; ++iter) {
		    newrow.emplace_back(iter->first, endlMi->first,
					FF.div(tmp, endlMi->second, iter->second));
		    PairMap[newrow.back()]++;
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


// Removing CSEs until no pair occurs twice
template<typename triple,typename _Mat>
void OneSub(std::ostream& sout, _Mat& M, std::vector<triple>& multiples,
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

    if (PairMap.size()==0) return;

	// The loop below removes CSEs until no pair occurs twice, updating
	// AllPairs and PairMap exactly as it goes (RemOneCSE); so when it
	// ends there is nothing left for another call to find, and it
	// returns false: Optimizer need not rebuild the pair lists to
	// learn that.
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

	if (maxfrq <=1) return;
	    // Factoring will gain something if maximal frequency > 1

	    // More than one pair with maximal frequency
	if (MaxCSE.size()>1) {
		// Tie breaking heuristics
#if defined(RANDOM_TIES) && !defined(DENSITY_OPTIMIZATION)
	    cse = MaxCSE[ generator()() % MaxCSE.size() ];
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

#if VERBATIM_PARSING > 1
#  if defined(RANDOM_TIES) && !defined(DENSITY_OPTIMIZATION)
	printEtriple(std::clog << "# Found: ", FF, cse) << '=' << maxfrq
			       << ", Saved: " << savings.first << "+|"
			       << savings.second << 'x' << std::endl;
#  else
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
    OneSub(sout, M, multiples, nbadd, nbmul, tev, rav);

	// No more useful CSE, thus
	// generate the rest of the program from the new M
    ProgramGen(sout, M, multiples, nbadd, nbmul, ouv, tev, rav);

    return Pair<size_t>(nbadd,nbmul);
}



	// Postcondition _Matrix A is nullified
template<typename outstream, typename _Mat>
Pair<size_t> nullspacedecomp(outstream& sout, _Mat& x, _Mat& A,
			     const bool KFI, const bool mostCSE) {
	std::vector<size_t> l;
    if (KFI) {
	const auto& FF(A.field());
	const size_t Ni(A.rowdim());
	const size_t Nj(A.coldim());

	    // Add identity to goals
	x.resize(x.rowdim()+A.rowdim(),x.coldim());
	A.resize(A.rowdim(),Nj+A.rowdim());
	for(size_t i=0; i<A.rowdim(); ++i)
	    A.setEntry(i,Nj+i,FF.one);

	auto Pops = nullspacedecomp(sout, x, A, l, 'q', 'i', mostCSE);

	std::vector<std::string> iVars(A.rowdim());
	for(size_t i=0; i<A.rowdim(); ++i) {
	    iVars[i]="q"+std::to_string(Nj+i);
	}

	std::stringstream ssin; ssin << sout.str();
	VProgram_t progV; programParser(progV, ssin);
	bool isreduced(false);
	while( OpOnVars(progV, iVars) ) {
	    isreduced = true;
	    RemoveVars(progV, Pops, iVars);
	}

	if (isreduced) {
	    sout.clear(); sout.str(std::string());
	    char next(0u); parenthesisExpand(progV, next);
		// Make output suitable for potential transpose
	    input2Temps(sout, Ni, 'i', next);
	    for(auto& line: progV) for(auto& word: line)
		if (word[0]=='i') word[0]=next;
	    sout << progV;
	}

	    // out all goals, except the added identity
	    // Warning: the identity goals could now use unnecessary ops
	for(size_t j=0; j<Nj; ++j)
	    sout << 'o' << j << ":=" << 'q' << j << ';' << std::endl;

	return Pops;
    } else {
	return nullspacedecomp(sout, x, A, l, 'o', 'i', mostCSE);
    }
}

	// Postcondition _Matrix A is nullified
template<typename outstream, typename _Mat>
Pair<size_t> nullspacedecomp(outstream& sout, _Mat& x, _Mat& A,
			     std::vector<size_t>& l, const char ouv,
			     const char inv, const bool mostCSE) {
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


		       std::default_random_engine(generator()()));
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
	    // Remove (randomly chosen) dependent rows from FreePart
#if defined(RANDOM_TIES)
	const bool restrictNotIndep(true);
#else
	const bool restrictNotIndep(KFI);
#endif
	size_t FreePartKept(nullity), NotIndep(0u);
	if (restrictNotIndep) {
	    NotIndep = (generator()() % nullity);
	    FreePartKept -= NotIndep;
	}

	for(size_t i=0; i<FreePartKept; ++i) {
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
#pragma omp critical
	{
	std::clog << "# FreePart dimensions:\t(" << Rank
		  << '+' << NotIndep << ')'
		  << 'x' << FreePart.coldim() << std::endl;
	std::clog << std::string(30,'#') << std::endl;
	}
#endif

	    // ============================================
	    // Optimize the set of independent rows
	input2Temps(sout, FreePart.coldim(), inv, 't');
	    // o <-- Free . i
	Pair<size_t> Fops;
	if (mostCSE) {
	    Fops = RecOptimizer(sout, FreePart, ouv, 't', 'r');
	} else {
	    Fops = Optimizer(sout, FreePart, ouv, 't', 'r');
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
	input2Temps(sout, Tx.coldim(), ouv, 'v', x);
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
	for(size_t i=0; i<FreePartKept; ++i) {
	    sout << ouv << Q[ P.getStorage()[ Rank+i ] ] << ":="
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


// ============================================================
// Exhaustive greedy CSE search, memoized over canonical states, with a budget
//   A state is the current matrix, the temporaries created so far and the
//   scalings already emitted.  Each temporary is a linear form over the
//   inputs; two paths creating the same temporaries in different orders
//   reach the same state up to their numbering, and the key below numbers
//   them canonically, so every state is searched once.  The memo keeps the
//   best savings from a state and the CSE achieving them; the program is
//   then replayed from the root along the best choices.
//   RecBudget bounds the number of states searched: past it, a state keeps
//   only its first candidate (a greedy completion) and the search is
//   reported as not exhaustive.
size_t RecBudget(100000);

template<typename _Mat>
struct RecMemo {
    typedef typename _Mat::Element Element;
    typedef std::vector<std::pair<size_t,Element>> Form;    // over the inputs
    struct Entry { long adds; long muls; std::string choice; };
    std::unordered_map<std::string, Entry> table;
    size_t inputs, states;
    bool exhaustive;
    RecMemo(const size_t n) : inputs(n), states(0), exhaustive(true) {}
};

template<typename Field, typename Form>
std::string formKey(const Field& F, const Form& v) {
    std::ostringstream o;
    for(const auto& [k, c]: v) { F.write(o << k << ':', c) << ';'; }
    return o.str();
}

    // t_a + c t_b, as a form over the inputs (sorted, without zeroes)
template<typename Field, typename Form>
Form formCombine(const Field& F, const Form& a, const Form& b,
		 const typename Field::Element& c) {
    std::map<size_t, typename Field::Element> m;
    for(const auto& [k, v]: a) {
	auto it(m.find(k));
	if (it == m.end()) m.emplace(k, v); else F.addin(it->second, v);
    }
    typename Field::Element t; F.init(t);
    for(const auto& [k, v]: b) {
	F.mul(t, c, v);
	auto it(m.find(k));
	if (it == m.end()) m.emplace(k, t); else F.addin(it->second, t);
    }
    Form r;
    for(const auto& [k, v]: m) if (!F.isZero(v)) r.emplace_back(k, v);
    return r;
}

    // The canonical key of a state
template<typename triple, typename _Mat, typename Memo>
std::string stateKey(const _Mat& M, const std::vector<triple>& multiples,
		     const std::vector<typename Memo::Form>& colvec,
		     const Memo& memo) {
    const auto& F(M.field());
	// temporaries in the order of their forms
    std::vector<std::pair<std::string,size_t>> temps;
    for(size_t j=memo.inputs; j<colvec.size(); ++j)
	temps.emplace_back(formKey(F, colvec[j]), j);
    std::sort(temps.begin(), temps.end());
    std::vector<size_t> canon(colvec.size());
    for(size_t j=0; j<memo.inputs; ++j) canon[j] = j;
    for(size_t k=0; k<temps.size(); ++k) canon[temps[k].second] = memo.inputs + k;
    std::ostringstream o;
    for(const auto& t: temps) o << t.first << '|';
    o << '#';
    for(auto it=M.rowBegin(); it != M.rowEnd(); ++it) {
	std::vector<std::pair<size_t,typename _Mat::Element>> row;
	for(const auto& e: *it) row.emplace_back(canon[e.first], e.second);
	std::sort(row.begin(), row.end(),
		  [](const auto& x, const auto& y){ return x.first < y.first; });
	for(const auto& [k, c]: row) { F.write(o << k << ':', c) << ','; }
	o << '/';
    }
    o << '#';
    std::vector<std::string> muls;
    for(const auto& m: multiples) {
	std::ostringstream mo; F.write(mo << canon[std::get<1>(m)] << '*', std::get<2>(m));
	muls.push_back(mo.str());
    }
    std::sort(muls.begin(), muls.end());
    for(const auto& m: muls) o << m << ',';
    return o.str();
}

template<typename Field, typename triple, typename Form>
std::string choiceKey(const Field& F, const triple& cse,
		      const std::vector<Form>& colvec) {
    std::ostringstream o;
    F.write(o << formKey(F, colvec[std::get<0>(cse)]) << '|'
	      << formKey(F, colvec[std::get<1>(cse)]) << '|', std::get<2>(cse));
    return o.str();
}

    // Best savings from a state (filling the memo)
template<typename triple, typename _Mat, typename Memo>
Pair<long> RecSearch(const _Mat& Mat, const std::vector<triple>& multiples,
		     const std::vector<typename Memo::Form>& colvec,
		     Memo& memo, const size_t lvl,
	    const char tev, const char rav) {
    const auto& FF(Mat.field());
    const std::string key(stateKey<triple>(Mat, multiples, colvec, memo));
    bool known(false); Pair<long> result{0,0};
#pragma omp critical(recmemo)
    {
	const auto it(memo.table.find(key));
	if (it != memo.table.end()) { known = true; result = Pair<long>{it->second.adds, it->second.muls}; }
    }
    if (known) return result;

	// The candidates: pairs occurring in at least two rows
    std::vector<std::vector<triple>> AllPairs;
    for(auto iter=Mat.rowBegin(); iter != Mat.rowEnd(); ++iter)
	AllPairs.push_back(listpairs(*iter, FF));
    std::map<triple,size_t> PairMap;
    for(const auto& rows: AllPairs) for (const auto& p: rows) PairMap[p]++;
    std::vector<triple> candidates;
    for(const auto& [p, frequency]: PairMap) if (frequency > 1) candidates.push_back(p);

    bool over(false);
#pragma omp critical(recmemo)
    {
	++memo.states;
	if (memo.states > RecBudget) { over = true; memo.exhaustive = false; }
    }
    if (over && candidates.size() > 1) candidates.resize(1);

    long bestadds(-1), bestmuls(-1); std::string bestchoice;
    if (!candidates.empty()) {
	bestadds = std::numeric_limits<long>::min(); bestmuls = bestadds;
	    // parallel at the root only (nested regions are not active)
#pragma omp parallel for shared(AllPairs,candidates,bestadds,bestmuls,bestchoice) if(lvl==0)
	for(size_t k=0; k<candidates.size(); ++k) {
	    const triple& cse(candidates[k]);
		_Mat lM(FF,Mat.rowdim(),Mat.coldim()); sparse2sparse(lM, Mat);
		std::vector<triple> lmultiples(multiples);
	    std::vector<typename Memo::Form> lcolvec(colvec);
	    std::vector<std::vector<triple>> lAllPairs(AllPairs);
	    std::ostringstream ssout; size_t moremul(0);
	    triple applied;
	    std::map<triple,size_t> lPairMap;
	    const Pair<long> savings = RemOneCSE(ssout, lM, moremul, lmultiples, cse,
						 lAllPairs, lPairMap, tev, rav, false, &applied);
	    lcolvec.push_back(formCombine(FF, lcolvec[std::get<0>(applied)],
					  lcolvec[std::get<1>(applied)], std::get<2>(applied)));
	    const Pair<long> rest(RecSearch(lM, lmultiples, lcolvec, memo, lvl+1, tev, rav));
	    const long adds(savings.first + rest.first), muls(savings.second + rest.second);
#pragma omp critical(recbest)
		{
		if ( (adds > bestadds) || ((adds == bestadds) && (muls > bestmuls)) ) {
		    bestadds = adds; bestmuls = muls;
		    bestchoice = choiceKey(FF, cse, colvec);
		}
		}
	    }
    } else {
	bestadds = 0; bestmuls = 0;
	}
#pragma omp critical(recmemo)
    {
	memo.table[key] = typename Memo::Entry{bestadds, bestmuls, bestchoice};
    }
    return Pair<long>{bestadds, bestmuls};
    }

    // Replays the best choices from the root, emitting the program
template<typename triple, typename _Mat, typename Memo>
Pair<long> RecReplay(std::ostream& sout, _Mat& Mat, std::vector<triple>& multiples,
		     std::vector<typename Memo::Form>& colvec, Memo& memo,
		     const char tev, const char rav) {
    const auto& FF(Mat.field());
    Pair<long> total{0,0};
    for(;;) {
	const std::string key(stateKey<triple>(Mat, multiples, colvec, memo));
	const auto it(memo.table.find(key));
	if ((it == memo.table.end()) || it->second.choice.empty()) break;
	    // the chosen CSE, by the forms of its two columns
	const std::string& ch(it->second.choice);
	const size_t s1(ch.find('|')), s2(ch.find('|', s1+1));
	const std::string ka(ch.substr(0,s1)), kb(ch.substr(s1+1, s2-s1-1));
	size_t a(colvec.size()), b(colvec.size());
	for(size_t j=0; j<colvec.size(); ++j) {
	    const std::string fj(formKey(FF, colvec[j]));
	    if ((a == colvec.size()) && (fj == ka)) a = j;
	    else if ((b == colvec.size()) && (fj == kb)) b = j;
	}
	triple cse; std::get<0>(cse) = a; std::get<1>(cse) = b;
	std::istringstream ratio(ch.substr(s2+1)); FF.read(ratio, std::get<2>(cse));
	std::vector<std::vector<triple>> AllPairs;
	for(auto iter=Mat.rowBegin(); iter != Mat.rowEnd(); ++iter)
	    AllPairs.push_back(listpairs(*iter, FF));
	size_t moremul(0); triple applied;
	std::map<triple,size_t> PairMap;
	const Pair<long> savings = RemOneCSE(sout, Mat, moremul, multiples, cse,
					     AllPairs, PairMap, tev, rav, false, &applied);
	colvec.push_back(formCombine(FF, colvec[std::get<0>(applied)],
				     colvec[std::get<1>(applied)], std::get<2>(applied)));
	total.first += savings.first; total.second += savings.second;
    }
    return total;
}

// Global exhaustive greedy CSE optimization function (pairs and factors)
template<typename outstream, typename _Mat>
Pair<size_t> RecOptimizer(outstream& sout, _Mat& M,
			  const char ouv, const char tev, const char rav) {
    using triple=std::tuple<size_t, size_t, typename _Mat::Element>;
    const auto& FF(M.field());

    Pair<size_t> nbops(naiveOps(M));

	// Factoring sums: the memoized search, then the replay of its choices
    std::vector<triple> multiples;
    RecMemo<_Mat> memo(M.coldim());
    std::vector<typename RecMemo<_Mat>::Form> colvec(M.coldim());
    for(size_t j=0; j<M.coldim(); ++j) colvec[j].emplace_back(j, FF.one);

    const Pair<long> best(RecSearch<triple>(M, multiples, colvec, memo, 0, tev, rav));
    const Pair<long> replayed(RecReplay<triple>(sout, M, multiples, colvec, memo, tev, rav));
    if (replayed != best) {
	std::clog << "# \033[1;31mWARNING: replayed savings " << replayed.first << '|'
		  << replayed.second << " differ from the searched " << best.first
		  << '|' << best.second << "\033[0m" << std::endl;
    }
    std::clog << "# Exhaustive CSE search: " << memo.states << " states"
	      << (memo.exhaustive ? "" : " (budget reached, not exhaustive)")
	      << std::endl;
    nbops.first -= replayed.first;
    nbops.second -= replayed.second;

    size_t alreadymuls(nbops.second);
    for(auto it = M.IndexedBegin(); it != M.IndexedEnd(); ++it)
	if (notAbsOne(FF,it.value())) --alreadymuls;

    size_t addcount(0);

    ProgramGen(sout, M, multiples, addcount, alreadymuls, ouv, tev, rav);

    return Pair<size_t>(nbops.first,alreadymuls);
}




// ============================================================
// Optimizing a linear program (Gaussian elimination method)
template<typename Field>
Pair<size_t>& LUOptimiser(Pair<size_t>& nbops, std::ostringstream& sout,
			  const Field& F, const Matrix& M, const Matrix& T,
			  Givaro::Timer& global, const size_t randomloops, const int verbose) {

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
    auto gops(nbops), sops(nbops);
    size_t besti(std::numeric_limits<size_t>::max());

#pragma omp parallel for shared(Q,L,U,P,gout,gops,sops,besti)
    for(size_t i=0; i<randomloops; ++i) {
	reseed('G', i);
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
	const bool empty(gout.tellp() == std::streampos(0));
	if ( empty || cmpOpCount(Uops,sops)
	     || (eqOpCount(Uops,sops) && (i < besti)) ) {
	    gout.clear(); gout.str(std::string());
	    gout.str(luout.str());
	    sops = Uops; besti = i;
	    if (better) {
		if (verbose>0) {
		    std::clog << "# Found G: "
			      << Uops.first << '|' << Uops.second
			      << " instead of "
			      << gops.first << '|' << gops.second << "\t["
			      << i << '/' << omp_get_thread_num() << ']'
			      << std::endl;
		}
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
			  Givaro::Timer& global, const size_t randomloops, const int verbose) {

	// ============================================================
	// Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;
    Givaro::Timer chrono; chrono.start();

    FMatrix FM(M, F);

    FMatrix Alt(F, FM.rowdim(), selectinnerdim);
    FMatrix CoB(F, selectinnerdim, FM.coldim());

    const size_t factoloops(1+(randomloops>>3));
    Factorizer(Alt, CoB, FM, factoloops, selectinnerdim, false);

    const size_t sa(density(Alt)), sc(density(CoB));

    std::ostringstream gout;
    auto gops(nbops), sops(nbops);
    size_t besti(std::numeric_limits<size_t>::max());

#pragma omp parallel for shared(Alt,CoB,sa,sc,gout,gops,sops,besti)
    for(size_t i=0; i<randomloops; ++i) {
	reseed('A', i);
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
	const bool empty(gout.tellp() == std::streampos(0));
	if ( empty || cmpOpCount(ABops,sops)
	     || (eqOpCount(ABops,sops) && (i < besti)) ) {
	    gout.clear(); gout.str(std::string());
	    gout << about.str();
	    sops = ABops; besti = i;
	    if (better) {
		if (verbose>0) {
		    std::clog << "# Found A: ("
			      << Alt.rowdim() << 'x' << Alt.coldim() << 'x'
			      << CoB.coldim() << ' ' << sa << '/' << sc
			      << ')' << '\t' << Bops.first << '+' << Aops.first
			      << '|' << Bops.second << '+' << Aops.second
			      << " instead of "
			      << gops.first << '|' << gops.second << "\t["
			      << i << '/' << omp_get_thread_num() << ']'
			      << std::endl;
		}
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
			   Givaro::Timer& global, const size_t randomloops, const int verbose) {
	// ============================================================
	// Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;

    Givaro::Timer chrono; chrono.start();
    std::ostringstream dout;
    auto dops(nbops), sops(nbops);
    size_t besti(std::numeric_limits<size_t>::max());

#pragma omp parallel for shared(M,T,dout,dops,sops,besti)
    for(size_t i=0; i<randomloops; ++i) {
	reseed('D', i);
	FMatrix lM(M, F);
	FMatrix lT(T, F);

	std::ostringstream ldout;
	    // Cancellation-free optimization
	input2Temps(ldout, lM.coldim(), 'i', 't', lT);
	auto lnbops( Optimizer(ldout, lM, 'o', 't', 'r') );

#pragma omp critical
	{
	const bool better(cmpOpCount(lnbops,dops));
	const bool empty(dout.tellp() == std::streampos(0));
	if ( empty || cmpOpCount(lnbops,sops)
	     || (eqOpCount(lnbops,sops) && (i < besti)) ) {
#ifdef VERBATIM_PARSING
	    std::clog << "# Found, direct: "
		      << lnbops.first << "\tadditions, "
		      << lnbops.second << "\tmultiplications." << std::endl;
#endif
	    dout.clear(); dout.str(std::string());
	    dout << ldout.str();
	    sops = lnbops; besti = i;
	    if (better) {
		if (verbose>0) {
		    std::clog << "# Found D: "
			      << lnbops.first << '|' << lnbops.second
			      << " instead of "
			      << dops.first << '|' << dops.second << "\t["
			      << i << '/' << omp_get_thread_num() << ']'
			      << std::endl;
		}
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
			      const Field& F, const Matrix& T, Givaro::Timer& global,
			      const size_t randomloops, const bool KFI, const int verbose) {
	// ============================================================
	// Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;
    Givaro::Timer chrono; chrono.start();
    std::ostringstream kout;
    auto kops(nbops), sops(nbops);
    size_t besti(std::numeric_limits<size_t>::max());
    bool nonzerokernel(true);

#pragma omp parallel for shared(T,kout,kops,sops,besti)
    for(size_t i=0; i<randomloops; ++i) {
	if (nonzerokernel) {
	    reseed('K', i);
	    FMatrix lT(T, F);
	    std::ostringstream lkout;
	    FMatrix NullSpace(F,lT.coldim(),lT.coldim());
	    auto lnbops( nullspacedecomp(lkout, NullSpace, lT, KFI) );

	    if (lnbops == Pair<size_t>{-1,-1}) {
#pragma omp critical
	    {
		nonzerokernel = false;
	    }
	    }

#pragma omp critical
	    {
	    const bool better(cmpOpCount(lnbops,kops));
	    const bool empty(kout.tellp() == std::streampos(0));
	    if ( empty || cmpOpCount(lnbops,sops)
		 || (eqOpCount(lnbops,sops) && (i < besti)) ) {
#ifdef VERBATIM_PARSING
		std::clog << "# Found, kernel: "
			  << lnbops.first << "\tadditions, "
			  << lnbops.second << "\tmultiplications." << std::endl;
#endif
		kout.clear(); kout.str(std::string());
		kout << lkout.str();
		sops = lnbops; besti = i;
		if (better) {
		    if (verbose>0) {
			std::clog << "# Found " << (KFI?'F':'K') << ": "
				  << lnbops.first << '|' << lnbops.second
				  << " instead of "
				  << kops.first << '|' << kops.second
				  << "\t[" << i << '/'
				  << omp_get_thread_num() << ']'
				  << std::endl;
		    }
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
			   Givaro::Timer& global, const size_t randomloops, const int verbose) {
	// ============================================================
	// Rebind matrix type over sub field matrix type
    typedef typename Matrix::template rebind<Field>::other FMatrix;
    Givaro::Timer chrono; chrono.start();
    const size_t m(T.coldim());
    std::vector<size_t> Fm { factorial(m) };
    std::ostringstream aout;
    auto aops(nbops), sops(nbops);
    size_t besti(std::numeric_limits<size_t>::max());

#pragma omp parallel for shared(Fm,T,aout,aops,sops,besti)
    for(size_t i=0; i<Fm.back(); ++i) {
	reseed('N', i);
	std::vector<size_t> l{kthpermutation(i,m,Fm)};
	FMatrix lT(T, F);
	std::ostringstream laout;
	FMatrix NullSpace(F,lT.coldim(),T.coldim());
	auto lkops( nullspacedecomp(laout, NullSpace, lT, l, 'o', 'i', mostCSE) );

#pragma omp critical
	{
	const bool better(cmpOpCount(lkops,aops));
	const bool empty(aout.tellp() == std::streampos(0));
	if ( empty || cmpOpCount(lkops,sops)
	     || (eqOpCount(lkops,sops) && (i < besti)) ) {
#ifdef VERBATIM_PARSING
	    std::clog << "# Found, kernel: "
		      << lkops.first << "\tadditions, "
		      << lkops.second << "\tmultiplications." << std::endl;
#endif
	    aout.clear(); aout.str(std::string());
	    aout << laout.str();
	    sops = lkops; besti = i;
	    if (better) {
		if (verbose>0) {
		    std::clog << "# Found E: "
			      << lkops.first << '|' << lkops.second
			      << " instead of "
			      << aops.first << '|' << aops.second << "\t["
			      << i << '/' << omp_get_thread_num() << ']'
			      << std::endl;
		}
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
			const bool tryLU, bool tryAB, const bool mostCSE,
			const bool allkernels, const bool KFI, const int verbose) {

    Pair<size_t> nbops(opsinit);

	// ============================================================
	// Optimize the Alternative factorization
    if (tryAB) {
//         for(int id(M.rowdim()-1); id >= M.coldim(); --id)
	for(size_t id(M.coldim()); id < M.rowdim(); ++id)
	    ABOptimiser(nbops, sout, F, M, T, id, global, randomloops,verbose);
    }


	// ============================================================
	// Otimize the whole matrix
    if (tryDirect) {
	CSEOptimiser(nbops, sout, F, M, T, global, randomloops,verbose);
    }

	// ============================================================
	// Separate independent and dependent rows
    if (tryKernel) {
	KernelOptimiser(nbops, sout, F, T, global, randomloops, KFI, verbose);
    }

	// ============================================================
	// Optimize the factorized matrix
    if (tryLU) {
	LUOptimiser(nbops, sout, F, M, T, global, randomloops,verbose);
    }

	// ============================================================
	// Exhaustive nullspace permutation search (if # is <= 12!)
    if (allkernels && (M.rowdim() < 13)) {
	AllKernelOpt(nbops, sout, F, T, mostCSE, global, randomloops,verbose);
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
