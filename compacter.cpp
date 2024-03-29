// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * Transposition of programs via Tellegen's principle
 * Program syntax: see main Tellegen function below
 * References:
 *   [ Tellegen's principle into practice.
 *     ISSAC'03:37-44, A. Bostan, G. Lecerf, É. Schost
 *     https://doi.org/10.1145/860854.860870 ]
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *     Strassen's algorithm is not optimally accurate
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
#include <streambuf>

// ============================================================
// Testing for a natural number within a string
bool isNatural(const std::string& s) {
    for(const auto& c: s) if (!std::isdigit(c)) return false;
    return true;
}
// ============================================================


// ============================================================
// Program is vector of lines
typedef std::vector<std::vector<std::string>> VProgram_t;

size_t progSize(const VProgram_t& P) {
    size_t PVs(0);
    for(const auto& line: P) PVs += line.size();
    return PVs-(P.size()<<1); // do not count ':=', nor ';'
}

std::ostream& printline(std::ostream& sout, const std::vector<std::string>& line) {
    sout << '#' << line.size() << ':';
    for(const auto& word: line)
        sout << word << ' ';
    return sout;
}

std::ostream& operator<<(std::ostream& sout, const VProgram_t& P) {
    for(const auto& line: P) {
        for(const auto& word: line)
            sout << word ;
        if (line.size()>0) sout << std::endl;
    }
    return sout;
}

// // ============================================================
// std::ostream& operator<<(std::ostream& sout, const std::vector<size_t>& v) {
//     sout << '[';
//     for(const auto& elt: v) sout << elt << ' ';
//     return sout << ']';
// }
// // ============================================================


// ============================================================
// Creating a vectror of lines from a text file
VProgram_t& programParser(VProgram_t& ProgramVector, std::stringstream& ssin) {

        // ============================================================
        // ====================
        // Working line by line
    std::string line;
    while(std::getline(ssin, line)) {
            // ========================================================
            // Finding the LHS
        const std::size_t postassign(line.find(":=",0));
        if (postassign != std::string::npos) {
            std::vector<std::string> wordVector;

            wordVector.push_back(line.substr(0,postassign));

            wordVector.push_back(":=");

                // ====================================================
                // Now the RHS with successive monomials
                //     until end of string (npos)
            size_t prev(postassign+2);
            for(size_t pos(0);
                (pos = line.find_first_of("+-*/;", prev)) != std::string::npos;
                prev = pos+1) {
                    // One of '+', '-', '*', '/', ';'
                std::string delimiter(line.substr(pos, 1));
                if (pos > prev) {
                        // delimited word
                    std::string node(line.substr(prev, pos-prev));

                    if ((delimiter == "/") && isNatural(node)) {
                            // Found rational number, get next natural
                        prev = pos+1;
                        pos = line.find_first_of("+-*;", prev);
                        const std::string nnod(line.substr(prev, pos-prev));
#ifdef VERBATIM_PARSING
                        std::clog << "# Rational found: "
                                  << node << '/' << nnod << std::endl;
#endif
                        node.append(delimiter);
                        node.append(nnod);
                        delimiter = line.substr(pos, 1);
                    }
                    wordVector.push_back(node);		// add found word

                }
                wordVector.push_back(delimiter);	// add delimiter
            }

                // =================
                // Add the rest of the line as is
                // will be ignored when after ';'
            if (prev < line.length())
                wordVector.push_back(line.substr(prev, std::string::npos));

#ifdef VERBATIM_PARSING
                // =================
                // Got a line
            std::clog << "# Parsed: ";
            for(const auto& iter: wordVector) std::clog << iter << ' ';
            std::clog << std::endl;
#endif
            ProgramVector.push_back(wordVector);

        }
    }

    return ProgramVector;
}
// ============================================================






// ============================================================
// Test for no-op operations like a0 := r0 ;
auto idempots {[](const std::vector<std::string>& s) {
    return ((s.size()==4) && (s[0]==s[2]));} };

// Comparing two vectors by their first value within a pair
typedef std::pair<std::string, std::vector<size_t>> pairSVs_t;
auto prgorder {[](const pairSVs_t& a, const pairSVs_t& b) {
    return a.second.front()<b.second.front();} };
// ============================================================



// ============================================================
// Rewrites:
//   [1] temporarily replaces output variables by temporary, except last ones
//   [2] replaces variable only assigned to temporary, directly by input
//   [*] Removes all no-op operations like ai:=ai;
//   [3] Backward reassingment of output variables
//   [4] (optional) rewrites singly used variables in-place
size_t variablesThiner(VProgram_t& P, const bool simplSingle=true,
                       const char inchar = 'i', const char outchar = 'o') {
        // ==================================
        // Find two unused variable name
    std::set<char> varsChar;
    for(const auto& line: P) {
        for(const auto& word: line) {
            varsChar.insert(word[0]);
        }
    }
    if (varsChar.size() > 50) {
        std::cerr << "\033[1;31m**** ERROR ****\033[0m"
                  << "not enough free single char variables." << std::endl;
        exit(-1);
    }
    char freechar('a');
    for( ; varsChar.find(freechar) != varsChar.end(); ++freechar){}
    if (freechar > 'z') {
        for(freechar='A'; varsChar.find(freechar)!=varsChar.end();++freechar){}
    }
    char tmpchar(freechar);
    for(++tmpchar ; varsChar.find(tmpchar) != varsChar.end(); ++tmpchar){}
    if (tmpchar > 'z') {
        for(tmpchar='A'; varsChar.find(tmpchar)!=varsChar.end();++tmpchar){}
    }


        // ==================================
        // [1] Outputs only at the end
    VProgram_t OutP;
    for(auto line(P.begin()); line != P.end(); ++line) {
        std::string outvar(line->front());
        if (outvar[0] == outchar) {
            std::string repvar(outvar); repvar[0] = freechar;
            for(auto next(line); next != P.end(); ++next) {
                for(auto& word: *next) {
                    if (word == outvar) word = repvar;
                }
            }
            OutP.emplace_back( std::vector{outvar, std::string(":="), repvar, std::string(";")} );
        }
    }

    P.insert(P.end(), OutP.begin(), OutP.end());

//     std::clog << std::string(20,'#') << " AFT outputs" << std::endl;
//     std::clog << P << std::endl;


        // ==================================
        // [2] Direct substitution of affectations
    for(auto line(P.begin()); line != P.end(); ++line) {
        if ((line->size() == 4) && (line->front()[0] != outchar) ) {
            std::string& outvar(line->front()), &invar(*(line->begin()+2));
            auto next(line);
            for(++next; next != P.end(); ++next) {
                    // pass output and ':='
                for(auto word(next->begin()+ 2); word != next->end(); ++word) {
                    if (*word == outvar) *word = invar;
                }
                if (next->front() == outvar) break;
            }
            outvar = invar;
        }
    }

        // ==================================
        // [*] Removes all idempotents
    P.erase(std::remove_if(P.begin(), P.end(), idempots), P.end());

        // ==================================
        // [3] Backward substitution of outputs
    for(auto line(P.rbegin()); line != P.rend(); ++line) {
        if ((line->size() == 4) && (line->front()[0] != inchar) ) {

            std::string& outvar(line->front()), &invar(*(line->begin()+2));
            for(auto next(line+1); next != P.rend(); ++next) {
                    // if variable is set, stop substituting
                if (next->front() == invar) {
                    next->front() = outvar;
                    invar = outvar;
                    break;
                }
                    // substitute in next line
                for(auto word(next->begin()+ 2); word != next->end(); ++word) {
                    if (*word == invar) *word = outvar;
                }
            }
        }
    }

        // ==================================
        // [*] Removes all idempotents
    P.erase(std::remove_if(P.begin(), P.end(), idempots), P.end());


    if (! simplSingle) return 0;

        // ==================================
        // [4] Rewriting singly used variables
    size_t tmpnum(0), merged(0);

        // Occurences of temporary variables
    std::map<std::string, std::vector<size_t> > varsSet, varsUse;
    for(size_t i=0; i<P.size(); ++i) {
        auto& line(P[i]);
        auto& variable(line.front());
        const auto prevar(line.front());
        if ( (prevar[0] != outchar)
             && (varsSet.find(variable) != varsSet.end()) ) {
                // Variable is overwritten, change it
            variable = tmpchar;
            variable += std::to_string(tmpnum);

            for(size_t j=i+1; j<P.size(); ++j) {
                auto& nextline(P[j]);
                for(auto& word:nextline) {
                    if (word == prevar) {
                        word = variable;
                    }
                }
            }

            ++tmpnum;
        }
        varsSet[variable].emplace_back(i);
        for(const auto& word: line) {
            if (varsSet.find(word) != varsSet.end()) {
                varsUse[word].emplace_back(i);
            }
        }
    }

        // Need to sort map by (front of) values (first occurence of that var)
        // Thus rewriting is done in the program order
    std::vector<std::pair<std::string, std::vector<size_t>>> ordsUse;
    for (const auto& [variable, occurences] : varsUse)
        ordsUse.emplace_back(variable,occurences);
    std::sort(ordsUse.begin(), ordsUse.end(), prgorder);


        // finding singly used variables and replacing them
    for (const auto& [variable, occurences] : ordsUse) {
        if ( (variable[0] != outchar) && (occurences.size()==2) ) {
            const size_t i(occurences.front()), j(occurences.back());
            auto& init(P[i]); auto& line(P[j]);
// printline(std::clog << "# P[" << i << "]:", init) << "  -->  ";
// printline(std::clog << "# P[" << j << "]:", line) << std::endl;

                // Single use of variable once set
            auto word(std::find(line.begin(), line.end(), variable));
            *word = '(';
                // start of replacement expression, from 2 to end-1
            const auto penultimate(init.end()-1);
            for(auto iter(init.begin()+2); iter != penultimate; ++iter)
                *word += *iter;
            *word += ')';
            merged += (init.size()-3);
            P[i].resize(0); // No need for that variable anymore
        }
    }

    return merged;
}
// ============================================================





// ============================================================
// Main compacting procedure, parsing then rewriting
std::ostream& Compacter(std::ostream& sout, std::istream& input,
                        const bool simplSingle=true) {
        // Files contains a program with the following SYNTAX
        // [+] input variables start with a character (default is 'i')
        // [+] output variables start with a character (default is 'o')
        // [+] lines are of the forms:
        //    xi := sum (yi op(li)), with sum: + or -, and: op * or / or empty
        //    One operation per line (';' stops the parsing of that line)
    std::stringstream ssin; ssin << input.rdbuf();

        // Line by line parsing
    VProgram_t ProgramVector; programParser(ProgramVector, ssin);
    const size_t PVs { progSize(ProgramVector) };
    std::clog << std::string(40,'#') << std::endl;

        // Semantic line removal
    size_t merged = variablesThiner(ProgramVector, simplSingle);
    sout << ProgramVector;
    std::clog << std::string(40,'#') << std::endl;

    const size_t PRs { merged + progSize(ProgramVector) };

        // Comparing number of elements in the programs
    std::clog << "# \033[1;32m" << PRs << "\telements\tinstead of "
              << PVs << "\033[0m" << std::endl;
    std::clog << std::string(40,'#') << std::endl;

    return sout;
}
// ============================================================




// ============================================================
// Main: select between file / std::cin
int main(int argc, char** argv) {
    bool simplSingle(false);
    std::string filename;

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0] << "[-s] [stdin|file.prg]\n"
                      << "  -s: also replace singly used variables\n";
            exit(-1);
        }
        else if (args == "-s") { simplSingle = true; }
        else { filename = args; }
    }

    if (filename == "") {
        Compacter(std::cout, std::cin, simplSingle);
    } else {
        std::ifstream ifile(filename);
        if ( ifile ) {
            Compacter(std::cout, ifile, simplSingle);
            ifile.close();
        }
    }

    return 0;
}
// ============================================================
