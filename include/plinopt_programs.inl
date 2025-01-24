// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet
// ==========================================================================

/****************************************************************
 * PLinOpt Library, program manipulations
 * Reference:
 *   [ J-G. Dumas, B. Grenet; Jul. 2023
 *     In-place accumulation of fast multiplication formulae
 *     (https://hal.science/hal-04167499) ]
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *     Strassen's algorithm is not optimally accurate
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/

#include "plinopt_library.h"

#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <numeric>
#include <streambuf>


// ============================================================
// Testing for a natural number within a string
inline bool isNatural(const std::string& s) {
    for(const auto& c: s) if (!std::isdigit(c)) return false;
    return true;
}
// ============================================================

// ============================================================
inline size_t progSize(const VProgram_t& P) {
    size_t PVs(0);
    for(const auto& line: P) PVs += line.size();
    return PVs-(P.size()<<1); // do not count ':=', nor ';'
}
// ============================================================


// ============================================================
inline std::ostream& operator<<(std::ostream& sout, const VProgram_t& P) {
    for(const auto& line: P) {
        for(const auto& word: line)
            sout << word ;
        if (line.size()>0) sout << std::endl;
    }
    return sout;
}
// ============================================================


// ============================================================
std::ostream& operator<<(std::ostream& out,
                         const std::map<std::string,size_t>& m) {
    out << '{';
    for(const auto& [variable, index] : m)
        out << variable << '(' << index << ')' << ',';
    return out << '}';
}
// ============================================================


// ============================================================
inline std::ostream& printline(std::ostream& sout,
                               const std::vector<std::string>& line) {
    sout << '#' << line.size() << ':';
    for(const auto& word: line)
        sout << word << ' ';
    return sout;
}
// ===============================================================



// ============================================================

// Test for no-op operations like a0 := r0 ;
auto idempots {[](const std::vector<std::string>& s) {
    return ((s.size()==4) && (s[0]==s[2]));} };

// Test for no-op operations like a0 := r0 ;
auto emptyline {[](const std::vector<std::string>& s) {
    return (s.size()==0);} };

// Comparing two vectors by their first (front) value, as second in a pair
typedef std::pair<std::string, std::vector<size_t>> pairSVs_t;
auto prgorder {[](const pairSVs_t& a, const pairSVs_t& b) {
    return a.second.front()<b.second.front();} };

// Test addition or substraction
auto isAddSub {[](const std::string& s) {
    return (s=="+") || (s=="-");} };

auto isMulDiv {[](const std::string& s) {
    return (s=="*") || (s=="/");} };

auto isParAff {[](const std::string& s) {
    return (s=="(") || (s==":=");} };

auto isVariable {[](const std::string& s) {
    return (s.find_first_of("+-*/;:=()") == std::string::npos);} };

// Echange '+' with '-' and vice-versa
void swapsign(std::string& s) {
    if (s=="+") { s.replace(0,1,1,'-'); }
    else if (s=="-") { s.replace(0,1,1,'+'); }
}
// ============================================================




// ============================================================
// Stack of lines (strings), to print in reverse order
struct stackbuf : public std::streambuf {
    stackbuf(): d_lines(1) {}
    std::string str() const {
        return std::accumulate(this->d_lines.rbegin(),
                               this->d_lines.rend(),
                               std::string());
    }

    auto rbegin() const { return this->d_lines.rbegin(); }
    auto rend() const { return this->d_lines.rend(); }

private:
    int overflow(int c) {
        if (c != std::char_traits<char>::eof()) {
            this->d_lines.back().push_back(c);
            if (c == '\n') {
                this->d_lines.push_back(std::string());
            }
        }
        return std::char_traits<char>::not_eof(c);
    }
    std::vector<std::string> d_lines;
};
// ============================================================



// ============================================================
// Find unused variable name
char unusedChar(const std::set<char>& C, const char cstart /* = 'a'-1 */) {
    if (C.size() > 50) {
        std::cerr << "\033[1;31m**** ERROR ****\033[0m"
                  << "not enough free single char variables." << std::endl;
        exit(-1);
    }
    char tmpchar(cstart);
    for(; C.find(++tmpchar) != C.end();){}
    if (tmpchar > 'z') {
        for(tmpchar='A'; C.find(tmpchar)!=C.end();++tmpchar){}
    }
    return tmpchar;
}
// ============================================================


// ============================================================
// Main parsing procedure, writing the ouput program
int Tellegen(std::istream& input,
             const char ichar /* = 'i'*/, const char ochar /*= 'o'*/,
             const char cchar /* = 'c'*/, const char zchar /* = 'z'*/) {
        // Files contains a program with the following SYNTAX
        // [+] constants name start with 'c'
        // [+] input variables start with 'i' --> transformed into 'o'
		//    supposed to be just a simple affectation "x0:=i0;"
        // [+] output variables start with 'o' --> transformed into 'i'
        // [+] lines are of the forms:
        //    xi := sum (yi op(li)), with sum: + or -, and: op * or / or empty
        //        xi, yi are supposed to be without any of "+-:=/*;"
        //        liin particular rational numbers a/b should be precomputed,
        //        e.g., into a constant 'cxx:=a/b;' when used in computation.
        //        e.g. replace: "x0:=x1*3/4;" by "c0:=3/4;x0:=x1*c0;"
        //    each yi occurs singly in a given line
        //    One operation per line (';' stops the parsing of that line)
    std::stringstream ssin; ssin << input.rdbuf();
    stackbuf sbuf;
    std::ostream sout(&sbuf);

        // Sets of variables
    std::set<std::string> varSet;   // List of found variables
    VProgram_t outSet;              // Final program with outputs
    std::set<std::string> inSet;    // List of input variables

    std::pair<size_t,size_t> Nops(0,0);

    std::clog << std::string(40,'#') << std::endl;

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

                // Constants are left unmodified
            if (wordVector.front()[0] == cchar) {
#ifdef VERBATIM_PARSING
                std::clog << "# Constant found, unmodified : "
                          << line << std::endl;
#endif
                std::cout << line << std::endl;
                continue;
            }
                // Outputs are transposed into inputs
            if (wordVector.front()[0] == ochar) {
                std::string oneout(wordVector.front());
                inSet.insert( oneout );
                oneout[0] = ichar;
#ifdef VERBATIM_PARSING
                std::clog << "# Output found, " << wordVector.front()
                          << ", becomes input: " << oneout<< std::endl;
#endif
            }

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

                // =================
                // Transposing that line
            const std::string fixvar(wordVector.front());
            for(auto iter(wordVector.begin()+2); iter != wordVector.end(); ) {
                    // Stop at ';'
                if (*iter == ";") break;

                    // Sign of monomial
                bool neg(false);
                if (*iter == "-") { ++iter; neg=true; }
                if (*iter == "+") { ++iter; }

                    // input variables are transposed into output
                const std::string variable(*iter); ++iter;
                if (variable == "0") continue;

                if (variable[0] == ichar) {
                    std::string onein(fixvar);
                    onein[0] = ochar;
                    outSet.push_back( { onein, ":=", fixvar, ";" });
#ifdef VERBATIM_PARSING
                    std::clog << "# Input found, " << variable
                              << ", becomes output: " << onein << std::endl;
#endif
                    continue;
                }

                    // Monomial indeterminate
                varSet.insert(variable);

                    // Monomial coefficient
                std::string multiplier;
                if ((*iter == "/") || (*iter == "*")) {
                    ++(Nops.second);
                    multiplier = *iter;
                    ++iter;
                    multiplier.append(*iter);
                    ++iter;
                    if (*iter == "/") {
                            // Rational found, get next natural
                        multiplier.append(*iter);
                        ++iter;
                        multiplier.append(*iter);
                    }
                } else {
                    multiplier = "";
                }

                    // Transposed monomial into the stack
                sout << variable << ":=" << variable << (neg?'-':'+') << fixvar << multiplier << ';' << '\n';
            }
        }
    }
        // End line by line parsing
        // ============================================================

    std::clog << std::string(40,'#') << std::endl;

        // ============================================================
        // Select only temporary variables (others are initialized)
    std::vector<std::string> varVector;
    for(const auto& iter: varSet) {
        if (inSet.find(iter) == inSet.end()) varVector.push_back(iter);
    }

        // ============================================================
        // Now Produce program reversed

        // ======================================
        // Outputs each line in reverse order
        // (and simplifies first accumulation
        //   instead of: "xi:=ZZ; xi:=xi + ... ;"
        //       directly write: "xi:=ZZ + ...;")

    std::set<std::string> modSet; // first modification of variable

    for(auto iter=sbuf.rbegin(); iter != sbuf.rend(); ++iter) {
        std::string line(*iter);
        if (line.size() == 0) continue;

            // ==================================
            // Instead of "xi:=0; xi:=xi + ... ;"
            //          directly write "xi:=...;"
        const auto varlen = line.find(":=",0);
        const std::string variable(line.substr(0,varlen));

            // 'o' is replaced by 'z' in computations
        if (line[0] == ochar) line[0]=zchar;

        const auto varposvec =
            std::find(varVector.begin(), varVector.end(), variable);

        if ( varposvec != varVector.end()) {
#ifdef VERBATIM_PARSING
            std::clog << "# First temporary, accumulation with " << variable
                      << "=0 is simplified in: " << line << std::endl;
#endif
            const auto accupos = line.find(variable, varlen);

#if DEBUG
            if (accupos != varlen+2) {
                std::cerr << "\033[1;31m**** ERROR ****\033[0m"
                          << " first occurence of temporary " << variable
                          << " is not an accumulation." << std::endl;
                exit(-1);
            }
#endif

            line.erase(accupos, varlen);
            if (line[accupos] == '+') {
                    // replace 0+x by x, but 0-x by -x
                line.erase(accupos, 1);
            }
            varVector.erase(varposvec);
        }

            // ==================================
            // Instead of "o1:=i1; t1:=o1 + ...;"
            //     directly write "t1:=i1 + ...;"
            //     parsing the RHS with successive monomials
        const auto lhslen = line.find(":=",0);
        for(size_t prev(lhslen+2), pos(0);
            (pos = line.find_first_of("+-*/;", prev)) != std::string::npos;
            prev = pos+1) {
            if (line[prev] == ochar) {
                if (modSet.find(line.substr(prev, pos-prev)) == modSet.end()){
#ifdef VERBATIM_PARSING
                    std::clog << "# Unmodified input usage of: "
                              << line.substr(prev, pos-prev)
                              << ", is simplified in RHS: " << line
                              << std::endl;
#endif
                    line[prev]=ichar;
                } else {
                    line[prev]=zchar; // 'o' is replaced by 'z' in computations
                }
            }
        }

            // ==================================
            // Now variable has been modified
            // --> will not be simplified anymore
        modSet.insert(variable);

            // ==================================
            // Just output the line

        if (line.find_first_of("+-", 0) != std::string::npos)
            ++(Nops.first);
        if (line.find("=-", 0) != std::string::npos)
            --(Nops.first);

        std::cout << line << std::flush;
    }

#if DEBUG
    if (varVector.size() > 0) {
        std::cerr << "\033[1;31m**** ERROR ****\033[0m"
                  << " First accumulation of temporaries failed: ";
        for(const auto& iter:varVector) std::cerr << iter << ' ';
        std::cerr << std::endl;
    }

    if (varSet != modSet) {
        std::cerr << "\033[1;31m**** ERROR ****\033[0m"
                  << " found varSet differ from modSet : \n"
                  << varSet << std::endl
                  << modSet << std::endl;
    }
#endif

    std::clog << std::string(40,'#') << std::endl;

        // Produce output results "oi:=xi;"
    for(const auto& iter: outSet) {
#if DEBUG
        if ((iter[1] != ":=") || (iter[3] != ";")) {
            std::cerr << "\033[1;31m**** ERROR ****\033[0m"
                      << " malformed output :" << iter << std::endl;
        }
#endif
        std::cout << iter.front() << ":="
                  << (modSet.find(iter[2]) == modSet.end() ? "0" : iter[2])
                  << ';' << std::endl;
    }


    const int dimOffset(inSet.size()-outSet.size());


    std::clog << std::string(40,'#') << std::endl;
    std::clog << "# \033[1;32m" << Nops.first << "\tadditions ("
              << (Nops.first-dimOffset)
              << (dimOffset<0?'-':'+') << abs(dimOffset)
              << ")\033[0m" << std::endl;
    std::clog << "# \033[1;32m" << Nops.second << "\tmultiplications"
              << "\033[0m" << std::endl;
    std::clog << std::string(40,'#') << std::endl;

    return 0;
}
// ============================================================



// ============================================================
// Creating a vector of lines from a text file program
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
                (pos = line.find_first_of(")(+-*/;", prev)) != std::string::npos;
                prev = pos+1) {
                    // One of '+', '-', '*', '/', ';'
                std::string delimiter(line.substr(pos, 1));
                if (pos > prev) {
                        // delimited word
                    std::string node(line.substr(prev, pos-prev));

                    if ((delimiter == "/") && isNatural(node)) {
                            // Found rational number, get next natural
                        prev = pos+1;
                        pos = line.find_first_of("()+-*;", prev);
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
// Rotates lines starting with a '-'
bool rotateMinus(std::vector<std::string>& line) {
    if (line[2] == "-") {
        size_t depth(0);
        for(size_t i=2; i<line.size(); ++i) {
            if (line[i] == "(") ++depth;
            if (line[i] == ")") --depth;
            if ((depth ==0) && (line[i] == "+")) {
                std::rotate(line.begin()+2,line.begin()+i,line.end()-1);
                line.erase(line.begin()+2); // no need for '+' anymore
                return true;
            }
        }
    }
    return false;
}
// ============================================================


// ============================================================
// Negating line after afect
std::vector<std::string> negateLine(const std::vector<std::string>& line) {
    size_t depth(0);
    std::vector<std::string> varline(line.size()+1); varline.resize(2);
    varline[0] = line[0];
    varline[1] = line[1];
    size_t j(2);
    if (line[j] != "-") {
        varline.push_back("-");
        ++j;
    }
    varline.insert(varline.end(), line.begin()+2, line.end());
    for(; j<varline.size(); ++j) {
        if (varline[j] == "(") ++depth;
        if (varline[j] == ")") --depth;
        if (depth == 0) {
            if (varline[j] == "-") {
                if (isParAff(varline[j-1])) {
                    std::rotate(varline.begin()+j,varline.begin()+j+1,
                                varline.end());
                    varline.pop_back();
                    --j;
                } else {
                    varline[j]= "+";
                }
            } else if (varline[j] == "+") {
                varline[j]= "-";
            }
        }
    }
// std::clog << "# Negate line: " << line << " ---> " << varline << std::endl;
    return std::move(varline);
}
// ============================================================


// ============================================================
// Negating variable in line
bool negatingVariable(std::vector<std::string>& line,
                      const std::vector<std::string>& oldline,
                      const std::string& variable) {
    bool linmod(false);
    line.assign(oldline.begin(),oldline.end());
    for(size_t j=2; j<line.size(); ++j) {
        if (line[j] == variable) {
            if (line[j-1] == "+") {
                line[j-1] = "-";
            } else if (line[j-1] == "-") {
                if (isParAff(line[j-2])) {
                    std::rotate(line.begin()+j-1,line.begin()+j,
                                line.end());
                    line.pop_back();
                } else {
                    line[j-1] = "+";
                }
            } else if (isParAff(line[j-1])) {
                line.push_back("-");
                std::rotate(line.begin()+j,line.end()-1,line.end());
                ++j;
            }
            linmod = true;
        }
    }
// std::clog << "# Negate variable (" << variable << "); " << oldline << " ---> " << line << std::endl;
    return linmod;
}
// ============================================================


// ============================================================
// Negate temporaries for output lines starting with a '-'
size_t endingMinus(VProgram_t & P, const char inchar, const char outchar) {
    if (P.empty()) return 0;
    for(int i=P.size()-1; i>=0; --i) {
        auto & varline(P[i]);
        if (varline.front()[0] == outchar) {
            rotateMinus(varline);
            if ( varline[2] == "-") {
                for(size_t j=3; j<varline.size(); ++j) {
                    const auto variable(varline[j]);
                    if (isVariable(variable) &&
                        (variable[0] != inchar) && (variable[0] != outchar) ) {
// std::clog << "# found variable: " << variable << std::endl;
                        for(int k=i-1; k>=0; --k) {
                            if (P[k].front() == variable) {
// std::clog << "# found line[" << k << "]: " << P[k] << std::endl;
                                auto refline(negateLine(P[k]));
                                rotateMinus(refline);
                                P[k] = std::move(refline);
// std::clog << "# repl. line[" << k << "]: " << P[k] << std::endl;

                                for(size_t l=k+1; l<P.size(); ++l) {
                                    std::vector<std::string> negline;
                                    bool linmod(negatingVariable(negline, P[l], variable));
                                    rotateMinus(negline);
                                    P[l] = std::move(negline);
// std::clog << "# impa. line[" << l << "]: " << P[l] << std::endl;
                                    if (P[l].front() == variable) break;
                                }
                            }
                        }
                        break;
                    }
                }
            }
        }
    }
    return 0;
}
// ============================================================





// ============================================================
size_t countMinus(const VProgram_t& P) {
    size_t cm(0);
    for(auto& line : P) if (line[2] == "-") ++cm;
    return cm;
}
// ============================================================


// ============================================================
// Negates some temporary vairables if it reduces the leading minuses
size_t swapMinus(VProgram_t & P, const char outchar) {
    for(auto& line : P) rotateMinus(line);
    size_t pm( countMinus(P) );
    if (pm == 0) return 0;


    for(size_t i=0; i<P.size(); ++i) {
        if ( (P[i])[2] == "-" && (P[i].front()[0] != outchar) ) {
            std::map<size_t, std::vector<std::string>> mP;
                // Negating line
            mP[i] = negateLine(P[i]);

                // Negating variable until reaffect
            const auto& variable(mP[i][0]);
            for( ; ++i < P.size() ; ) {
                std::vector<std::string> line; line.reserve(P[i].size());
                bool linmod(negatingVariable(line, P[i], variable));
                if (linmod) {
                    rotateMinus(line);
                    mP[i] = line;
                }
                if (P[i].front() == variable) break;
            }

            size_t cm(pm);
            for (const auto& [k, line] : mP) {
                if ( (P[k][2] == "-") && (mP[k][2] != "-") ) --cm;
                if ( (mP[k][2] == "-") && (P[k][2] != "-") ) ++cm;
            }
            if (cm < pm) {
#ifdef VERBATIM_PARSING
                std::clog << "# Reduced leading minus from " << pm
                          << " to " << cm << std::endl;
#endif
                for (const auto& [k, line] : mP) P[k] = std::move(line);
                pm = cm;
            }
        }
    }
    return pm;
}
// ============================================================




// ============================================================
// Rewrites:
//   [1] temporarily replaces output variables by temporary, except last ones
//   [2] replaces variable only assigned to temporary, directly by input
//   [*] Removes all no-op operations like ai:=ai;
//   [3] Backward reassingment of output variables
//   [4] (optional) rewrites singly used variables in-place
size_t variablesTrimer(VProgram_t& P, const bool simplSingle,
                       const char inchar, const char outchar) {
        // ==================================
        // Find two unused variable names
    std::set<char> varsChar;
    for(const auto& line: P) for(const auto& word: line)
        varsChar.insert(word[0]);

    char freechar(unusedChar(varsChar));
    char tmpchar(unusedChar(varsChar,freechar));

        // ==================================
        // [1] Use output variables only at the end of the program
        //     --> replace them temporarily with another name (freechar)
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

        // ==================================
        // [2] Direct substitution of simple affectations "x := a ;"
        //     --> variable x is replaced by variable a thereafter
    for(auto line(P.begin()); line != P.end(); ++line) {
        if ((line->size() == 4) && (line->front()[0] != outchar)) {
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
        // [*] Removes all idempotents "x := x ;"
    P.erase(std::remove_if(P.begin(), P.end(), idempots), P.end());

        // ==================================
        // [3] Backward substitution of outputs
        //     --> from the end of the program, look for output variables
    for(auto line(P.rbegin()); line != P.rend(); ++line) {
        std::string& outvar(line->front()), &invar(*(line->begin()+2));
        if ((line->size() == 4) && (invar != "0") &&
            (invar[0] != inchar) && (invar[0] != outchar) ) {
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
        // [*] Removes all idempotents "x := x ;"
    P.erase(std::remove_if(P.begin(), P.end(), idempots), P.end());

        // ==================================
        // [*] Rotates lines starting with a '-'
    for(auto& line : P) rotateMinus(line);

        // ==================================
        // [*] Negate temporary variables so that output do not start with a '-'
    endingMinus(P, inchar, outchar);

    if (! simplSingle) return 0;

        // ==================================
        // [4] Rewriting singly used variables
        //     --> "x:= ...;" folllowed by "y:=... x ...;";
    size_t tmpnum(0);

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
            variable += std::to_string(++tmpnum);

            for(size_t j=i+1; j<P.size(); ++j) {
                auto& nextline(P[j]);
                for(auto& word:nextline) {
                    if (word == prevar) {
                        word = variable;
                    }
                }
            }
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
            rotateMinus(init); // less sign changes if possible

#ifdef VERBATIM_PARSING
            std::clog << "# variable: " << variable
                      << ", P[" << i << "]: " << init
                      << ", P[" << j << "]: " << line << std::endl;
#endif

                // Single use of variable once set, not as storage
            size_t varloc(2);
            for(; varloc<line.size(); ++varloc) {
                if (line[varloc] == variable) break;
            }
            if (varloc<line.size()) {
                    // Whether replacement requires sign-change
                bool changesign(false), multimonomial(false);
                if ( (init[2] == "-") && (isAddSub(line[varloc-1])) ) {
                    swapsign(line[varloc-1]);   // Change dest sign
                    init.erase(init.begin()+2);	// remove src '-'
                    changesign = true;
                }

                    // Whether replacement requires parenthesis
                auto penultimate(init.end()-1);
                size_t depth(0);
                for(auto iter(init.begin()+3); iter != penultimate; ++iter) {
                    if (*iter == "(") ++depth;
                    else if (*iter == ")") --depth;
                    else if (isAddSub(*iter)) {
                        if (changesign && (depth ==0)) swapsign(*iter);
                        multimonomial = true;
                    }
                }

                    // Replacement expression
                auto replength(init.end()-init.begin()-2);
                bool tobeneg(line[varloc-1] == "-");
                const bool tobemul(isMulDiv(line[varloc+1]));

                if ((!tobemul) && tobeneg && multimonomial) {
                        // Parenthesis not needed
                        // Just negate the variable in place
                    init = negateLine(init);
                    rotateMinus(init);
                    replength = init.end()-init.begin()-2;
                    penultimate = init.end()-1;
                    if ( (varloc == 3) || (init[2] == "-") ) {
                            // If starting cmd, remove the minus sign
                            // Or if replacement starts also with minus
                        --varloc;
                        line.erase (line.begin()+varloc);
                    } else {
                            // Otherwise transform it into a plus sign
                        line[varloc-1] = "+";
                    }
                    tobeneg = false;
                }

                line[varloc] = std::move(init[2]); // replace previous variable
                line.insert(line.begin()+varloc+1,init.begin()+3,penultimate);
                if (multimonomial && (tobemul || tobeneg)) {
                    line.insert(line.begin()+varloc,"(");
                    line.insert(line.begin()+varloc+replength,")");
                }

                P[i].resize(0); // No need for that variable (& line) anymore
            }
#ifdef VERBATIM_PARSING
            std::clog << "#         : " << variable
                      << " --> P[" << j << "]: " << line << std::endl;
#endif

        }
    }

        // ==================================
        // [*] Removes empty lines
    P.erase(std::remove_if(P.begin(), P.end(), emptyline), P.end());

        // ==================================
        // [*] Tries to reduce lines starting with a '-'
        //     by negating some temporary variables
    swapMinus(P, outchar);

        // ==================================
        // [*] Negate temporary variables so that output do not start with a '-'
    endingMinus(P, inchar, outchar);

        // ==================================
        // [*] Rotates lines starting with a '-'
    for(auto& line : P) rotateMinus(line);

    return tmpnum;
}
// ============================================================



// ============================================================
template<typename _Mat>
_Mat& matrixBuilder(_Mat& A, const VProgram_t& P, const char outchar /* ='o'*/) {
    using Field=typename _Mat::Field;
    const Field& F(A.field());
    _Mat M(F);
    std::map<std::string,size_t> inputs;
    std::map<std::string,size_t> variables;

    for(const auto& line: P) {
#ifdef VERBATIM_PARSING
        printline(std::clog << "# line: ", line) << std::endl;
#endif
        const auto& output(line.front());

        if (variables.find(output) == variables.end()) {
            const size_t m(M.rowdim()), n(M.coldim());
            variables[output] = m;
            M.resize(m+1,n);
        }
        const size_t i(variables[output]);
        for(auto word=line.begin()+2; word !=line.end(); ++word) {
// std::clog << "## word: " << *word << " (output:" << output << ')' << std::endl;
            typename Field::Element coeff; F.init(coeff); F.assign(coeff, F.one);
            if (*word == ";") break;
            if (*word == "-") {
                F.mulin(coeff, F.mOne); ++word;
            } else if (*word == "+") {
                ++word;
            }
            if (*word == output) {
                ++word;
                typename Field::Element tmp; F.init(tmp);
                if (*word == "*") {
                    ++word;
                    std::stringstream ssin; ssin << word->c_str();
                    F.read(ssin, tmp);
                    F.mulin(coeff,tmp);
                } else if (*word == "/") {
                    ++word;
                    std::stringstream ssin; ssin << word->c_str();
                    F.read(ssin, tmp);
                    F.divin(coeff, tmp);
                } else {
                    --word;
                }
                if (F.isMOne(coeff)) negRow(M, i);
                else if (! F.isOne(coeff)) mulRow(M, i, coeff);
// std::clog << "## mul by " << coeff << std::endl;
                continue;
            }
            if (std::find_if(word->begin(),word->end(), [](const char&c) { return std::isalpha(c); } ) == word->end()) {
                    // Constant found: should be zero
                assert(std::stoi(*word)==0);
                continue;
            }
            if (variables.find(*word) == variables.end()) {
                    // Input found
                if (inputs.find(*word) == inputs.end()) {
                    const size_t m(M.rowdim()), n(M.coldim());
                    inputs[*word] = n;
                    M.resize(m,n+1);
                }
                const size_t j = inputs[*word]; ++word;
                typename Field::Element tmp; F.init(tmp);
                if (*word == "*") {
                    ++word;
                    std::stringstream ssin; ssin << word->c_str();
                    F.read(ssin, tmp);
                    F.mulin(coeff,tmp);
                } else if (*word == "/") {
                    ++word;
                    std::stringstream ssin; ssin << word->c_str();
                    F.read(ssin, tmp);
                    F.divin(coeff, tmp);
                }  else {
                    --word;
                }
                F.assign(tmp, M.getEntry(i,j)); // might be zero (no refEntry)
                M.setEntry(i,j, F.addin(tmp, coeff) );
            } else {
                const size_t j(variables[*word]);
// std::clog << "### row[" << i << "] : " << output << std::endl;
// std::clog << "### & row[" << j << "] : " << *word << std::endl;

                ++word;

                if (*word == "*") {
                    ++word;
                    typename Field::Element tmp;
                    std::stringstream ssin; ssin << word->c_str();
                    F.read(ssin, tmp);
                    F.mulin(coeff,tmp);
// std::clog << "### mul : " << coeff << std::endl;
                } else if (*word == "/") {
                    ++word;
                    typename Field::Element tmp;
                    std::stringstream ssin; ssin << word->c_str();
                    F.read(ssin, tmp);
                    F.divin(coeff, tmp);
// std::clog << "### div : " << coeff << std::endl;
                } else {
                    --word;
                }

// M.write(std::clog << "# BEF opR M:", FileFormat::Pretty) << std::endl;
                opRow(M, i, M[j], coeff);
// M.write(F.write(std::clog << "# AFT ", coeff) << " opR M:", FileFormat::Pretty) << std::endl;
            }
        }

#ifdef VERBATIM_PARSING
        std::clog << "# I: " << inputs << std::endl;
        std::clog << "# V: " << variables << std::endl;
        M.write(std::clog << "# M:", FileFormat::Pretty) << std::endl;
        std::clog << std::string(40,'#') << std::endl;
#endif

    }

    const size_t n(inputs.size());
    A.resize(0,n);
    for(const auto& [variable, index] : variables) {
        if(variable[0] == outchar) {
            const size_t i(std::stoi(variable.substr(1,std::string::npos)));
            if (i >= A.rowdim()) A.resize(i+1,n);
            setRow(A,i,M,index);
        }
    }

    _Mat T(F, A.coldim(), A.rowdim()); Transpose(T,A);
    _Mat B(F, 0, A.rowdim());
    for(const auto& [input, index] : inputs) {
        const size_t j(std::stoi(input.substr(1,std::string::npos)));
        if (j >= B.rowdim()) B.resize(j+1,T.coldim());
        setRow(B,j,T,index);
    }

    Transpose(A,B);


// A.write(std::clog, FileFormat::Pretty) << std::endl;

    return A;
}
// ============================================================



// ============================================================
// Recursive extraction of parenthesis into a new variable
size_t extractParenthesis(VProgram_t& newP, std::vector<std::string>& line,
                          const char freechar, size_t& tmpnum) {
#ifdef VERBATIM_PARSING
    std::clog << "# initial    line: " << line  << std::endl;
#endif
    size_t ep(0);
    auto openp = std::find(line.begin(),line.end(),"(");
    if (openp != line.end()) {
        size_t depth(1);
        const auto ocp = { "(", ")" };
        auto closp = openp;
        do {
            closp = std::find_first_of(closp+1,line.end(),ocp.begin(),ocp.end());

            if (*closp == ")") --depth;
            else ++depth;
        } while(depth>0);
        std::vector<std::string> newline;
        std::string newvar;
        newvar += freechar; newvar += std::to_string(++tmpnum);
        newline.push_back(newvar);
        newline.emplace_back(":=");
        newline.insert(newline.end(), openp+1, closp);
        newline.emplace_back(";");
        *openp = newvar;
        line.erase(openp+1, closp+1);
#ifdef VERBATIM_PARSING
        std::clog << "# new created line: " << newline  << std::endl;
#endif
        ep = extractParenthesis(newP, newline, freechar, tmpnum);
        newP.push_back(newline);
#ifdef VERBATIM_PARSING
        std::clog << "# replaced    line: " << line
                  << " (" << (ep+1) << ')' << std::endl;
#endif
        ep += extractParenthesis(newP, line, freechar, tmpnum);
        return ep+1;
    } else
        return 0;
}
// ============================================================






// ============================================================
// Transform program by creating temporaries for parenthesis blocks
VProgram_t& parenthesisExpand(VProgram_t& P) {
        // ==================================
        // Find two unused variable names
    std::set<char> varsChar;
    for(const auto& line: P) for(const auto& word: line)
        varsChar.insert(word[0]);
    char freechar(unusedChar(varsChar));
    size_t tmpnum(0);

    VProgram_t nProgram;
    for(auto& line: P) {
        extractParenthesis(nProgram, line, freechar, tmpnum);
        nProgram.push_back(std::move(line));
    }
    return P=std::move(nProgram);
}
// ============================================================
