// ==========================================================================
// PLinOpt: a collection of C++ routines handling linear & bilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * Transposition of programs via Tellegen's principle
 * Program syntax: see main Tellegen function below
 * Reference: [ Tellegen's principle into practice.
 *              ISSAC'03:37-44, A. Bostan, G. Lecerf, É. Schost
 *              https://doi.org/10.1145/860854.860870 ]
 ****************************************************************/

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <numeric>
#include <streambuf>

// ============================================================
// Define to print comments during parsing
//#define VERBATIM_PARSING
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
// Testing for a natural number within a string
bool isNatural(const std::string& s) {
    for(const auto& c: s) if (!std::isdigit(c)) return false;
    return true;
}
// ============================================================


// ============================================================
// Main parsing procedure, writing the ouput program
int Tellegen(std::istream& input) {
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
    std::set<std::string> varSet;
    std::set<std::string> outSet;
    std::set<std::string> inSet;

    std::pair<size_t, size_t> Nops(0,0);

    std::clog << std::string(20,'#') << std::endl;

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
            if (wordVector.front()[0] == 'c') {
#ifdef VERBATIM_PARSING
                std::clog << "# Constant found, unmodified : "
                          << line << std::endl;
#endif
                std::cout << line << std::endl;
                continue;
            }
                // Outputs are transposed into inputs
            if (wordVector.front()[0] == 'o') {
                std::string oneout(wordVector.front());
                inSet.insert( oneout );
                oneout[0] = 'i';
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

                if (variable[0] == 'i') {
                    std::string onein(fixvar);
                    onein[0] = 'o';
                    onein += ":="; onein += fixvar; onein += ';';
                    outSet.insert(onein);
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

    std::clog << std::string(20,'#') << std::endl;

        // ============================================================
        // Select only temporary variables (others are initialized)
    std::vector<std::string> varVector;
    for(const auto& iter: varSet) {
        if ((outSet.find(iter) == outSet.end()) &&
            (inSet.find(iter) == inSet.end())) {
            varVector.push_back(iter);
        }
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

            // ==================================
            // Instead of "xi:=0; xi:=xi + ... ;"
            //          directly write "xi:=...;"
        const auto varlen = line.find(":=",0);
        const std::string variable(line.substr(0,varlen));
        const auto varposvec =
            std::find(varVector.begin(), varVector.end(), variable);

        if ( varposvec != varVector.end()) {
#ifdef VERBATIM_PARSING
            std::clog << "# First temporary, accumulation with " << variable << "=0 is simplified in: " << line << std::endl;
#endif
            const auto accupos = line.find(variable, varlen);

#if DEBUG
            if (accupos != varlen+2) {
                std::cerr << "Problem, first occurence of temporary is not an accumulation." << std::endl;
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
            (pos = line.find_first_of("+-*/;", prev)) !=
                std::string::npos; prev = pos+1) {
            if ((line[prev] == 'o') &&
                (modSet.find(line.substr(prev, pos-prev)) == modSet.end())){
#ifdef VERBATIM_PARSING
                std::clog << "# Unmodified input usage of: " << line.substr(prev, pos-prev) << ", is simplified in RHS: " << line << std::endl;
#endif
                line[prev]='i';
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
#endif

    std::clog << std::string(20,'#') << std::endl;
        // Produce output results "oi:=xi;"
    for(const auto& iter: outSet) std::cout << iter << std::endl;

    const int dimOffset(inSet.size()-outSet.size());


    std::clog << std::string(20,'#') << std::endl;
    std::clog << "# " << Nops.first << "\tadditions ("
              << (Nops.first-dimOffset)
              << (dimOffset<0?'-':'+') << abs(dimOffset)
              << ')' << std::endl;
    std::clog << "# " << Nops.second << "\tmultiplications" << std::endl;
    std::clog << std::string(20,'#') << std::endl;

    return 0;
}
// ============================================================



// ============================================================
// Main: select between file / std::cin
int main(int argc, char** argv) {
    if ( argc > 1 ) {
        std::string args(argv[1]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0] << "[stdin|matrixfile.sms]\n";
            exit(-1);
        }

        std::ifstream ifile(argv[1]);
        if ( ifile ) {
            int rt=Tellegen(ifile);
            ifile.close();
            return rt;
        }
    } else {
        return Tellegen(std::cin);
    }
    return -1;
}
// ============================================================
