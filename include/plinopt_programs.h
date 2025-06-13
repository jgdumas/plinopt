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
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/

#include "plinopt_library.h"

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <givaro/givprint.h>

#ifndef _PLINOPT_LIBRARY_PROGRAMS_H_
#define _PLINOPT_LIBRARY_PROGRAMS_H_

// ============================================================
// Program is vector of lines
typedef std::vector<std::vector<std::string>> VProgram_t;

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

// ============================================================



// ============================================================
// Main parsing procedure, writing the ouput program
int Tellegen(std::istream& input,
             const char ichar = 'i', const char ochar = 'o',
             const char cchar = 'c', const char zchar = 'z');
// ============================================================


// ============================================================
// Creating a vector of lines from a text file
VProgram_t& programParser(VProgram_t& ProgramVector, std::stringstream& ssin);
// ============================================================


// ============================================================
// Rewrites:
//   [1] temporarily replaces output variables by temporary, except last ones
//   [2] replaces variable only assigned to temporary, directly by input
//   [*] Removes all no-op operations like ai:=ai;
//   [3] Backward reassingment of output variables
//   [4] (optional) rewrites singly used variables in-place
size_t variablesTrimer(VProgram_t& P, const bool simplSingle=true,
                       const char inchar = 'i', const char outchar = 'o');
// ============================================================


// ============================================================
// Recreate matrix from Program
template<typename _Mat>
_Mat& matrixBuilder(_Mat& A, const VProgram_t& P, const char outchar = 'o');
// ============================================================


// ============================================================
// Transform program by creating temporaries for parenthesis blocks
VProgram_t& parenthesisExpand(VProgram_t& P);
// ============================================================


// ============================================================
// Testing for a natural number within a string
bool isNatural(const std::string& s);
// ============================================================

// ============================================================
std::ostream& printline(std::ostream& sout,
                        const std::vector<std::string>& line);
// ===============================================================

// ===============================================================
size_t progSize(const VProgram_t& P) ;
Pair<size_t> progOperations(const VProgram_t& P) ;
// ===============================================================

// ===============================================================
std::ostream& operator<<(std::ostream& sout, const VProgram_t& P);
std::ostream& operator<<(std::ostream& out,
                         const std::map<std::string,size_t>& m);
// ===============================================================

// ============================================================
// Stack of lines (strings), to print in reverse order
struct stackbuf;
// ============================================================

// ============================================================
// Find unused variable name
char unusedChar(const std::set<char>& C, const char cstart = 'a'-1);
// ============================================================


#include "plinopt_programs.inl"

#endif
