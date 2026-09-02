// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * Transposition of programs via Tellegen's principle
 *
 * Syntax:
 *    Files contains a program with the following SYNTAX
 *    [+] constants name start with 'c'
 *    [+] input variables start with 'i' --> transformed into 'o'
 *             supposed to be just a simple affectation "x0:=i0;"
 *    [+] output variables start with 'o' --> transformed into 'i'
 *    [+] lines are of the forms:
 *       xi := sum (yi op(li)), with sum: + or -, and: op * or / or empty
 *           xi, yi are supposed to be without any of "+,-,:=,/,*,;"
 *           liin particular rational numbers a/b should be precomputed,
 *           e.g., into a constant 'cxx:=a/b;' when used in computation.
 *           e.g. replace: "x0:=x1*3/4;" by "c0:=3/4;x0:=x1*c0;"
 *       each yi occurs singly in a given line
 *       One operation per line (';' stops the parsing of that line)
 *
 *
 * References:
 *   [ Tellegen's principle into practice.
 *     ISSAC'03:37-44, A. Bostan, G. Lecerf, É. Schost
 *     https://doi.org/10.1145/860854.860870 ]
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/

// ============================================================
// Define to print comments during parsing
//#define VERBATIM_PARSING
// ============================================================

#include "plinopt_programs.h"

// ============================================================
int Transpozer(std::istream& input, bool compactOutput=true,
               size_t numloops=0, bool simplSingle=true) {
    std::stringstream ssin; ssin << input.rdbuf();
    PLinOpt::VProgram_t progV; PLinOpt::programParser(progV, ssin);
    char next('d'); PLinOpt::parenthesisExpand(progV, next);
    std::clog << std::string(40,'#') << std::endl;

    const size_t ps(PLinOpt::progSize(progV));
    std::clog << "# Transpozing " << ps << " elements" << std::endl;

    if (ps == 0) return -1;

    std::stringstream ssout;
    for(const auto& line: progV) {
        for(const auto& word: line) {
            ssout << word;
        }
        ssout << std::endl;
    }

#ifdef VERBATIM_PARSING
#  if VERBATIM_PARSING >= 4
    std::clog << progV << std::endl;
    std::clog << std::string(40,'#') << std::endl;
    std::clog << "# Expanded program:\n" << ssout.str() << std::endl;
    std::clog << std::string(40,'#') << std::endl;
#  endif
#endif

    std::stringstream rawslp;
    int rt=PLinOpt::Tellegen(rawslp,ssout);

    if (compactOutput) {
        PLinOpt::Compacter(std::cout, rawslp, numloops, simplSingle);
    } else {
        std::cout << rawslp.str() << std::flush;
    }
    return rt;
//     return PLinOpt::Tellegen(std::cout,ssout);
}




// ============================================================
// Main: select between file / std::cin for Transpozition
int main(int argc, char** argv) {
    bool compactOutput(true);
    size_t numloops(0);
    bool simplSingle(true);
    std::string filename;

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << "[-r|-c [-O #] [-s/-n]] [stdin|file.prg]\n"
                      << "  -r/-c: raw/compacted input (default compact)\n"
                      << "  if compact:\n"
                      << "    -s/-n: replace or not singly used variables\n"
                      << "    -O #: number of trim loops (default until stable)"
                      << std::endl;
            exit(-1);
        }
        else if (args == "-r") { compactOutput = false; }
        else if (args == "-c") { compactOutput = true; }
        else if (args == "-s") { simplSingle = true; }
        else if ((args == "-n") || (args == "-ns")) { simplSingle = false; }
        else if (args == "-O") { numloops = atoi(argv[++i]); }
        else { filename = args; }
    }

    if (filename == "") {
        return Transpozer(std::cin, compactOutput, numloops, simplSingle);
    } else {
        std::ifstream ifile(filename);
        if ( ifile ) {
            int rt=Transpozer(ifile, compactOutput, numloops, simplSingle);
            ifile.close();
            return rt;
        }
    }
    return -1;
}
// ============================================================
