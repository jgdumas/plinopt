// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * Compacting straight-line programs
 * Program syntax: see Compacter function below
 *                 see also: optimizer.cpp and transpozer.cpp
 * - Removes no-op
 * - Replaces variable only assigned to temporary, directly by input
 * - Rewrites singly used variables in-place
 * - Reduces leading minus usage '-'
 * Reference:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic; Feb. 2024
 *     Strassen's algorithm is not optimally accurate
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/


#include "plinopt_programs.h"




// ============================================================
// Main compacting procedure, parsing then rewriting
std::ostream& Compacter(std::ostream& sout, std::istream& input,
                        const size_t numloops,
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
    size_t prevPRs(0), currPRs(PVs);

    int iter(numloops); // decreasing 0 will nver be == 0
    do {
        prevPRs = currPRs;
        variablesTrimer(ProgramVector, simplSingle);
        currPRs = progSize(ProgramVector);
        std::clog << "# \033[1;32m" << currPRs << "\telements\tinstead of "
                  << prevPRs << "\033[0m" << std::endl;
#ifdef VERBATIM_PARSING
        std::clog << ProgramVector;
        std::clog << std::string(40,'#') << std::endl;
#endif
    } while ( (currPRs < prevPRs) && (--iter != 0) ) ;

    sout << ProgramVector;

       // Comparing number of elements in the programs
    std::clog << "# \033[1;32m" << currPRs << "\telements\tinstead of "
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
    size_t numloops(0);

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << "[-s] [-O #] [stdin|file.prg]\n"
                      << "  -s: also replace singly used variables\n"
                      << "  -O #: number of trim loops (default until stable)"
                      << std::endl;
            exit(-1);
        }
        else if (args == "-s") { simplSingle = true; }
        else if (args == "-O") { numloops = atoi(argv[++i]); }
        else { filename = args; }
    }

    if (filename == "") {
        Compacter(std::cout, std::cin, numloops, simplSingle);
    } else {
        std::ifstream ifile(filename);
        if ( ifile ) {
            Compacter(std::cout, ifile, numloops, simplSingle);
            ifile.close();
        }
    }

    return 0;
}
// ============================================================
