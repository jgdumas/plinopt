// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
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


#include "plinopt_programs.h"




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
    variablesTrimer(ProgramVector, simplSingle);
    sout << ProgramVector;
    const size_t PRs { progSize(ProgramVector) };
    std::clog << std::string(40,'#') << std::endl;

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
