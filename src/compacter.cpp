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
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/


#include "plinopt_programs.h"

// ============================================================
// Main: select between file / std::cin
int main(int argc, char** argv) {
    bool simplSingle(true);
    std::string filename;
    size_t numloops(0);

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << "[-s/-n] [-O #] [stdin|file.prg]\n"
                      << "  -s/-n: replace/not-replace singly used variables\n"
                      << "  -O #: number of trim loops (default until stable)"
                      << std::endl;
            exit(-1);
        }
        else if (args == "-s") { simplSingle = true; }
        else if ((args == "-n") || (args == "-ns")) { simplSingle = false; }
        else if (args == "-O") { numloops = atoi(argv[++i]); }
        else { filename = args; }
    }

    if (filename == "") {
        PLinOpt::Compacter(std::cout, std::cin, numloops, simplSingle);
    } else {
        std::ifstream ifile(filename);
        if ( ifile ) {
            PLinOpt::Compacter(std::cout, ifile, numloops, simplSingle);
            ifile.close();
        }
    }

    return 0;
}
// ============================================================
