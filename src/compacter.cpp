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
    bool inPlace(false);
    std::string filename;
    size_t numloops(0);

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0]
                      << "[-s/-n] [-O #] [stdin|file.slp]\n"
                      << "  -i: overwrites file.slp with output\n"
                      << "  -s/-n: replace/not-replace singly used variables\n"
                      << "  -O #: number of trim loops (default until stable)"
                      << std::endl;
            exit(-1);
        }
        else if (args == "-s") { simplSingle = true; }
        else if (args == "-i") { inPlace = true; }
        else if ((args == "-n") || (args == "-ns")) { simplSingle = false; }
        else if (args == "-O") { numloops = atoi(argv[++i]); }
        else { filename = args; }
    }

    if (filename == "") {
        PLinOpt::Compacter(std::cout, std::cin, numloops, simplSingle);
    } else {
        std::fstream file(filename);
        if ( file ) {
            if (inPlace) {
                std::ostringstream sout;
                PLinOpt::Compacter(sout, file, numloops, simplSingle);
                file.close();
                file.open(filename, std::ofstream::out | std::ofstream::trunc);
                file << sout.str();
            } else {
                PLinOpt::Compacter(std::cout, file, numloops, simplSingle);
            }
            file.close();
        }
    }

    return 0;
}
// ============================================================
