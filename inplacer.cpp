// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet
// ==========================================================================

/****************************************************************
 * Returns an in-place program computing a linear function
 * Usage: L.sms R.sms P.sms [expansion]
 *
 * Reference:
 *      [ In-place accumulation of fast multiplication formulae
 *        J-G. Dumas, B. Grenet
 *        https://hal.science/hal-04167499 ]
 *
 * Matrix syntax: SMS format, see:
 *      [Sparse Integer Matrix Collection](https://hpac.imag.fr)
 *		- Starts with: `m n 'R'`
 *		- then: `i j value`
 *		- ends with: `0 0 0`
 ****************************************************************/

#include "plinopt_inplace.h"


void usage(int argc, char ** argv) {
    std::clog << "Usage: " << argv[0] << " L.sms [-O #]\n"
              << "  -O #: randomized search with that many loops\n";
    exit(-1);
}


std::ostream& FindProgram(std::ostream& out, std::istream& input,
                          const size_t randomloops, const bool transposed) {
    QRat QQ; QMstream qs(QQ, input); Matrix M(qs);
    const char inchar(transposed?'t':'i'), tmpchar('z'), outchar('o');
    const size_t outdim(transposed?M.rowdim():M.coldim());
    const size_t innerdim(transposed?M.coldim():M.rowdim());

#ifdef VERBATIM_PARSING
    M.write(std::clog << "M:=Matrix(",FileFormat::Maple) << ");" << std::endl;
#endif
    std::clog << std::string(40,'#') << std::endl;

    Tricounter opcount; AProgram_t Program;
    LinBox::Permutation<QRat> P(QQ, innerdim);

    if (transposed) {
        Matrix T(QQ, M.coldim(), M.rowdim()); Transpose(T, M);
        opcount = SearchLinearAlgorithm(Program, P, T, tmpchar, randomloops);
    } else {
        opcount = SearchLinearAlgorithm(Program, P, M, tmpchar, randomloops);
    }

        // Print the chosen algorithm
    input2Temps(out, outdim, inchar, tmpchar);
    printwithOutput(out, outchar, Program, P) << std::flush;

#ifdef INPLACE_CHECKER
    std::clog << std::string(40,'#') << std::endl;
    std::clog << '<';
    for(size_t i=0; i< outdim; ++i) {
        if (i != 0) std::clog << '|';
        std::clog << tmpchar << i << '-' << inchar << i;
    }
    std::clog << '>' << ';' << std::endl;
#endif
        // =================================
        // Resulting program operation count
    std::clog << std::string(40,'#') << std::endl;
    std::clog << "# \033[1;32m" << std::get<0>(opcount) << "\tADD\033[0m\n";
    std::clog << "# \033[1;32m" << std::get<1>(opcount) << "\tSCA\033[0m\n";
    std::clog << "# \033[1;32m" << std::get<2>(opcount) << "\tROWS\033[0m\n";
    std::clog << std::string(40,'#') << std::endl;

    return out;
}



// ===============================================================
// argv[1]: L.sms
//   -O # randomized search of that many loops
//        looks for reduced number of additions, then multiplications
int main(int argc, char ** argv) {

    size_t randomloops(DORANDOMSEARCH?DEFAULT_RANDOM_LOOPS:1);
    std::vector<std::string> filenames;
    bool transposed(false);

    for (int i = 1; i<argc; ++i) {
        std::string args(argv[i]);
        if (args == "-h") { usage(argc,argv); }
        else if (args == "-t") { transposed=true; }
        else if (args == "-O") {
            randomloops = atoi(argv[++i]);
            if ( (randomloops>1) && (!DORANDOMSEARCH) ) {
                randomloops = 1;
                std::cerr << "#  \033[1;36mWARNING: RANDOM_TIES not defined,"
                          << " random loops disabled\033[0m." << std::endl;
            }
        } else { filenames.push_back(args); }
    }
        // =================================
        // Reading matrices
    if (filenames.size() > 0) {
        for(size_t i(0); i<filenames.size(); ++i) {
            std::ifstream input (filenames[i]);
            FindProgram(std::cout, input, randomloops, transposed) << std::flush;
        }
    } else {
        FindProgram(std::cout, std::cin, randomloops, transposed) << std::flush;
    }

    return 0;
}
