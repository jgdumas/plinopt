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

std::ostream& operator<<(std::ostream& out, const std::map<std::string,size_t>& m) {
    out << '{';
    for(const auto& [variable, index] : m)
        out << variable << '(' << index << ')' << ',';
    return out << '}';
}


// template<typename _Row, typename Field>
// inline void negRow(_Row& r, const Field& F) {
//     for(auto e: r) F.negin(e.second);
// }

template<typename _Mat, typename Field>
inline void opRow(_Mat& M, const size_t i, const typename _Mat::Row& s,
                  const typename Field::Element& c, const Field& F) {
    for(auto e: s) {
        const size_t j(e.first);
        typename Field::Element t; F.init(t);
        F.assign(t,M.getEntry(i,j)); // might be zero
        F.axpyin(t,c,e.second);
        M.setEntry(i,j,t);
    }
}


// ============================================================
template<typename _Mat, typename Field>
_Mat& matrixBuilder(_Mat& A, const VProgram_t& P, const Field& F, const char outchar = 'o') {

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
// std::clog << "## word: " << *word << std::endl;
            typename Field::Element coeff; F.init(coeff); F.assign(coeff, F.one);
            if (*word == ";") break;
            if (*word == "-") {
                F.mulin(coeff, F.mOne); ++word;
            } else if (*word == "+") {
                ++word;
            }
            if (*word == output) {
                if (F.isMOne(coeff)) negRow(M, i, F);
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
                if (*word == "*") {
                    ++word;
                    typename Field::Element tmp;
                    std::stringstream ssin; ssin << word->c_str();
                    F.read(ssin, tmp);
                    F.mulin(coeff,tmp);
                } else if (*word == "/") {
                    ++word;
                    typename Field::Element tmp;
                    std::stringstream ssin; ssin << word->c_str();
                    F.read(ssin, tmp);
                    F.divin(coeff, tmp);
                }  else {
                    --word;
                }
                M.setEntry(i,j,coeff);
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
                opRow(M, i, M[j], coeff, F);
//  M.write(F.write(std::clog << "# AFT ", coeff) << " opR M:", FileFormat::Pretty) << std::endl;
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
            setRow(A,i,M,index,F);
        }
    }
    
    _Mat T(F,A.coldim(), A.rowdim()); Transpose(T,A);
    _Mat B(F,0, A.rowdim()); 
    for(const auto& [input, index] : inputs) {
        const size_t j(std::stoi(input.substr(1,std::string::npos)));
        if (j >= B.rowdim()) B.resize(j+1,T.coldim());
        setRow(B,j,T,index,F);
    }
    
    Transpose(A,B);


// A.write(std::clog, FileFormat::Pretty) << std::endl;

    return A;
}




// ============================================================




// ============================================================
template<typename _Mat>
_Mat& PMcheck(_Mat& A, std::istream& input) {
    std::stringstream ssin; ssin << input.rdbuf();

        // Line by line parsing
    VProgram_t ProgramVector; programParser(ProgramVector, ssin);
    const size_t PVs { progSize(ProgramVector) };
    std::clog << std::string(40,'#') << std::endl;

    QRat QQ;


    return matrixBuilder(A, ProgramVector, QQ);
}
// ============================================================




// ============================================================
// Main: select between file / std::cin
int main(int argc, char** argv) {
    std::string prgname, matname;

    for (int i = 1; argc>i; ++i) {
        std::string args(argv[i]);
        if (args == "-h") {
            std::clog << "Usage: " << argv[0] << "[-M file.sms] [stdin|file.prg] \n";
            exit(-1);
        } else if (args == "-M") { matname = std::string(argv[++i]); }
        else { prgname = args; }
    }

    QRat QQ;
    Matrix A(QQ);

    if (prgname == "") {
        PMcheck(A, std::cin);
    } else {
        std::ifstream ifile(prgname);
        if ( ifile ) {
            PMcheck(A, ifile);
            ifile.close();
        }
    }

    if (matname != "") {
        std::ifstream matfile(matname);
        QMstream ms(QQ, matfile);
        Matrix B(ms); const size_t m(B.rowdim()), n(B.coldim());
        Matrix R(QQ,m,n);
        LinBox::MatrixDomain<QRat> BMD(QQ);
        BMD.sub(R,B,A);

        if (BMD.isZero (R))
            std::clog <<"# \033[1;32mSUCCESS: correct "
                      << m << 'x' << n
                      << " Matrix-Vector multiplication!\033[0m" << std::endl;
        else {
            std::cerr << "# \033[1;31m****** ERROR, not a "
                      << m << 'x' << n
                      << " m-v algorithm******\033[0m"
                      << std::endl;

            A.write(std::clog, FileFormat::Pretty) << std::endl;
            B.write(std::clog, FileFormat::Pretty) << std::endl;
        }
    } else
        A.write(std::clog, FileFormat::Pretty) << std::endl;


    return 0;
}
// ============================================================
