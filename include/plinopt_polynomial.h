// ==========================================================================
// PLinOpt: C++ routines handling linear, bilinear & trilinear programs
// Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
// ==========================================================================

/****************************************************************
 * LinBox::PolynomialRing with intuitive polynomial interface
 ****************************************************************/

#ifndef _PLINOPT_POLYNOMIAL_H_
#define _PLINOPT_POLYNOMIAL_H_

#include <iostream>
#include <givaro/givpoly1.h>

namespace PLinOpt {
// ============================================

#  if 0

// Legacy reading interface is LinBox::PolynomialRing
template<typename Ring> using PRing = Givaro::Poly1Dom<Ring,Givaro::Dense>;

#  else

// PolynomialDomain with reading interface matching the writing one
template<typename Ring> class PRing : public Givaro::Poly1Dom<Ring,Givaro::Dense> {
    typedef Givaro::Poly1Dom<Ring,Givaro::Dense> Parent_t;
public:
    using Parent_t::Element;
    typedef typename Parent_t::Element* Element_ptr;
    typedef const typename Parent_t::Element* ConstElement_ptr;
    using Parent_t::Type_t;
    using Parent_t::Parent_t; // using inherited constructors

    std::istream& read ( std::istream& i, typename Parent_t::Element& P) const {
        static std::ostringstream Xstream; Xstream << this->getIndeter();
        static std::string X(Xstream.str()); // Indeterminate

        this->init(P); // reset P
        std::string s; std::getline(i, s);
        s.erase(std::remove_if(s.begin(), s.end(), ::isspace),s.end());

            // Break s into units at every '+' or '-'; the limits will be [p,q)
            // Units should be of the generic form c*X^n:
            //   --> coefficient before X, degree afterwards
        size_t p = 0, q = p;
        while ( q < s.size() ) {
            for ( q = p + 1; q < s.size() && s[q] != '+' && s[q] != '-'; ++q );
            std::string unit(s.substr( p, q - p ));
// std::clog << "unit: " << unit << std::endl;
                // Identify coefficient (c) and exponent (n)
                // First comes the coefficient c
            typename Ring::Element c;
            this->getDomain().assign(c,this->getDomain().one);
            Givaro::Degree n(0);
            const size_t pos(unit.find( X )); // position of char X
            if ( pos == std::string::npos ) { // X not found; pure number
                std::stringstream( unit ) >> c;
            } else {
                if ( pos != 0 ) { // pos == 0 would mean default c = 1
                    const std::string first = unit.substr( 0, pos );
// std::clog << "first: " << first << std::endl;
                    if ( first != "+" ) { // just "+" would mean default c = +1
                        if ( first == "-" ) {
                                // just "-" means -1
                            this->getDomain().assign(c,this->getDomain().mOne);
                        } else {
                                // get the coefficient
                            std::stringstream( first ) >> c;
                        }
                    }
                }

                    // Second, X, followed by the degree n, thus n >= 1
                n = 1;
                if ( pos != unit.size() - 1 ) { // Last would mean just X, i.e. degree 1
                        // monomial c * X ^ n
                        // in fact X followed by any non-number, then number
                    const std::size_t exp(unit.find_first_of("0123456789",pos));
                    if (exp != std::string::npos) {
                                // get the degree
                        std::stringstream( unit.substr( exp ) ) >> n;
                    }
                }
            }

// std::clog << "# coeff.: " << c << std::endl;
// std::clog << "# degree: " << n << std::endl;

                // Accumulate the found monomial c*X^n, to the current polynomial
            typename Parent_t::Element monomial;
            this->init(monomial, Givaro::Degree(n), c);
            this->addin(P,monomial);
            p = q;
        }

// this->write(std::clog << "# current parent-read " << P << ':', P) << std::endl;
        return i;
    }
};
#  endif


} // End of namespace PLinOpt
// ============================================

#endif
