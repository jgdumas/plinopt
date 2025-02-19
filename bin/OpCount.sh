#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Counting addition/subtraction and multiplication/division signs in an SLP
#   (not including leading negations, like 'a:=-t')
# ==========================================================================


file=$*

operations=`sed -E 's/:=-/:=/g;s/[*/][0-9]*[*/]/*/g;s/([^-+*/]*)([-+*/]*)/\2/g' ${file} | tr -d '\n'`

echo -e "`echo -n ${operations} | sed -E 's/[*/]//g'| wc -m`\tadditions/subtractions"
echo -e "`echo -n ${operations} | sed -E 's/[+-]//g'| wc -m`\tmultiplications/divisions"
