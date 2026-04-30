#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Counting addition/subtraction and multiplication/division signs in an SLP
#   (not including leading negations/additions, like 'a:=-t')
# ==========================================================================

file=$*

GRE='\033[0;92m'
RED='\033[0;91m'
BLU='\033[0;96m'
NC='\033[0m'    # No Color

operations=`sed 's/\#.*//' ${file} | sed -E 's/:=-/:=/g;s/:=+/:=/g;s/[*/]([^1-9c])/.\1/;s/[*/][0-9]*[*/]/*/g;s/([^-+*/&.]*)([-+*/&.]*)/\2/g' | tr -d '\n'`


btf=`echo -n ${operations} | sed -E 's/[^&]//g'| wc -m`
add=`echo -n ${operations} | sed -E 's/[^+-]//g'| wc -m`
sca=`echo -n ${operations} | sed -E 's/[^*/]//g'| wc -m`
trk=`echo -n ${operations} | sed -E 's/[^\.]//g'| wc -m`

if [[ $((trk)) -gt 0 ]]; then
    echo -e "${BLU}${trk}\trank${NC}"
fi
echo -e "${add}\tadditions/subtractions"
echo -e "${sca}\tscalings (mul/div)"
if [[ $((btf)) -gt 0 ]]; then
    echo -e "${btf}\tbutterflies"
fi

# Total
((add+=sca))
((add+=btf*2))
if [[ $((trk)) -gt 0 ]]; then
	echo -en "${GRE}${trk}M+${add}A="
	((trk+=add))
	echo -n ${trk}
else
	echo -en "${GRE}${add}"
fi
echo -e "\tTotal${NC}"
