#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Rotater: build the rotated Matrix Multiplication Algorithm
#     left = R ; P^T_s ; L_s^T
#     right= P^T_s ; L ; R_s^T
# ==========================================================================

DIR=`dirname $0`
SUFF=left
FILS=()
MFLG="" # "-b 3"

while [[ $# -gt 0 ]]; do
  case $1 in
    -d|--direction)
    SUFF="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--left)
    SUFF=left
    shift # past argument
    ;;
    -r|--right)
    SUFF=right
    shift # past argument
    ;;
    -h|--h|-help|--help|-*|--*)
    echo "Usage: $0 [-d left/right] [-l|-r] L.sms R.sms P.sms"
    echo "  -l|-d left: rotation to the left (default)"
    echo "  -r|-d right: rotation to the right"
    exit 1
    ;;
    *)
    FILS+=($1)
    shift # past argument
    ;;
    esac
done


TR=${DIR}/matrix-transpose
CS=${DIR}/columns-swap

TS=${DIR}/transpozer
CI=${DIR}/input-swap.sh

MMMCHK=${DIR}/MMchecker
SLPCHK=${DIR}/SLPchecker

L=${FILS[0]}
R=${FILS[1]}
P=${FILS[2]}

echo "${MMMCHK} ${L} ${R} ${P}: "
MMC=$(${MMMCHK} ${L} ${R} ${P} 2>&1)
echo -e "${MMC}"

MMD=$(sed 's/x/ /g' <<< ${MMC} | cut -d' ' -f4,5,6)
m=$(cut -d' ' -f1 <<< ${MMD})
k=$(cut -d' ' -f2 <<< ${MMD})
n=$(cut -d' ' -f3 <<< ${MMD})

Lr=`basename ${L}| sed -e "s/_[LRP]\./_${SUFF}&/"`
Rr=`basename ${R}| sed -e "s/_[LRP]\./_${SUFF}&/"`
Pr=`basename ${P}| sed -e "s/_[LRP]\./_${SUFF}&/"`

GROT=0
if [ -f ${Lr} ] || [ -f ${Rr} ] || [ -f ${Pr} ]; then
    read -p "# Overwrite ${Lr}, ${Rr} and ${Pr}? [y/N] " RESP
    if [[ "${RESP}" == "y" ]]; then
    GROT=1
    fi
else
    GROT=1
fi

Ls=$(echo ${L} | sed 's/sms/slp/')
Rs=$(echo ${R} | sed 's/sms/slp/')
Ps=$(echo ${P} | sed 's/sms/slp/')

SLPT=0
if [ -f ${Ls} ] && [ -f ${Rs} ] && [ -f ${Ps} ]; then
    Lp=`basename ${Ls}| sed -e "s/_[LRP]\./_${SUFF}&/"`
    Rp=`basename ${Rs}| sed -e "s/_[LRP]\./_${SUFF}&/"`
    Pp=`basename ${Ps}| sed -e "s/_[LRP]\./_${SUFF}&/"`
    if [ -f ${Lp} ] || [ -f ${Rp} ] || [ -f ${Pp} ]; then
        read -p "# Overwrite ${Lp}, ${Rp} and ${Pp}? [y/N] " RESP
        if [[ "${RESP}" == "y" ]]; then
            SLPT=1
        fi
    else
        SLPT=1
    fi
fi

if [[ "$GROT" -eq 1 ]]; then
    # Do generate the rotation
    if [[ "${SUFF}" == "right" ]]; then
        ${TR} ${P} | ${CS} -m ${n} > ${Lr}
        cat ${L} > ${Rr}
        ${CS} -m ${n} ${R} | ${TR} > ${Pr}
        if [[ "$SLPT" -eq 1 ]]; then
            ${TS} ${Ps} | ${CI} -m ${n} > ${Lp}
            cat ${Ls} > ${Rp}
            ${CI} -m ${n} ${Rs} | ${TS} > ${Pp}
       fi
    else
        cat ${R} > ${Lr}
        ${TR} ${P} | ${CS} -m ${n} > ${Rr}
        ${CS} -m ${k} ${L} | ${TR} > ${Pr}
        if [[ "$SLPT" -eq 1 ]]; then
            cat ${Rs} > ${Lp}
            ${TS} ${Ps} | ${CI} -m ${n} > ${Rp}
            ${CI} -m ${k} ${Ls} | ${TS} > ${Pp}
        fi
    fi

    # Check the results
    if [[ "$SLPT" -eq 1 ]]; then
        ${SLPCHK} -M ${Lr} ${Lp}
        ${SLPCHK} -M ${Rr} ${Rp}
        ${SLPCHK} -M ${Pr} ${Pp}
    fi
    ${MMMCHK} ${MFLG} ${Lr} ${Rr} ${Pr}
fi
