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

TR=${DIR}/matrix-transpose
CS=${DIR}/columns-swap

L=$1
R=$2
P=$3

SUFF=left

if [ "$#" -ge 4 ]; then
    SUFF=$4
fi


MMchecker $L $R $P

Lr=`basename $L| sed -e "s/_[LRP]\./_${SUFF}&/"`
Rr=`basename $R| sed -e "s/_[LRP]\./_${SUFF}&/"`
Pr=`basename $P| sed -e "s/_[LRP]\./_${SUFF}&/"`

GROT=0
if [ -f ${Lr} ] || [ -f ${Rr} ] || [ -f ${Pr} ]; then
    read -p "# Overwrite ${Lr}, ${Rr} and ${Pr}? [y/N] " RESP
    if [[ "${RESP}" == "y" ]]; then
	GROT=1
    fi
else
    GROT=1
fi

if [[ "$GROT" -eq 1 ]]; then
    # Do generate the rotation
    if [[ "${SUFF}" == "right" ]]; then
	${TR} $P | ${CS} > ${Lr}
	cat $L > ${Rr}
	${CS} $R | ${TR} > ${Pr}
    else
	cat $R > ${Lr}
	${TR} $P | ${CS} > ${Rr}
	${CS} $L | ${TR} > ${Pr}
    fi

    MMchecker ${Lr} ${Rr} ${Pr}
fi
