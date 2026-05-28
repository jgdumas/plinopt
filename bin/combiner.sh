#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Wrapper testing combination of variables
# ==========================================================================

FIL=$1
SMS=/tmp/m_$$.sms
SLP=/usr/local/soft/plinopt/bin/SLPchecker
${SLP} ${FIL} > ${SMS}

VARS=$(sed 's/:=.*//g' ${FIL})
LVARS=(`echo "${VARS}"`)
OVARS=$(sed -r "s/([:+-])/ /g" ${FIL}|awk '{if (NF>3) print $1}')
LOARS=(`echo "${OVARS}"`)
# echo "LVARS ${#LVARS[@]}: ${LVARS[@]}"
# echo "LOARS ${#LOARS[@]}: ${LOARS[@]}"


echo "# [combiner] ${#LOARS[@]} x ${#LVARS[@]} x ${#LVARS[@]}: [${LOARS[@]}] ${LVARS[@]}"

function Show() {
    local NAM=$1
    echo -e "----- BEG: ${NAM} -----\n${!NAM}\n----- END: ${NAM} -----"
}

GRE='\033[1;32m'
RED='\033[0;41m'
BLU='\033[0;36m'
NC='\033[0m'    # No Color

for var in ${OVARS[@]}; do
    OTH=$(egrep -v "(${var})" <<< "${VARS[@]}")
    echo -n "${var}:"
    for rav in ${OTH[@]}; do
	THD=$(egrep -v "(${rav})" <<< "${OTH[@]}")
	for avr in ${THD[@]}; do
	    CMD="${var}:=${rav}+${avr};";
	    # echo "${var} : ${rav} : ${avr} : ${CMD}"
	    TSR=$(egrep -v "${var}:=" ${FIL})
	    TST="${TSR}"`echo -e "\n${CMD}"`
	    # Show TST

	    if ("${SLP}" -M "${SMS}" <<< "${TST}" 2> /dev/null); then
		echo -e "\n${GRE}/!\ IMPROVEMENT /!\:\t${NC} ${CMD} \n"
		"${SLP}" -M "${SMS}" <<< "${TST}"
	    fi
	    CMD="${var}:=${rav}-${avr};";
	    # echo "${var} : ${rav} : ${avr} : ${CMD}"
	    TST="${TSR}"`echo -e "\n${CMD}"`
	    # Show TST

	    if ("${SLP}" -M "${SMS}" <<< "${TST}" 2> /dev/null); then
		echo -e "\n${GRE}/!\ IMPROVEMENT /!\:\t${NC} ${CMD} \n"
		"${SLP}" -M "${SMS}" <<< "${TST}"
	    fi
	done
	stdbuf -o0 printf '+'
    done
    echo ""
done

\rm -rf ${SMS}