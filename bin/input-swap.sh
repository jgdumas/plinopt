#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Swaping the input variables (all i#) of an SLP:
#   considering the vector as an m x n vectorization
#   produces the vectorization of the transposed n x m
# ==========================================================================


IM=0
while [[ $# -gt 0 ]]; do
  case $1 in
    -m|--dim)
    IM="$2"
    shift # past argument
    shift # past value
    ;;
    *)
    IFIL=$1
    shift # past argument
    ;;
    esac
done

FIL=$(cat ${IFIL})

VARS=$(sed 's/:=/ /;s/+/ /g;s/-/ /g;s/;.*/ /;s/\*[0-9]* / /g;s/\/[0-9]* / /g;s/)//g;s/(//g' <<< ${FIL} | tr ' ' '\n')
VINPS=($(grep i <<< ${VARS}| sort -r -u| tr '\n' ' '))
CHARS=($(sed 's/[0-9]//g;/^$/d' <<< ${VARS}| sort -u| tr '\n' ' '))
# >&2 echo "# VINPS: ${VINPS[@]}"
# >&2 echo "# CHARS: ${CHARS[@]}"

if [[ ${IM} -eq 0 ]]; then
    m=$(echo "sqrt(${#VINPS[@]})" | bc)
else
    m=${IM}
fi

n=$(( ${#VINPS[@]} / ${m} ))

>&2 echo "# SLP input columns swap: ${m}x${n} --> (${n}x${m})^T"


NCHAR="a"
while grep -q ${NCHAR} <<< ${CHARS[@]}; do
    NCHAR=`echo ${NCHAR} | tr "a-z" "b-za"`
done
# >&2 echo "# NCHAR: ${NCHAR}"

sedf=""
for var in ${VINPS[@]}; do
    k=$(sed 's/i//' <<< ${var})
    i=$(( ${k} % ${m} ))
    j=$(( (${k}-${i})/${m} ))
    t=$(( ${i}*${n}+${j} ))
    # >&2 echo "# k: ${k}, i: ${i}, j: ${j}, t: ${t}"
    sedf="${sedf}s/${var}/${NCHAR}${t}/g;"
done
sedf="${sedf}s/${NCHAR}/i/g"

OUTP=$(sed "${sedf}" <<< ${FIL})

echo "${OUTP}"
