#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Tests: optimizer, compacter, transpozer, factorizer, sparsifier, inplacer
# ==========================================================================

if [ "$#" -ge 1 ]; then
    fics=("$@")
else
    fics=(data/*slp)
    # Do not test polynomial files yet
    fics=(${fics[@]//*-X_*})
    # Do not test large files automatically
    fics=(${fics[@]//*32x32x32*})
fi

BINDRS=(./ ./bin)

for dir in ${BINDRS[@]}
do
  if [ -x ${dir}/optimizer ]; then
      echo "Found executables in ${dir}"
      BINDIR=${dir}
  fi
done

SLPCK=${BINDIR}/SLPchecker
OPTCT=${BINDIR}/OpCount.sh
CMPCT=${BINDIR}/compacter
TELLG=${BINDIR}/transpozer


let tries=${#fics[@]}
let current=0
let success=0

for fic in ${fics[@]}
do
if [ -e ${fic} ]; then
    echo -n "${fic}: "

    BEF=(`(${SLPCK} ${fic} |& egrep '(additions|multiplications)' | sed 's/\x1b\[[0-9;]*[a-zA-Z]//g'| awk '{print $2}') 2> /dev/null`)
# echo ${BEF[*]}
    AFT=(`${OPTCT} ${fic} | sed 's/\x1b\[[0-9;]*[a-zA-Z]//g'| awk '{print $1}'`)
# echo ${AFT[*]}

    let current=1+${current}
    DADD=$((BEF[0]-AFT[0]))
    DMUL=$((BEF[1]-AFT[1]))
# echo ${DADD} ${DMUL}

    if [ ${DADD} -eq 0 ] && [ ${DMUL} -eq 0 ]; then
	let success=1+${success}
    fi
    echo -ne "\033[1;93mSUCCESS: "
    if [ ${success} -ne ${current} ]; then
	echo -ne "\033[1;91m <${DADD}&${DMUL}> "
    fi
    echo -e "${success} \033[1;93m/ ${tries}\033[0m"
fi
done

if [ ${success} -ne ${current} ]; then
    echo -e "\033[1;91m*** ERROR *** $((current -${success}))/${tries}\033[0m"
else
    echo -e "\033[1;92mSUCCESS: ${success}/${tries}\033[0m"
fi