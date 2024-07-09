#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# Tests: trilinear in-placer
# Requires: maple from maplesoft
# ==========================================================================

MAPLEPROG=maple
tmpfile=/tmp/fdt_plinopt.$$
numopt=100


if [ "$#" -ge 1 ]; then
    fics=("$@")
else
    fics=(data/L*sms)
fi

BINDRS=(./ ./bin)

for dir in ${BINDRS[@]}
do
  if [ -x ${dir}/trilplacer ]; then
      echo "Found executables in ${dir}"
      BINDIR=${dir}
  fi
done

TRIPL=${BINDIR}/trilplacer

for fic in ${fics[@]}
do
  if [ -e ${fic} ]; then
      threefics=`./bin/3.sh ${fic}`
      mdims=(`fgrep -v '#' ${threefics} | egrep '( M| R)' | cut -d':' -f2 | tr '\n' ' '`)
      Lm=${mdims[0]}; let Ldm=1+${Lm}
      Ln=${mdims[1]}; let Ldn=2*${Ln}
      Rm=${mdims[3]}; let Rdm=1+${Rm}
      Rn=${mdims[4]}; let Rdn=2*${Rn}
      Pm=${mdims[6]}; let Pdm=2*${Pm}
      Pn=${mdims[7]}; let Pdn=1+${Pn}
      echo -n "${fic} ${Lm}x${Ln} ${Rm}x${Rn} ${Pm}x${Pn} :"
      (${TRIPL} -O ${numopt} `./bin/3.sh $fic`  >& /dev/stdout) | ${MAPLEPROG} | egrep '(errors :=|ADD|SCA|[0-9][^#]*AXPY)' | tr '\n' '\t'| tee -a ${tmpfile}
      echo | tee -a ${tmpfile}
      echo -n "${fic} ${Ldm}x${Ldn} ${Rdm}x${Rdn} ${Pdm}x${Pdn} :"
      (${TRIPL} -e -O ${numopt} `./bin/3.sh $fic`  >& /dev/stdout) | ${MAPLEPROG} | egrep '(errors :=|ADD|SCA|[0-9][^#]*AXPY)' | tr '\n' '\t'| tee -a ${tmpfile}
      echo | tee -a ${tmpfile}
  fi
done


let tries=2*${#fics[@]}
success=`grep 'errors := \[0, 0, 0\]' ${tmpfile} | wc -l`

echo -ne "\033[1;93mSUCCESS: "
if [ ${success} -ne ${tries} ]; then
    echo -ne "\033[1;91m"
fi
echo -e "${success} \033[1;93m/ ${tries}\033[0m"
\rm -rf ${tmpfile}
