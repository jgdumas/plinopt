#!/bin/bash
tmpfile=/tmp/fdt_plinopt.$$
numopt=10
modulus=7
coeffs=5

if [ "$#" -ge 1 ]; then
    fics=("$@")
else
    fics=(data/*sms)
fi

BINDRS=(./ ./bin)

for dir in ${BINDRS[@]}
do
  if [ -x ${dir}/optimizer ]; then
      echo "Found executables in ${dir}"
      BINDIR=${dir}
  fi
done

MTRSP=${BINDIR}/matrix-transpose
OPTIM=${BINDIR}/optimizer
CMPCT=${BINDIR}/compacter
PMCHK=${BINDIR}/PMchecker
TELLG=${BINDIR}/transpozer
FCTZR=${BINDIR}/factorizer
INPLR=${BINDIR}/inplacer
SPSFR=${BINDIR}/sparsifier


for fic in ${fics[@]}
do
    if [ -e ${fic} ]; then
    mdims=(`fgrep -v '#' $fic | head -1`)
    m=${mdims[0]}
    n=${mdims[1]}
    echo "${fic} ${m}x${n}:"

    ((${OPTIM} -O ${numopt} $fic | ${CMPCT} -s | ${PMCHK} -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}

    ((${MTRSP} $fic | ${OPTIM} -q ${modulus} -O ${numopt} | ${TELLG} | ${CMPCT} -s | ${PMCHK}  -q ${modulus} -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}

    ((${MTRSP} $fic | ${OPTIM} -O ${numopt} | ${TELLG} | ${CMPCT} -s | ${PMCHK} -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}

    ((${SPSFR} -c ${coeffs} $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}

    if [ "$m" -ge "$n" ]; then
	((${FCTZR} -q ${modulus} $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}

	((${FCTZR} -O ${numopt} $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
    else
	((${MTRSP} $fic |${FCTZR} -q ${modulus}) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}

	((${MTRSP} $fic |${FCTZR} -O ${numopt}) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
    fi

    ((${INPLR} -O ${numopt} $fic | ${PMCHK} -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
fi
done


let tries=7*${#fics[@]}
success=`grep SUCCESS ${tmpfile} | wc -l`

echo -ne "\033[1;93mSUCCESS: "
if [ ${success} -ne ${tries} ]; then
    echo -ne "\033[1;91m"
fi
echo -e "${success} \033[1;93m/ ${tries}\033[0m"
\rm -rf ${tmpfile}
