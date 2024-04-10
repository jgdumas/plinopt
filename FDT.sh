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

for fic in ${fics[@]}
do
    mdims=(`fgrep -v '#' $fic | head -1`)
    m=${mdims[0]}
    n=${mdims[1]}
    echo "${fic} ${m}x${n}:"

    ((./optimizer -O ${numopt} $fic | ./compacter -s | ./PMchecker -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
    ((./matrix-transpose $fic | ./optimizer -q ${modulus} -O ${numopt} | ./transpozer | ./compacter -s | ./PMchecker  -q ${modulus} -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
    ((./matrix-transpose $fic | ./optimizer -O ${numopt} | ./transpozer | ./compacter -s | ./PMchecker -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
    ((./sparsifier -c ${coeffs} $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
    if [ "$m" -ge "$n" ]; then
	((./factorizer -q ${modulus} $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
	((./factorizer -O ${numopt} $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
    else
	((./matrix-transpose $fic |./factorizer -q ${modulus}) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
	((./matrix-transpose $fic |./factorizer -O ${numopt}) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
    fi
done


let tries=6*${#fics[@]}
success=`grep SUCCESS ${tmpfile} | wc -l`

echo -ne "\033[1;93mSUCCESS: "
if [ ${success} -ne ${tries} ]; then
    echo -ne "\033[1;91m"
fi
echo -e "${success} \033[1;93m/ ${tries}\033[0m"
\rm -rf ${tmpfile}
