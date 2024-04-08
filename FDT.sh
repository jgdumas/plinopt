#!/bin/bash
tmpfile=/tmp/fdt_plinopt.$$
numopt=10
modulus=7

if [ "$#" -ge 1 ]; then
    fics=("$@")
else
    fics=(data/*sms)
fi

for fic in ${fics[@]}
do
    echo "${fic}:"
    ((./optimizer -O ${numopt} $fic | ./compacter -s | ./PMchecker -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
    ((./matrix-transpose $fic | ./optimizer -q ${modulus} -O ${numopt} | ./transpozer | ./compacter -s | ./PMchecker  -q ${modulus} -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
    ((./matrix-transpose $fic | ./optimizer -O ${numopt} | ./transpozer | ./compacter -s | ./PMchecker -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
done


let tries=3*${#fics[@]}
success=`grep SUCCESS ${tmpfile} | wc -l`

echo -ne "\033[1;93mSUCCESS: "
if [ ${success} -ne ${tries} ]; then
    echo -ne "\033[1;91m"
fi
echo -e "${success} \033[1;93m/ ${tries}\033[0m"
\rm -rf ${tmpfile}
