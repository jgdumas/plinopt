#!/bin/bash
tmpfile=/tmp/fdt_plinopt.$$

if [ "$#" -ge 1 ]; then
    fics=("$@")
else
    fics=("data/*sms")
fi

for fic in ${fics[@]}
do
    echo "${fic}:"
    ((./optimizer -O 10 $fic | ./compacter | ./PMchecker -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
#     ((./matrix-transpose $fic | ./optimizer -O 10 | ./transpozer | ./PMchecker -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
    ((./matrix-transpose $fic | ./optimizer -O 10 | ./transpozer | ./compacter | ./PMchecker -M $fic) >& /dev/stdout) | egrep '(SUCCESS|ERROR)' | tee -a ${tmpfile}
done

tries=`wc -l ${tmpfile}|cut -d' ' -f1`
success=`grep SUCCESS ${tmpfile} | wc -l`

echo -ne "\033[1;93mSUCCESS: "
if [ ${success} -ne ${tries} ]; then
    echo -ne "\033[1;91m"
fi
echo -e "${success} \033[1;93m/ ${tries}\033[0m"
\rm -rf ${tmpfile}
