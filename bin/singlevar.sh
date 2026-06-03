#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Replaces all temporary variables names by variables with the same char x
# ==========================================================================

FIL=
CHR=
STR=10
INP=0

while [[ $# -gt 0 ]]; do
  case $1 in
    -c|-z|--char)
    CHR="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--start)
    STR="$2"
    shift # past argument
    shift # past value
    ;;
    -i|--inplace)
    INP=1
    shift # past argument
    ;;
    -h)
      echo "Usage: $0 [-c|-z x] [-s #] [-i] f.slp"
      exit 1
      ;;
    *)
    FIL=$1
    shift # past argument
    ;;
  esac
done

function Show() {
    local NAM=$1
    >&2 echo -e "----- BEG: ${NAM} -----\n${!NAM}\n----- END: ${NAM} -----"
}

SLP=
if [[ "$FIL" == "" ]]; then
  INP=0
  SLP=$(cat)
else
  SLP=`cat ${FIL}`
fi

RDS=`echo "${SLP}" | egrep -v '(^o|^A|^B|^C)'`
VARS=(`echo "${RDS}" | sed -r 's/(.).*/\1/' | sort -u | tr '\n' ' '`)
# echo "VARS: ${VARS[@]}"

if [[ "$CHR" == "" ]]; then

  CHR='d'
  while grep -q ${CHR} <<< ${VARS[@]}; do
    CHR=`echo ${CHR} | tr "a-z" "b-za"`
  done

fi

# echo "CHR: ${CHR}"
>&2 echo -n "# [SingleVars] variables ${VARS[@]}"
VARS+=(${CHR})
ZHR='d'
while grep -q ${ZHR} <<< ${VARS[@]}; do
    ZHR=`echo ${ZHR} | tr "a-z" "b-za"`
done

# Show RDS
NBD=`echo "${RDS}" | wc -l`
END=$((STR+NBD-1))

SED=`echo "${RDS}" | sed -r "s/${CHR}/${ZHR}/g;s/:=.*//" | sort -u | tac | awk -v chr=${CHR} -v str=${END} 'BEGIN{s=str;ORS=""} {print "s/"$1"/"chr""s"/g;";s--}'`

# Show SED
>&2 echo " (${NBV})  -->  ${CHR}${STR} .. ${CHR}${END}"

# echo "ZHR: ${ZHR}"
# echo "SED: ${SED}"
# echo "CHG: s/${CHR}/${ZHR}/g;${SED}"

if [[ ${INP} -eq 1 ]]; then
  sed -i "s/${CHR}/${ZHR}/g;${SED}" ${FIL}
else
  echo "${SLP}" | sed "s/${CHR}/${ZHR}/g;${SED}"
fi
