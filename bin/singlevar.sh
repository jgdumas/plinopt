#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Replaces all temporary variables names by variabbles with the same char x
# ==========================================================================

FIL=
CHR=

while [[ $# -gt 0 ]]; do
  case $1 in
    -c)
    CHR="$2"
    shift # past argument
    shift # past value
    ;;
    -h)
      echo "Usage: $0 [-c x] f.slp"
      exit 1
      ;;
    *)
    FIL=$1
    shift # past argument
    ;;
  esac
done

SLP=
if [[ "$FIL" == "" ]]; then
  SLP=$(cat)
else
  SLP=`cat ${FIL}`
fi

RDS=`echo "${SLP}" | egrep -v '(^o|^A|^B|^C)'`

if [[ "$CHR" == "" ]]; then
  VARS=(`echo "${RDS}" | sed -r 's/(.).*/\1/' | sort -u | tr '\n' ' '`)
#   echo "VARS: ${VARS[@]}"

  CHR='d'
  while grep -q ${CHR} <<< ${VARS[@]}; do
    CHR=`echo ${CHR} | tr "a-z" "b-za"`
  done

#   echo "CHR: ${CHR}"
fi

VARS+=(${CHR})
ZHR='d'
while grep -q ${ZHR} <<< ${VARS[@]}; do
    ZHR=`echo ${ZHR} | tr "a-z" "b-za"`
done

#   echo "ZHR: ${ZHR}"

SED=`echo "${RDS}" | sed -r "s/${CHR}/${ZHR}/g;s/:=.*//" | sort -u | tac | awk -v chr=${CHR} 'BEGIN{s=10;ORS=""} {print "s/"$1"/"chr""s"/g;";s++}'`

# echo "SED: ${SED}"

echo "${SLP}" | sed "s/${CHR}/${ZHR}/g;${SED}"
