#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Wrapper selecting variables to create subprograms
#         calls mirabelle.sh with each found variable
# ==========================================================================


DIR=`dirname $0`
OPTFLAGS=""
ALLFLAGS=""
VARS=()
MOD=""
LIM=""
ALL=0

while [[ $# -gt 0 ]]; do
  case $1 in
    -O|--Optflags)
    OPTFLAGS="-O $2"
    shift # past argument
    shift # past value
    ;;
    -r|-m|--modular)
    MOD="-m $2"
    shift # past argument
    shift # past value
    ;;
    -l|--limits)
    LIM="$2"
    shift # past argument
    shift # past value
    ;;
    -c)
    LIM="o|c_"
    shift # past argument
    ;;
    -n|-nl|--no-limit)
    LIM=""
    ALL=0
    shift # past argument
    ;;
    -a|--all|--all-vars)
    ALL=1
    ALLFLAGS="-O $2"
    shift # past argument
    shift # past value
    ;;
    -h|--h|-help|--help|-*|--*)
    echo "Usage: $0 [-O|-a #] [-q|-r #] [-l pattern] [-n|-nl] P.slp"
      exit 1
      ;;
    *)
    FIL=$1
    shift # past argument
    ;;
    esac
done


VARS=(`egrep -v '(^o)' ${FIL} | egrep "(${LIM})" | cut -d':' -f1|sort -u`)
VARV="`echo ${VARS[@]}| sed 's/ / -v /g'`"
echo "VARS: ${VARS[@]}"
if [[ ${#VARS[@]} -gt 0 ]]; then
  LOG=prune_$$.log
  for var in ${VARS[@]}; do
    ${DIR}/mirabelle.sh -v $var ${OPTFLAGS} ${MOD} ${FIL} |& tee -a ${LOG}
  done

  if [[ ${ALL} -gt 0 ]]; then
    ${DIR}/mirabelle.sh -v ${VARV} ${ALLFLAGS} ${MOD} ${FIL} |& tee -a ${LOG}
  fi

  echo "-------- prune ${FIL} --------"
  egrep -i '(improv|other)' ${LOG}
  \rm ${LOG}
fi
