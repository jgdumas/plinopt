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
MOD=""	# modular search
LIM=""	# subset of variables containing LIM
ALL=0	# Mirabelle on all variables at once, (with reduced -O $2)
SQR=0	# Try all *pairs* of variables

while [[ $# -gt 0 ]]; do
  case $1 in
    -O|--Optflags)
    OPTFLAGS="-O $2"
    shift # past argument
    shift # past value
    ;;
    -r|-q|-m|--modular)
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
    -p|-s|--pairs)
    SQR=1
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
    echo "Usage: $0 [-O|-a #] [-q|-r|-m #] [-l pattern] [-n|-nl] [-p] f.slp"
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
  if [[ ${SQR} -gt 0 ]]; then
    for var in ${VARS[@]}; do
      for rav in ${VARS[@]}; do
	${DIR}/mirabelle.sh -v $var -v $rav ${OPTFLAGS} ${MOD} ${FIL} |& tee -a ${LOG}
      done
    done
  else
    for var in ${VARS[@]}; do
      ${DIR}/mirabelle.sh -v $var ${OPTFLAGS} ${MOD} ${FIL} |& tee -a ${LOG}
    done
  fi

  if [[ ${ALL} -gt 0 ]]; then
    ${DIR}/mirabelle.sh -v ${VARV} ${ALLFLAGS} ${MOD} ${FIL} |& tee -a ${LOG}
  fi

  echo "-------- prune ${FIL} --------"
  grep -i less ${LOG}
  grep -i impr ${LOG}
  \rm ${LOG}
fi
