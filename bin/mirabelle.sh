#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Wrapper optimizing subprograms that contain the given variables
# ==========================================================================

DIR=`dirname $0`
OPTFLAGS=""
VARS=()
MOD=""

while [[ $# -gt 0 ]]; do
  case $1 in
    -O|--Optflags)
    OPTFLAGS=" -O $2"
    shift # past argument
    shift # past value
    ;;
    -v|--var)
    VARS+=(`echo "$2" | sed 's/|/ /g'`)
    shift # past argument
    shift # past value
    ;;
    -r|-m|--modular)
    MOD=" -q $2"
    shift # past argument
    shift # past value
    ;;
    -h|--h|-help|--help|-*|--*)
    echo "Usage: $0 [-O #] [-v var[|var]*] P.slp"
      exit 1
      ;;
    *)
    FIL=$1
    shift # past argument
    ;;
    esac
done

SLPCHK="${DIR}/SLPchecker${MOD}"
OPTMZR="${DIR}/optimizer${OPTFLAGS}${MOD}"
MATTRP="${DIR}/matrix-transpose"
TRSPZR="${DIR}/transpozer"

GRE='\033[0;32m'
RED='\033[0;41m'
BLU='\033[0;36m'
NC='\033[0m'    # No Color

#############################################################
## Temporary file names

VAR=${VARS[0]}
VARP=`echo ${VARS[@]}|sed 's/ /|/g'`
# echo "VARS: ${VARS[@]}"

NAM="mirabelle_${VAR}"
HEA="${NAM}_inp-$$.slp"
BOD="${NAM}_bod-$$.slp"
OPT="${NAM}_opt-$$.slp"
COM="${NAM}-$$.com"
RUN="${NAM}-$$.run"
SDO="${NAM}_sdo-$$.sed"
SDI="${NAM}_sdi-$$.sed"
RES="${NAM}-$$.slp"
FND="${NAM}-$$.log"

#############################################################
## Recover all related variables to form the subprogram

FROS="("
for vari in ${VARS[@]}; do
    FRO=`grep ${vari} ${FIL}| head -1 | sed 's/.*:=+//;s/.*:=-//;s/.*:=//;s/+/ /g;s/-/ /g;s/;.*/ /;s/\*[0-9]* / /g;s/\/[0-9]* / /g;s/ /\[\^0-9\]|/g'`
    FROS="${FROS}${FRO}"
done
FROS="${FROS}${VARP})"
# echo "FROS: ${FROS}"

egrep ${FROS} ${FIL} > ${BOD}

#############################################################
## Define output and input variables of the subprogram
##     together with replacement names

CHARS=(`sed 's/:=/ /;s/+/ /g;s/:=/ /;s/+/ /g;s/-/ /g;s/;.*/ /;s/\*[0-9]* / /g;s/\/[0-9]* / /g;s/)//g;s/(//g' ${BOD} | tr ' ' '\n' | sed 's/[0-9]//g;/^$/d' | sort -u| tr '\n' ' '`)
# echo "CHARS: ${CHARS[@]}"

NCHAR="a"
while grep -q ${NCHAR} <<< ${CHARS[@]}; do
    NCHAR=`echo ${NCHAR} | tr "a-z" "b-za"`
done
# echo "NCHAR: ${NCHAR}"

CHARS+=(${NCHAR})
OCHAR=${NCHAR}
while grep -q ${OCHAR} <<< ${CHARS[@]}; do
    OCHAR=`echo ${OCHAR} | tr "a-z" "b-za"`
done
# echo "OCHAR: ${OCHAR}"


sed -i "s/i/${NCHAR}/g;s/o/${OCHAR}/g" ${BOD}


INP=`sed 's/:=.*/\[\^0-9\]|/g' ${BOD} | tr '\n' ' '|sed 's/ //g;s/|$//'`
# echo "INP: ${INP}"

sed 's/.*:=//;s/+/ /g;s/-/ /g;s/;.*/ /;s/\*[0-9]*/ /g;s/\/[0-9]*/ /g;s/)//g;s/(//g' ${BOD} | tr -s '[:space:]' | tr ' ' '\n'|sort -u| sed 's/$/;/'|egrep -v "(^;$|^i|${VAR}|${INP})"|sed 's/;.*//'|awk 'BEGIN {s=0} {print $1":=i"s";";s++}' > ${HEA}

cat ${HEA} ${BOD} > ${RES}


egrep -v "(^${VAR})" ${BOD} | cut -d':' -f1 | sort -r| awk 'BEGIN {s=0} {print "s/"$1"/o"s"/g";s++}' |tac > ${SDO}
egrep -v "(^${VAR})" ${BOD} | cut -d':' -f1 | sort -r| awk 'BEGIN {s=0} {print "o"s":="$1";";s++}' >> ${RES}

###### sed -i -f ${SDO} ${RES}
sed -i 's/s\/\([^\/]*\)\/\([^\/]*\)\/g/s\/\2\/\1\//g' ${SDO}
echo "s/${OCHAR}/o/g;s/${NCHAR}/i/g" >> ${SDO}

#############################################################
## Compute original subprogram number of operations

BEF=(`(${SLPCHK} ${RES} |& egrep '(additions|multiplications)' | sed 's/\x1b\[[0-9;]*[a-zA-Z]//g'| awk '{print $2}') 2> /dev/null`)
#echo ${BEF[*]}

#############################################################
## Function comparing the subprogram and an optimized version
function Compare() {
  AFT=(`egrep '(additions|multiplications)' ${COM} | tail -2 | sed 's/\x1b\[[0-9;]*[a-zA-Z]//g'| awk '{print $2}'`)
#echo ${AFT[*]}

  DIF=$((BEF[0]+BEF[1]-AFT[0]-AFT[1]))
#echo $DIF

  if [[ "$DIF" -gt 0 ]]; then
      >&2 echo -e "${GRE}> ${AFT[*]}\t\t/!\\ IMPROVEMENT /!\ ${NC}"

      tac ${HEA} | sed 's/:=/ /;s/;.*//' | awk '{print "s/"$2"/"$1"/g;"}' > ${SDI}
      sed "s/${OCHAR}/o/g;s/${NCHAR}/i/g" ${BOD} > ${FND}
      uniq ${COM} &>> ${FND}
      ((compacter ${OPT} | egrep -v '(:=0;)' | sed -f ${SDI} | sed -f ${SDO}) >> ${FND}) 2> /dev/null
  else
      ADD=$((BEF[0]-AFT[0]))
      MUL=$((BEF[1]-AFT[1]))
      MSG=
      if [[ "$ADD" -gt 0 ]]; then
	  MSG="additions"
      fi
      if [[ "$MUL" -gt 0 ]]; then
	  MSG="multiplications"
      fi
      if [[ "$MSG" != "" ]]; then

	  >&2 echo -e "== ${BLU}${AFT[*]}\t\t less ${MSG} ...${NC}"

	  tac ${HEA} | sed 's/:=/ /;s/;.*//' | awk '{print "s/"$2"/"$1"/g;"}' > ${SDI}
	  sed "s/${OCHAR}/o/g;s/${NCHAR}/i/g" ${BOD} > ${FND}
	  uniq ${COM} &>> ${FND}
	  ((compacter ${OPT} | egrep -v '(:=0;)' | sed -f ${SDI} | sed -f ${SDO}) >> ${FND}) 2> /dev/null
      else
	  if [[ "$DIF" -eq 0 ]]; then
	      >&2 echo "== ${AFT[*]}"
	  else
	      >&2 echo "≤ ${AFT[*]}"
	  fi
      fi
  fi
}
#############################################################


#############################################################
## Optimize program

echo -n "${VARS[@]}: ${BEF[*]} "
echo "((${SLPCHK} ${RES} | ${OPTMZR}) > ${OPT}) 2> ${COM}" > ${RUN}
chmod +x ${RUN}
${RUN}
Compare

#############################################################
## Optimize its transposition

echo -n "${VARS[@]}t: ${BEF[*]} "
echo "((${SLPCHK} ${RES} | ${MATTRP} | ${OPTMZR} | ${TRSPZR} ) > ${OPT}) 2> ${COM}" > ${RUN}
chmod +x ${RUN}
${RUN}
FND="${NAM}t-$$.log"
Compare


#############################################################
## Clean-up tyemporary files

\rm -rf ${RES} ${HEA} ${BOD} ${SDO} ${SDI} ${OPT} ${COM} ${RUN}
