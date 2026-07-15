#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Testing linear combination of variables of SLPS
# ==========================================================================

DIR=`dirname $0`
MOD=""
COE=2
LVL=4
VEC=""

while [[ $# -gt 0 ]]; do
  case $1 in
    -c|-n|--coeffs)
    COE="$2"
    shift # past argument
    shift # past value
    ;;
    -v|--vector)
    VEC="$(printf ' %q' $2)" # use printf %q, to pass the string as is to DEPND
    shift # past argument
    shift # past value
    ;;
    -l|--level)
    LVL="$2"
    shift # past argument
    shift # past value
    ;;
    -r|-m|-q|--modular)
    MOD=" -q $2"
    shift # past argument
    shift # past value
    ;;
    -h|--h|-help|--help|-*|--*)
    echo "Usage: $0 [-c|-q|-l #] [-v \"# ... #\"] P.slp"
    exit 1
    ;;
    *)
    FIL=$1
    shift # past argument
    ;;
    esac
done

GRE='\033[1;32m'
RED='\033[0;41m'
BLU='\033[0;96m'
NC='\033[0m'    # No Color

#############################################################
## PLinOpt programs

SLPCHK="${DIR}/SLPchecker${MOD}"
DEPND="${DIR}/dependency${MOD}"
#############################################################
## Show here-strings

function Show() {
    local NAM=$1
    >&2 echo -e "# ----- BEG: ${NAM} -----\n${!NAM}\n# ----- END: ${NAM} -----"
}
function ShowQuietUniq() {
    local NAM=$1
    echo -e "${!NAM}"| egrep -v '($^)'| sort -k1,1 -t';' --stable --unique
}

#############################################################
## Variables to optimize

OVARS=$(sed -r "s/([:=+-])/ /g" ${FIL}|awk '{var=$1;nbm=NF;b=nbm>3;if (!b){gsub(/[^*]/,""); b=length}; if (b) print "<"var","nbm"> "}')
LOARS=(`echo "${OVARS}"`)
>&2 echo "# [CHTRS] VARS>3 (${#LOARS[@]}): ${LOARS[@]}"

#############################################################
## Temporary file names

NAM="chartreuse"

#############################################################
## Define output and input variables of the subprogram
##     together with replacement names

CHARS=(`sed 's/:=/ /;s/+/ /g;s/:=/ /;s/+/ /g;s/-/ /g;s/;.*/ /;s/\*[0-9]* / /g;s/\/[0-9]* / /g;s/)//g;s/(//g' ${FIL} | tr ' ' '\n' | sed 's/[0-9]//g;/^$/d' | sort -u| tr '\n' ' '`)
# >&2 echo "# CHARS: ${CHARS[@]}"

NCHAR="a"
while grep -q ${NCHAR} <<< ${CHARS[@]}; do
    NCHAR=`echo ${NCHAR} | tr "a-z" "b-za"`
done
# >&2 echo "# NCHAR: ${NCHAR}"

CHARS+=(${NCHAR})
OCHAR=${NCHAR}
while grep -q ${OCHAR} <<< ${CHARS[@]}; do
    OCHAR=`echo ${OCHAR} | tr "a-z" "b-za"`
done
# >&2 echo "# OCHAR: ${OCHAR}"


#############################################################
## Build the subprogram in input 'i' and output 'o'

BOD=$(sed "s/i/${NCHAR}/g;s/o/${OCHAR}/g" ${FIL})
# Show BOD
TSDO=$(cat <<< "${BOD}" | cut -d':' -f1 | sort -r| awk 'BEGIN {s=0} {print "s/"$1"/o"s"/g";s++}' |tac|tr '\n' ';')
# Show TSDO

SINP=$(sed 's/:=.*/\[\^0-9\]|/g' <<< ${BOD} | tr '\n' ' '|sed 's/ //g;s/|$//')
INP=`echo "${SINP}"`
# >&2 echo "# INP: ${INP}"

HEA=$(sed 's/.*:=//;s/+/ /g;s/-/ /g;s/;.*/ /;s/\*[0-9]*/ /g;s/\/[0-9]*/ /g;s/)//g;s/(//g' <<< "${BOD}" | tr -s '[:space:]' | tr ' ' '\n'|sort -u| sed 's/$/;/'|egrep -v "(^;$|^i|${INP})"|sed 's/;.*//'|awk 'BEGIN {s=0} {print $1":=i"s";";s++}')

REST=$(echo -e "${HEA}")
REST+=$(echo -e "\n${BOD}")
REST+=$(echo -e "${BOD}" | cut -d':' -f1 | sort -r| awk 'BEGIN {s=0;print} {print "o"s":="$1";";s++}')
# Show REST


#############################################################
## Replacements to go back to original variable names

SDO=$(sed 's/s\/\([^\/]*\)\/\([^\/]*\)\/g/s\/\2\/\1\/g/g' <<< "${TSDO}")
SDO=${SDO}"s/${OCHAR}/o/g;s/${NCHAR}/i/g"
# Show SDO

SDI=$(tac <<< "${HEA}" | sed 's/:=/ /;s/;.*//' | awk '{print "s/"$2"/"$1"/g"}'|tr '\n' ';')
# Show SDI

#############################################################
## Compute the dependencies
DPNDL=$(${SLPCHK} <<< "${REST}" | ${DEPND} -c ${COE} -v "${VEC}" -l ${LVL})

function SortLine() {
    while read line
      do
      echo "$line" | sed 's/[+-]/ &/g;s/;//' | tr ' ' '\n' | sort | tr '\n' ' ' | sed -E 's/ //g;s/$/;\n/'
    done
}

## Rewrite back the original variable names and sort the dependencies
COMBR=$(echo ${DPNDL} | sed "${SDI};${SDO}"|tr ' ' '\n')
Combinations=$(SortLine <<< "${COMBR}" |sort -u)
>&2 echo "# [CHRTS] All linear combinations:"
Show Combinations

## Write and sort the dependencies within the original program
RELPL=$(sed -r 's/(.*):=(-.*);/\2-\1;/;s/(.*):=(.*);/+\2-\1;/' "${FIL}")
RELPS=$(SortLine <<< "${RELPL}" | awk '{orig=$0;gsub(/\+/,"PLUSPLUS");gsub(/\-/,"+");gsub(/PLUSPLUS/,"-");print orig; print}' |sort -u)

## Select only not known dependencies
ONLYN=$(comm -23 <(cat <<< "${Combinations}") <(cat <<< "${RELPS}"))

## Select improving dependencies (counting only add/sub)

LESAD="${BLU}\t\t/!\\\ less additions /!\\\\${NC}"
IMPRO="${GRE}\t\t/!\\\  IMPROVEMENTS  /!\\\\${NC}"

TOTAL=()
for ovr in ${OVARS[@]}; do
    OVRC=(`echo "${ovr}" | sed 's/[<,>]/ /g'`)
    OFND=$(egrep "(${OVRC[0]}[^0-9])" <<< "${ONLYN}")
    SFND=$(awk -v onr="${OVRC[1]}" -v onv="${OVRC[0]}" -v lea="${LESAD}" -v imp="${IMPRO}" '{orig=$0;gsub("[^+-]",""); l=length+1;if (length>0 && l<=onr) {$0=orig;gsub(/[^*/]/,""); if (length) {msg=lea} else {msg=imp}; print orig,"\t# "onv" "onr" --> "(l-1),msg}}' <<< "${OFND}")
#     Show SFND
    TOTAL+='\n'
    TOTAL+=${SFND}
done

ShowQuietUniq TOTAL
