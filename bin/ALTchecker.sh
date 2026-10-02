#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# Tests: CoB/ALT coherency
#        --> creates the matrix computed by the composition of $2 o $1
#   If   C.sms is computed by LeftSLP=$1 and A.sms by RightSLP=$2
#   Then SLPchecker will pass iff (LeftSLP;RightSLP) passes (C.sms . A.sms)
# ==========================================================================

#example: ./bin/ALTchecker.sh data/4x4x4_48_204-16{CoB,ALT}_R.slp -M data/4x4x4_48_204_R.sms

DIR=`dirname $0`
SLPCHK="${DIR}/SLPchecker"
OPSCNT="${DIR}/OpCount.sh"

if [[ $# -lt 2  ]]; then
  echo "Usage: $0 l.slp r.slp [-q #] [-M file.sms]"
  echo "  (all parameters after l.slp r.slp are passed to ${SLPCHK})"
  echo "  -q #: search modulo"
  echo "  -M file.sms: checks the resulting SLP with file.sms"
  exit 1
fi

LeftSLP=$1
RightSLP=$2
shift
shift

function SLPSMSchk() {
    local SLP=$1
    local SMS=`dirname ${SLP}`'/'`basename ${SLP} .slp`'.sms'
    if [ -e "${SMS}" ]; then
	${SLPCHK} ${SLP} -M ${SMS}
    else
	${OPSCNT} ${SLP}
    fi
}

SLPSMSchk ${LeftSLP}
SLPSMSchk ${RightSLP}

cat <((compacter ${LeftSLP} 2>/dev/null) | sed 's/o/j/g') <((compacter ${RightSLP} 2>/dev/null) | sed 's/i/j/g' ) | ${SLPCHK} $*
