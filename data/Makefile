####################################################################
# PLinOpt: a collection of C++ routines handling linear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
####################################################################

# Change FILE for algorithm prefix to check
FILE=4x4x4_48_rational
#PM=1
##################"
MATS=$(shell echo ${FILE}_{L,R,P}.sms)
##################"
PLODIR=../bin
ifdef PM
  CHKPRG=${PLODIR}/PMchecker
else
  CHKPRG=${PLODIR}/MMchecker
endif

SHELL:=/bin/bash
PREFIX=${patsubst %_L.sms,%,${firstword ${notdir ${MATS}}}}

OPT=1000
DIM= #-f 62 -a
#MOD=

ALGO=${MATS:%.sms=%.slp}
CHKS=${MATS:%.sms=%.chk}

all: check

%.sms:%.slp
	${PLODIR}/SLPchecker ${MOD} $< > $@

count: ${ALGO}
	@ ${PLODIR}/OpCount.sh ${ALGO}

%.chk:%.sms %.slp
	@ ${PLODIR}/SLPchecker ${MOD} -M $^

check: docheck count

docheck: ${CHKS}
	@ ${CHKPRG} ${MOD} ${MATS}

prune_%.log: ${FILE}_%.slp
	(time ${PLODIR}/prune.sh -O ${OPT} $< ) &> $@

prune: $(shell echo prune_{L,R,P}.log)

