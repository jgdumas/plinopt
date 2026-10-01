####################################################################
# PLinOpt: a collection of C++ routines handling linear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
####################################################################

# Change FILE for algorithm prefix to check
FILE=4x4x4_48_204
ALT=16
#PM=1
##################"
MATS=$(shell echo ${FILE}_{L,R,P}.sms)
ALTS=$(shell echo ${FILE}-${ALT}ALT_{L,R,P}.sms)
COBS=$(shell echo ${FILE}-${ALT}CoB_{L,R,P}.sms)
##################"
PLODIR=../bin
ifdef PM
  CHKPRG=${PLODIR}/PMchecker
else
  CHKPRG=${MMMCHK}
endif

MMMCHK=${PLODIR}/MMchecker
ALTCHK=${PLODIR}/ALTchecker.sh
SLPCHK=${PLODIR}/SLPchecker
OPSCNT=${PLODIR}/OpCount.sh
PRUNER=${PLODIR}/prune.sh
CHRTRS=${PLODIR}/chartreuse.sh


SHELL:=/bin/bash
PREFIX=${patsubst %_L.sms,%,${firstword ${notdir ${MATS}}}}

OPT=1000
DIM= #-f 62 -a
#MOD=

ALGO=${MATS:%.sms=%.slp}
ALGA=${ALTS:%.sms=%.slp}
ALGC=${COBS:%.sms=%.slp}
CHKS=${MATS:%.sms=%.chk}

all: check

%.sms:%.slp
	${SLPCHK} ${MOD} $< > $@

count: ${ALGO}
	@ ${OPSCNT} ${ALGO}

%.chk:%.sms %.slp
	@ ${SLPCHK} ${MOD} -M $^

check: docheck count

docheck: ${CHKS}
	@ ${CHKPRG} ${MOD} ${MATS}

prune_%.log: ${FILE}_%.slp
	(time ${PRUNER} -O ${OPT} $< ) &> $@

prune: $(shell echo prune_{L,R,P}.log)

alt: ${MATS} ${ALGO} ${ALTS} ${COBS} ${ALGA} ${ALGC}
	$(foreach mat, L R, ${ALTCHK} ${FILE}-${ALT}{CoB,ALT}_${mat}.slp -M ${FILE}_${mat}.sms;)
	${ALTCHK} ${FILE}-${ALT}{ALT,CoB}_P.slp -M ${FILE}_P.sms
	${OPSCNT} ${ALGC}
	${OPSCNT} ${ALGA}