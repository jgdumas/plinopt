####################################################################
# PLinOpt: a collection of C++ routines handling linear programs
# Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
####################################################################

OPTFLAGS = -O3 -ffast-math

# OPTFLAGS += -D INPLACE_CHECKER	# adds Maple checks (needed by ./bin/TTP.sh)
# OPTFLAGS += -D VERBATIM_PARSING=1	# Verbose output
# OPTFLAGS += -D DEFAULT_RANDOM_LOOPS=30u	# Default # of loops
# OPTFLAGS += -D COEFFICIENT_SEARCH=20u		# Default # sparsifier coeffs
# OPTFLAGS += -D DENSITY_OPTIMIZATION	# Non-random optimizer
#######

RNDFLAGS = -D RANDOM_TIES			# Default randomized search

#######
CXXFLAGS += ${OPTFLAGS} ${RNDFLAGS} -I`pwd`/include/ `pkg-config linbox --cflags`
LOADLIBES+= `pkg-config linbox --libs |sed 's/-liml//;s/-lfplll//;s/-lflint//'`

#######

EXE  = optimizer
EXE += sparsifier factorizer
EXE += inplacer trilplacer
EXE += transpozer compacter SLPchecker
EXE += sms2pretty MMchecker
EXE += matrix-transpose columns-swap negater

SRC=${EXE:%=src/%.cpp}

BIN=${EXE:%=bin/%}

#######

all: ${BIN}

bin/%: src/%.cpp
	$(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	- \rm ${BIN}


SHELL=/bin/bash
check: ${BIN} mmcheck slpcheck

mmcheck:
	./bin/MMchecker data/2x2x2_7_Strassen_{L,R,P}.sms
	./bin/MMchecker data/2x2x2_7_DPS-accurate_{L,R,P}.sms -m 513083
	./bin/MMchecker data/2x2x2_7_DPS-accurate_{L,R,P}.sms -r 1013 2 3
	./bin/MMchecker data/2x2x2_7_DPS-accurate-X_{L,R,P}.sms -P "X^2-3"

slpcheck:
	./bin/optimizer data/2x2x2_7_DPS-accurate_L.sms -E -N | compacter | SLPchecker -M data/2x2x2_7_DPS-accurate_L.sms
	./bin/FDT.sh

largecheck:
	./bin/MMchecker -b 5 data/32x32x32_15096_{L,R,P}.sms
	./bin/SLPchecker -M data/32x32x32_15096_L.s{ms,lp}
	./bin/SLPchecker -M data/32x32x32_15096_R.s{ms,lp}
	./bin/SLPchecker -M data/32x32x32_15096_P.s{ms,lp}
