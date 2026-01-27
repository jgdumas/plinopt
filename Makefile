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
EXE += sms2pretty MMchecker PMchecker
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
check: ${BIN} pmcheck mmcheck opcheck slpcheck

mmcheck:
	./bin/MMchecker data/2x2x2_7_Strassen_{L,R,P}.sms
	./bin/MMchecker data/2x2x2_7_DPS-accurate_{L,R,P}.sms -m 513083
	./bin/MMchecker data/2x2x2_7_DPS-accurate_{L,R,P}.sms -r 1013 2 3
	./bin/MMchecker data/2x2x2_7_DPS-accurate-X_{L,R,P}.sms -P "X^2-3"

pmcheck:
	./bin/PMchecker data/4o4o8_Toom5_{L,R,P}.sms
	./bin/PMchecker data/4o4o8_Montgomery-13-58_{L,R,P}.sms
	./bin/PMchecker data/4o4o4_F243-11-44_{L,R,P}.sms -q 3 -P "1-X+X^5"
	./bin/PMchecker data/4o4o4_F243-Montgomery-13-42_{L,R,P}.sms -q 3 -P "X^5+X^4-X^3-X^2-1"


opcheck:
	./bin/GDT.sh

slpcheck:
	./bin/optimizer data/2x2x2_7_DPS-accurate_L.sms -E -N | ./bin/compacter | ./bin/SLPchecker -M data/2x2x2_7_DPS-accurate_L.sms
	./bin/FDT.sh

largecheck:
	./bin/MMchecker -b 5 data/32x32x32_15096_{L,R,P}.sms
	./bin/SLPchecker -M data/32x32x32_15096_L.s{ms,lp}
	./bin/SLPchecker -M data/32x32x32_15096_R.s{ms,lp}
	./bin/SLPchecker -M data/32x32x32_15096_P.s{ms,lp}
