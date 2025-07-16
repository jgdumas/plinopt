####################################################################
# PLinOpt: a collection of C++ routines handling linear programs
# Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
####################################################################

OPTFLAGS = -O3 -ffast-math

# OPTFLAGS += -D INPLACE_CHECKER		# adds Maple checks
# OPTFLAGS += -D VERBATIM_PARSING		# Verbose output
# OPTFLAGS += -D DEFAULT_RANDOM_LOOPS=30u	# Default # of loops
# OPTFLAGS += -D COEFFICIENT_SEARCH=20u		# Default # sparsifier coeffs

#######

CXXFLAGS += -D RANDOM_TIES # -D INPLACE_CHECKER

CXXFLAGS += ${OPTFLAGS} -I`pwd`/include/ `pkg-config linbox --cflags`
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
check: ${BIN}
	./bin/MMchecker data/2x2x2_7_DPS-accurate_{L,R,P}.sms -m 513083
	./bin/MMchecker data/2x2x2_7_Strassen_{L,R,P}.sms
	./bin/MMchecker data/2x2x2_7_DPS-accurate-X_{L,R,P}.sms -P "X^2-3"
	./bin/MMchecker data/2x2x2_7_DPS-accurate_{L,R,P}.sms -r 1013 2 3
	./bin/FDT.sh
