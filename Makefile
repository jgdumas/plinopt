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
EXE += sms2pretty matrix-transpose columns-swap MMchecker
EXE += PMMchecker

SRC=${EXE:%=src/%.cpp}

BIN=${EXE:%=bin/%}

#######

all: ${BIN}

bin/%: src/%.cpp
	$(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	- \rm ${BIN}


check: ${BIN}
	./bin/FDT.sh
