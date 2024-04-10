####################################################################
# PLinOpt: a collection of C++ routines handling linear programs
# Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
####################################################################

OPTFLAGS=-Ofast

# OPTFLAGS += -D INPLACE_CHECKER		# adds Maple checks
# OPTFLAGS += -D VERBATIM_PARSING		# Verbose output
# OPTFLAGS += -D DEFAULT_RANDOM_LOOPS=30u	# Default # of loops
# OPTFLAGS += -D COEFFICIENT_SEARCH=20u		# Default # sparsifier coeffs

CXXFLAGS += -D RANDOM_TIES -D INPLACE_CHECKER

CXXFLAGS += ${OPTFLAGS} `pkg-config linbox --cflags`
LOADLIBES+= `pkg-config linbox --libs |sed 's/-liml//;s/-lfplll//;s/-lflint//'`

EXE = sms2pretty matrix-transpose MMchecker
EXE += transpozer compacter PMchecker
EXE += optimizer
EXE += sparsifier factorizer
EXE += inplacer
SRC=${EXE:%=%.cpp}

all: ${EXE}

clean:
	- \rm ${EXE}


check: ${EXE}
	./FDT.sh