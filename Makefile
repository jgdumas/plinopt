####################################################################
# PLinOpt: a collection of C++ routines handling linear programs
# Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
####################################################################

OPTFLAGS=-Ofast -D RANDOM_TIES

CXXFLAGS += ${OPTFLAGS} `pkg-config linbox --cflags`
LOADLIBES+= `pkg-config linbox --libs |sed 's/-liml//;s/-lfplll//;s/-lflint//'`

EXE = sms2pretty matrix-transpose MMchecker 
EXE += transpozer optimizer sparsifier
EXE += factorizer
EXE += inplacer 
SRC=${EXE:%=%.cpp}

all: ${EXE}

clean:
	- \rm ${EXE}
