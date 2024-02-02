####################################################################
# PLinOpt: a collection of C++ routines handling linear programs
# Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
####################################################################

OPTFLAGS=-Ofast -D RANDOM_TIES

CXXFLAGS += ${OPTFLAGS} `pkg-config linbox --cflags`
LOADLIBES+= `pkg-config linbox --libs`

EXE = sms2pretty matrix-transpose MMchecker 
EXE += transpozer optimizer sparsifier
EXE += inplacer 
SRC=${EXE:%=%.cpp}

all: ${EXE}

clean:
	- \rm ${EXE}
