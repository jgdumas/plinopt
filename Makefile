####################################################################
# PLinOpt: a collection of C++ routines handling linear programs
# Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
####################################################################

OPTFLAGS=-Ofast

CXXFLAGS += ${OPTFLAGS} `pkg-config linbox --cflags`
LOADLIBES+= `pkg-config linbox --libs`

EXE=sms2pretty matrix-transpose transpozer optimizer inplacer
SRC=${EXE:%=%.cpp}

all: ${EXE}

clean:
	- \rm ${EXE}