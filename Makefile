####################################################################
# PLinOpt: a collection of C++ routines for straight-line programs
# Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
####################################################################

OPTFLAGS = -O3 -ffast-math

# CTAD, std::gcd: g++ defaults to C++17, Apple clang to C++14
STDFLAGS = -std=c++17

############################
# CXXFLAGS += -D DEFAULT_RANDOM_LOOPS=30u	# Default # of loops
# CXXFLAGS += -D COEFFICIENT_SEARCH=20u		# Default # sparsifier coeffs

# CXXFLAGS += -D VERBATIM_PARSING=1	# Verbose output
# CXXFLAGS += -D INPLACE_CHECKER	# Maple checks (needed by ./bin/TTP.sh)

# CXXFLAGS += -D DENSITY_OPTIMIZATION	# Non-random optimizer

# CXXFLAGS += -D KERNEL_FULL_IDENTITY	# Identity goals in kernel method

# CXXFLAGS += -D ACTION_FULL_PLUQ	# Random PLUQ action
# CXXFLAGS += -D ACTION_HOUSEHOLDER	# Random Householder action
# CXXFLAGS += -D NEGATER_ONLY_SIGN	# Disable common divisors to matrix P
############################

RNDFLAGS = -D RANDOM_TIES		# Default randomized search

############################
# OpenMP: Apple clang has no -fopenmp; it uses Homebrew's libomp
#         (brew install libomp)
ifeq ($(shell uname -s),Darwin)
  LIBOMP  := $(shell brew --prefix libomp 2>/dev/null || echo /opt/homebrew/opt/libomp)
  OMPFLAGS = -Xpreprocessor -fopenmp -I${LIBOMP}/include
  OMPLIBS  = -L${LIBOMP}/lib -lomp
else
  OMPFLAGS = -fopenmp
endif
# clang rejects givaro's specialization of std::make_signed
#       (recint/rrint.h)
ifneq (,$(findstring clang,$(shell ${CXX} --version 2>/dev/null)))
  WNOFLAGS = -Wno-invalid-specialization
endif

############################

CXXFLAGS += ${OPTFLAGS} ${STDFLAGS} ${OMPFLAGS} ${WNOFLAGS} ${RNDFLAGS} -I`pwd`/include/ `pkg-config linbox --cflags`
LOADLIBES+= `pkg-config linbox --libs |sed 's/-liml//;s/-lfplll//;s/-lflint//'` ${OMPLIBS}

############################

EXE  = optimizer
EXE += sparsifier factorizer
EXE += inplacer trilplacer
EXE += transpozer compacter SLPchecker
EXE += sms2pretty MMchecker PMchecker
EXE += matrix-transpose columns-swap negater orbiter
EXE += growthfactor dependency

SRC=${EXE:%=src/%.cpp}

BIN=${EXE:%=bin/%}

############################

all: ${BIN}

bin/%: src/%.cpp $(wildcard include/*.inl)
	$(LINK.cpp) $< $(LOADLIBES) $(LDLIBS) -o $@

clean:
	- \rm ${BIN}


############################

SHELL=/bin/bash
check: ${BIN} trspcheck modcheck pmcheck mmcheck altcheck opcheck slpcheck

mmcheck: ./bin/MMchecker
	$< data/2x2x2_7_Strassen_{L,R,P}.sms
	$< data/2x2x2_7_DPS-accurate_{L,R,P}.sms -m 513083
	$< data/2x2x2_7_DPS-accurate_{L,R,P}.sms -r 1013 2 3
	$< data/2x2x2_7_DPS-accurate-X_{L,R,P}.sms -P "X^2-3"

pmcheck: ./bin/PMchecker
	$< data/4o4o8_Toom5_{L,R,P}.sms
	$< data/4o4o8_Montgomery-13_{L,R,P}.sms
	$< data/4o4o4_F243-11-44_{L,R,P}.sms -q 3 -P "1-X+X^5"
	$< data/4o4o4_F243-Montgomery-13-42_{L,R,P}.sms -q 3 -P "X^5+X^4-X^3-X^2-1"

modcheck: AALG=2x2x2_7_Strassen
modcheck: ./bin/sparsifier
	$< data/${AALG}_L.sms
	$(foreach mod, 2 3 5, $< -q ${mod} data/${AALG}_L.sms;)

trspcheck: ./bin/compacter ./bin/transpozer ./bin/SLPchecker data/test-prg.sms
	./bin/transpozer data/test-prg.slp | ./bin/compacter | ./bin/transpozer | ./bin/compacter | ./bin/SLPchecker -M data/test-prg.sms


%.sms:%.slp
	./bin/compacter $< | ./bin/SLPchecker > $@

opcheck: ./bin/GDT.sh ${BIN}
	$<

slpcheck: ./bin/FDT.sh ./bin/optimizer
	./bin/optimizer data/2x2x2_7_DPS-accurate_L.sms -E -N | ./bin/compacter | ./bin/SLPchecker -M data/2x2x2_7_DPS-accurate_L.sms
	$<

altcheck: AALG=4x4x4_48_204
altcheck: ./bin/ALTchecker.sh
	$(foreach dim, 16 24, $(foreach mat, L R, $< data/${AALG}-${dim}{CoB,ALT}_${mat}.slp -M data/${AALG}_${mat}.sms;))
	$(foreach dim, 16 24 32, $< data/${AALG}-${dim}{ALT,CoB}_P.slp -M data/${AALG}_P.sms;)

largecheck: AALG=32x32x32_15096
largecheck: ./bin/SLPchecker ./bin/MMchecker
	./bin/MMchecker -b 5 data/${AALG}_{L,R,P}.sms
	$(foreach mat, L R P, $< -M data/${AALG}_${mat}.s{ms,lp};)
