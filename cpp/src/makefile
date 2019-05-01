#
# makefile
# Libao Jin, 2019-04-08 01:50
#

CXX = mpicxx

# Intel MKL for mac
#CXXFLAGS=-std=c++11 -O2 -Wall -m64 -I$(MKL_INCLUDE) -I$(ARMADILLO_INCLUDE) -I$(OPENMPIPATH)/include
#LINK=-I$(OPENMPIPATH)/lib -L$(MKL_LIB) -lmkl_rt -lpthread -lm -ldl -L$(ARMADILLO_LIB) -larmadillo -lmpi

# For Linux
# Intel MKL 1
CXXFLAGS=-std=c++11 -O2 -Wall -m64 -I$(MKL_INCLUDE) -I$(ARMADILLO_INCLUDE) -I$(CPATH)
LINK=-L$(MKL_LIB) -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl -L$(ARMADILLO_LIB) -larmadillo -lmpi

# Intel MKL 2
#CXXFLAGS=-std=c++11 -O2 -Wall -DMKL_ILP64 -m64 -I$(MKL_INCLUDE) -I$(ARMADILLO_INCLUDE) -I$(CPATH)
#LINK=-L$(MKL_LIB) -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -L$(ARMADILLO_LIB) -larmadillo -lmpi

# OpenBLAS
#CXXFLAGS = -std=c++11 -O2 -Wall -m64 -I$(ARMADILLO_INCLUDE) -I/home/ljin1/.linuxbrew/include -I$(CPATH) -L/home/ljin1/.linuxbrew/lib -L$(ARMADILLO_LIB) -larmadillo -lopenblas -lm -lmpi

BIN=../bin

RM = rm -rf

ALLPROGS = main

.DEFAULT_GOAL := all

all: $(ALLPROGS)

main: main.cpp ellipse.o matmul.o tweets_spatial_analysis.o
	$(CXX) $? $(CXXFLAGS) -o $(BIN)/$@ $(LINK)

clean:
	$(RM) $(ALLPROGS:=.dSYM) $(ALLPROGS:=.o) *.o

distclean: clean
	$(RM) $(addprefix $(BIN)/,$(ALLPROGS))
