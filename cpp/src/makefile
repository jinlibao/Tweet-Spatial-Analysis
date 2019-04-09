#
# makefile
# Libao Jin, 2019-04-08 01:50
#

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -I$(MKL_INCLUDE) -L$(MKL_LIB) -I$(ARMADILLO_INCLUDE) -L$(ARMADILLO_LIB) -larmadillo -lmkl_core
BIN=../bin

RM = rm -rf

ALLPROGS = main

.DEFAULT_GOAL := all

all: $(ALLPROGS)

main: main.cpp apd.o
	$(CXX) $(CXXFLAGS) -o $(BIN)/$@ $?

clean:
	$(RM) $(ALLPROGS:=.dSYM) $(ALLPROGS:=.o) *.o

distclean: clean
	$(RM) $(addprefix $(BIN)/,$(ALLPROGS))
