#
# makefile
# Libao Jin, 2019-04-08 01:50
#

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -I$(MKL_INCLUDE) -L$(MKL_LIB) -I$(ARMADILLO_INCLUDE) -L$(ARMADILLO_LIB) -larmadillo -lmkl_core
BIN=../bin

RM = rm -rf

ALLPROGS = example build_distance_matrix

.DEFAULT_GOAL := all

all: $(ALLPROGS)

build_distance_matrix: build_distance_matrix.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/$@ $<

example: example.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/$@ $<

clean:
	$(RM) $(ALLPROGS:=.dSYM) $(ALLPROGS:=.o)

distclean: clean
	$(RM) $(addprefix $(BIN)/,$(ALLPROGS))
