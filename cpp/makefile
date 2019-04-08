#
# makefile
# Libao Jin, 2019-04-08 01:50
#

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -I/opt/intel/mkl/include -L/opt/intel/mkl/lib  -I/Users/libao/.armadillo/armadillo-9.300.2-mkl/include -L/Users/libao/.armadillo/armadillo-9.300.2-mkl/lib -larmadillo -lmkl_core

RM = rm -rf

ALLPROGS = example build_distance_matrix

.DEFAULT_GOAL := all

all: $(ALLPROGS)

build_distance_matrix: build_distance_matrix.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

example: example.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	$(RM) $(ALLPROGS:=.dSYM) $(ALLPROGS:=.o)

distclean: clean
	$(RM) $(ALLPROGS)
