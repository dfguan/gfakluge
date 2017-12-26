CXX?=g++
CXXFLAGS:=-O3 -std=c++11 -g
LD_LIB_FLAGS=-L./
MERGE_OBJS = $(patsubst %.cpp, %.o, $(wildcard merge_*.cpp))
EXE=merge_gfa

all: $(EXE) 

$(EXE): $(MERGE_OBJS) lib 
	$(CXX) $(CXXFLAGS) $(LD_LIB_FLAGS) -o $@ $(MERGE_OBJS) -lgfakluge 

lib: gfakluge.o 
	ar -rs libgfakluge.a $<

gfakluge.o: gfakluge.cpp gfakluge.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

.PHONY: clean
clean:
	$(RM) -r *.o *.a $(EXE) 
