CXX?=g++
CXXFLAGS:=-O3 -std=c++11 -g
LD_LIB_FLAGS=-L./
MERGE_OBJS = $(patsubst %.cpp, %.o, $(wildcard merge_*.cpp))
CONV_OBJS = convert_gfa.o
MERGE_EXE=merge_gfa
CONV_EXE=conv_gfa

all: $(MERGE_EXE) $(CONV_EXE) 

$(MERGE_EXE): $(MERGE_OBJS) lib 
	$(CXX) $(CXXFLAGS) $(LD_LIB_FLAGS) -o $@ $(MERGE_OBJS) -lgfakluge 

$(CONV_EXE): $(CONV_OBJS) lib 
	$(CXX) $(CXXFLAGS) $(LD_LIB_FLAGS) -o $@ $(CONV_OBJS) -lgfakluge 

lib: gfakluge.o 
	ar -rs libgfakluge.a $<

gfakluge.o: gfakluge.cpp gfakluge.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

.PHONY: clean
clean:
	$(RM) -r *.o *.a $(MERGE_EXE) $(CONV_EXE) 
