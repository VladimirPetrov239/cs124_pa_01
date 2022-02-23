.PHONY: all clean

CXX=g++
CXXFLAGS=-std=c++17 -Wall -pedantic

all: randmst

randmst: src/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ -Iinclude $<

clean:
	rm -rf randmst