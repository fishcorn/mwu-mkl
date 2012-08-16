CC=gcc
CPP=g++
CFLAGS=-c -Wall -fPIC -m32 -O3
CXXFLAGS=
LDFLAGS=-shared -m32
CSOURCES=
CXXSOURCES=											\
	mwu_main.cpp
OBJECTS=$(CXXSOURCES:%.cpp=%.o) $(CSOURCES:%.c=%.o)
LIBRARYNAME=libmwu.so

library: $(LIBRARYNAME) 

$(LIBRARYNAME): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

.cpp.o:
	$(CPP) $(CFLAGS) $(CXXFLAGS) $< -o $@
