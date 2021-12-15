CXX := `root-config --cxx`
ROOTCLING=rootcling
MARCH := `root-config --arch`
LD:=$(CXX)
SRC=./src/

UNAME := $(shell uname)

CFLAGS += $(shell root-config --ldflags --cflags --glibs) -I$(ROOTSYS)/include

default: Alignment

.PHONY: Alignment

all: Alignment

Alignment: Alignment.cpp
	$(CXX) -o$@ $< $(CFLAGS)

Viewer: Viewer.cpp
	$(CXX) -o$@ $< $(CFLAGS)

clean:
	rm -f Alignment
