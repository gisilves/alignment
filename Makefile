CXX := `root-config --cxx`
ROOTCLING=rootcling
MARCH := `root-config --arch`
LD:=$(CXX)
SRC=./src/

UNAME := $(shell uname)

CFLAGS += $(shell root-config --ldflags --cflags --glibs) -I$(ROOTSYS)/include

default: Alignment Resolution Viewer

.PHONY: Alignment

all: Alignment

Alignment: Alignment.cpp
	$(CXX) -o$@ $< $(CFLAGS)

Resolution: Resolution.cpp
	$(CXX) -o$@ $< $(CFLAGS)

Viewer: Viewer.cpp
	$(CXX) -o$@ $< $(CFLAGS)

clean:
	rm -f Alignment
	rm -f Resolution
	rm -f Viewer
