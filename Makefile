CC=g++
INCDIR=include
CXXFLAGS=-march=native -Ofast -std=c++20 -I$(INCDIR)
LDFLAGS=-lpthread
VPATH=src:include

all: swmain

swmain: swmain.o swcnt.o swexhaust.o swutils.o swmeta.o fibogen.o

swmain.o: swmain.cpp

swcnt.o: swcnt.cpp swcnt.hpp

swexhaust.o: swexhaust.cpp swexhaust.hpp swcnt.hpp fibogen.cpp

swutils.o: swutils.cpp swutils.hpp swexhaust.hpp swcnt.hpp

swmeta.o: swmeta.cpp swmeta.hpp swcnt.hpp

fibogen.o: fibogen.cpp fibogen.hpp
