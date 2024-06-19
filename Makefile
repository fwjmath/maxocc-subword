CC=g++
CXXFLAGS=-march=native -Ofast -std=c++20
LDFLAGS=

all: swmain

swmain: swmain.o swcnt.o swexhaust.o swutils.o swmeta.o

swmain.o: swmain.cpp

swcnt.o: swcnt.cpp swcnt.hpp

swexhaust.o: swexhaust.cpp swexhaust.hpp swcnt.hpp

swutils.o: swutils.cpp swutils.hpp swexhaust.hpp swcnt.hpp

swmeta.o: swmeta.cpp swmeta.hpp swcnt.hpp


