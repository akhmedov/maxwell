#
#  Makefile
#  Maxwell
#
#  Created by Rolan Akhmedov on 08.07.17.
#  Copyright © 2017 Rolan Akhmedov. All rights reserved.
#

CXX = gcc
CXX_FLAGS = -Wall -Wextra -Wformat=2 -Wold-style-definition
CXX_LIB = -lstdc++ -lm -pthread 
CXX_STD = -std=c++14 -O2

PROJECT_DIR = $(shell pwd)
SOURCES = $(wildcard source/*.cpp)
#SOURCES = $(filter-out source/main.cpp, $(SOURCES))
INCLUDE = $(wildcard include/*.hpp)
OBJECTS = $(patsubst source/%.cpp, build/%.o, $(SOURCES))

GMP_LIBS = -L $(PROJECT_DIR)/gnump/lib -lgmp -lgmpxx
GMP_FLAGS = -I $(PROJECT_DIR)/gnump/include

.PHONY: gnuplot gnump clean list

maxwell: $(OBJECTS)
	@rm -f build/test.o
	$(CXX) $(CXX_STD) build/*.o $(CXX_LIB) $(GMP_LIBS) -o build/$@

test: $(OBJECTS)	
	@rm -f build/main.o
	$(CXX) $(CXX_STD) build/*.o $(CXX_LIB) $(GMP_LIBS) -o build/$@

list:
	@find . -name '*.cpp' -o -name '*.hpp' \
	-o -name '.gitignore' -o -name 'Makefile' | xargs wc -l

build/%.o: source/%.cpp $(INCLUDE) $(PROJECT_DIR)/gnump/include/*
	@mkdir -p build
	$(CXX) $(CXX_STD) $(CXX_FLAGS) $(GMP_FLAGS) -Iinclude -c $< -o $@

gnuplot:
	tar xf archive/gnuplot-5.0.6.tar
	mkdir -p gnuplot && cd gnuplot-5.0.6; \
	./configure --with-x --without-qt --prefix=$(PROJECT_DIR)/gnuplot; \
	make && make install; \
	rm -fr $(PROJECT_DIR)/gnuplot-5.0.6

gnump:
	tar xvjf archive/gmp-6.1.2.tar.bz2
	mkdir -p gnump && cd gmp-6.1.2; \
	./configure --enable-cxx --prefix=$(PROJECT_DIR)/gnump; \
	make && make check && make install; \
	rm -fr $(PROJECT_DIR)/gmp-6.1.2

clean:
	rm -f build/*.o *.gnp build/maxwell build/test
