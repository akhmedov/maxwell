#
#  Makefile
#  Maxwell
#
#  Created by Rolan Akhmedov on 08.07.17.
#  Copyright © 2017 Rolan Akhmedov. All rights reserved.
#

MAKE = make
CXX = gcc
CXX_LIB = -lstdc++ -lm -pthread 
CXX_STD = -std=c++14 -O2

PROJECT_DIR = $(shell pwd)

GMP_LIBS = -L $(PROJECT_DIR)/gnump/lib -lgmp -lgmpxx
GMP_FLAGS = -I $(PROJECT_DIR)/gnump/include

MAXWELL_LIB = -L $(PROJECT_DIR)/build/core -lmaxwell
MAXWELL_FLAGS = -I $(PROJECT_DIR)/core/include

MYSQL_LIBS = `mysql_config --libs`
MYSQL_FLAGS = `mysql_config --cflags`

.PHONY: core module interface

core: 
	$(MAKE) CXX='$(CXX)' CXX_STD='$(CXX_STD)' CXX_LIB='$(CXX_LIB)' GMP_FLAGS='$(GMP_FLAGS)' GMP_LIBS='$(GMP_LIBS)' MYSQL_FLAGS='$(MYSQL_FLAGS)' MYSQL_LIBS='$(MYSQL_LIBS)' -C core
	$(MAKE) CXX='$(CXX)' CXX_STD='$(CXX_STD)' CXX_LIB='$(CXX_LIB)' GMP_FLAGS='$(GMP_FLAGS)' GMP_LIBS='$(GMP_LIBS)' MAXWELL_LIB='$(MAXWELL_LIB)' -C core/test

module: 
	$(MAKE) CXX='$(CXX)' CXX_LIB='$(CXX_LIB)' CXX_STD='$(CXX_STD)' -C module/*

interfce:
	echo "Not implemented!"

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
	rm -f  build/core/{*.so,*.o,test/*}
	rm -f  build/interface/{maxwell,*.o,test/*,example/*}
	rm -fr build/module/*
