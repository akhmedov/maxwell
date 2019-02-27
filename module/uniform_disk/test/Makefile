#
#  Makefile
#  test.uniform_disk.module.maxwell
#
#  Created by Rolan Akhmedov on 27.02.19
#  Copyright © 2019 Rolan Akhmedov. All rights reserved.
#

CXX = gcc
CXX_LIB = -lstdc++
CXX_STD = -std=c++14 -O2

BUILD_DIR = ../../build/core/test
SOURCES = $(wildcard *.cpp)
EXECUTABLE = $(patsubst %.cpp, $(BUILD_DIR)/%, $(SOURCES))

MAXWELL_LIB = -L $(BUILD_DIR)/.. -lmaxwell
GMP_LIBS = -L ../../gnump/lib -lgmp -lgmpxx

all: $(EXECUTABLE)

$(BUILD_DIR)/%: %.cpp
	$(CXX) $(CXX_STD) $(CXX_LIB) $(MAXWELL_LIB) $(GMP_LIBS) -I../include $< -o $@
	$@ && echo

clean:
	rm -f $(BUILD_DIR)/*