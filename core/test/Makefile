#
#  Makefile
#  test.core.maxwell
#
#  Created by Rolan Akhmedov on 26.02.19
#  Copyright © 2019 Rolan Akhmedov. All rights reserved.
#

BUILD_DIR = ../../build/core/test
SOURCES = $(wildcard *.cpp)
EXECUTABLE = $(patsubst %.cpp, $(BUILD_DIR)/%, $(SOURCES))

all: $(EXECUTABLE)

$(BUILD_DIR)/%: %.cpp
	$(CXX) $(CXX_STD)  $(GMP_FLAGS) -I../include $< $(CXX_LIB) $(MAXWELL_LIB) $(GMP_LIBS) -o $@
	$@ && echo

clean:
	rm -f $(BUILD_DIR)/*
