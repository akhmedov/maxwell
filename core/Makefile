#
#  Makefile
#  core.maxwell
#
#  Created by Rolan Akhmedov on 25.02.19
#  Copyright © 2019 Rolan Akhmedov. All rights reserved.
#

CXX_FLAGS = -fPIC
CXX_FLAGS += $(if ${DEBUG},-Wall -Wextra -Wformat=2 -Wold-style-definition -DDEBUG,)
CXX_FLAGS += $(if ${PDISK_LINEAR_INT},-DNUMERIC_PDISK_LINEAR_INT,)
CXX_FLAGS += $(if ${AUW_NOICE},-DAUW_NOICE,)
CXX_FLAGS += $(if ${MGW_NOICE},-DMGW_NOICE,)

BUILD_DIR = ../build/core
SOURCES = $(wildcard source/*.cpp)
INCLUDE = $(wildcard include/*.hpp)
OBJECTS = $(patsubst source/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES))

all: $(OBJECTS)
	$(CXX) $(CXX_STD) $(BUILD_DIR)/*.o $(CXX_LIB) $(GMP_LIBS) $(MYSQL_LIBS) -shared -o $(BUILD_DIR)/libmaxwell.so

$(BUILD_DIR)/%.o: source/%.cpp $(INCLUDE)
	$(CXX) $(CXX_STD) $(CXX_FLAGS) $(MYSQL_FLAGS) $(GMP_FLAGS) -Iinclude -c $< -o $@

list:
	@find . -name '*.cpp' -o -name '*.hpp' -o -name 'Makefile' | xargs wc -l

clean:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/libmaxwell.so
