#!/usr/bin/env make
MKDIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
BIN = beer
BIN_DBG = $(BIN).d
BUILD_DIR = build
DOXYGEN_DIR = doxygen
SRC_DIR = src
SOURCES := $(wildcard $(SRC_DIR)/*.cpp $(SRC_DIR)/*/*.cpp)
SUBDIRS := $(sort $(dir ${SOURCES}))
OBJS = $(SOURCES:%=$(BUILD_DIR)/%.o)
OBJS_DBG = $(SOURCES:%=$(BUILD_DIR)/%.do)
DEPS = $(OBJS:%=%.d)
DEPS_DBG = $(OBJS_DBG:%=%.d)
CC = g++
LD = g++
LIB_DIR = lib
# SHLIB_EXT=dylib
SHLIB_EXT=so
INCLUDE_DIRS = $(SRC_DIR) $(LIB_DIR) $(LIB_DIR)/eigen_3.3.7 $(LIB_DIR)/libtp $(LIB_DIR)/cxxopts $(LIB_DIR)/rapidjson $(LIB_DIR)/z3-4.8.9/include
CFLAGS_OPT = -g -Wfatal-errors -Werror -Wall -O3 $(INCLUDE_DIRS:%=-I%) -std=c++11 -pthread
CFLAGS_DBG = -ggdb -pg -no-pie -Wfatal-errors -Werror -Wall -O0 $(INCLUDE_DIRS:%=-I%) -std=c++11 -pthread
LDFLAGS = -pthread $(LIB_DIR)/z3-4.8.9/lib/libz3.$(SHLIB_EXT) -Wl,-rpath $(LIB_DIR)/z3-4.8.9/lib
LDFLAGS_DBG = -pg -no-pie -pthread $(LIB_DIR)/z3-4.8.9/lib/libz3.$(SHLIB_EXT) -Wl,-rpath $(LIB_DIR)/z3-4.8.9/lib

# $(info SOURCES : $(SOURCES))
# $(info OBJS    : $(OBJS))
# $(info DEPS    : $(DEPS))
# $(info SOURCES : $(SOURCES))

.SUFFIXES:

.PHONY: default jall all release debug clean doc
default: release

release: $(BIN)

debug: $(BIN_DBG)

all: release debug

doc:
	doxygen

jall: 
	$(MAKE) -j 8 all

-include $(DEPS) 
-include $(DEPS_DBG)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
	mkdir -p $(SUBDIRS:%=$(BUILD_DIR)/%)

$(BUILD_DIR)/%.cpp.do: %.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS_DBG) -c -MMD -MF $@.d -o $@ $<

$(BUILD_DIR)/%.cpp.o: %.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS_OPT) -c -MMD -MF $@.d -o $@ $<

$(BIN): $(OBJS)
	$(LD) $^ $(LDFLAGS) -o $@

$(BIN_DBG): $(OBJS_DBG)
	$(LD) $(LDFLAGS_DBG) $^ -o $@

clean:
	rm -f $(BIN)
	rm -f $(BIN_DBG)
	rm -rf $(BUILD_DIR)
	rm -rf $(DOXYGEN_DIR)
