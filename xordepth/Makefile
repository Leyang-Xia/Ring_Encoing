CXX ?= c++
CORE ?= ../src
INCLUDE ?= ../include
BUILD ?= build

# This experiment reuses the ringlab core primitives (exact JSON integers, masks,
# circuit construction, Plank-style CSE) and nothing from its experiment layer.
CPPFLAGS += -I$(INCLUDE) -Iinclude
CXXFLAGS ?= -O3 -DNDEBUG
override CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

SHARED = $(CORE)/json.cpp $(CORE)/bitvector.cpp $(CORE)/circuit.cpp $(CORE)/cse.cpp
LOCAL = src/xd_depth.cpp src/xd_frontier.cpp src/xd_search.cpp

.PHONY: all test clean
all: $(BUILD)/xordepth

$(BUILD):
	mkdir -p $@

$(BUILD)/xordepth: apps/xordepth.cpp $(LOCAL) $(SHARED) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/xd-test: tests/xd_test.cpp $(LOCAL) $(SHARED) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/xd-search-test: tests/xd_search_test.cpp $(LOCAL) $(SHARED) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

test: all $(BUILD)/xd-test
	$(BUILD)/xd-test
	bash tests/cli_test.sh $(BUILD)/xordepth

clean:
	rm -rf $(BUILD)
