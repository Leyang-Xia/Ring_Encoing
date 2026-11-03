CXX ?= c++
INCLUDE ?= include
BUILD ?= build

CPPFLAGS += -I$(INCLUDE)
CXXFLAGS ?= -O3 -DNDEBUG
override CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

CORE = src/xd_json.cpp src/xd_core_mask.cpp src/xd_core_circuit.cpp src/xd_core_network.cpp
LOCAL = src/xd_depth.cpp src/xd_frontier.cpp src/xd_search.cpp

.PHONY: all test clean
all: $(BUILD)/xordepth

$(BUILD):
	mkdir -p $@

$(BUILD)/xordepth: apps/xordepth.cpp $(LOCAL) $(CORE) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/xd-test: tests/xd_test.cpp $(LOCAL) $(CORE) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/xd-search-test: tests/xd_search_test.cpp $(LOCAL) $(CORE) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/trace-replay-test: tests/trace_replay_test.cpp $(CORE) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

test: all $(BUILD)/xd-test $(BUILD)/xd-search-test $(BUILD)/trace-replay-test
	$(BUILD)/xd-test
	$(BUILD)/xd-search-test
	$(BUILD)/trace-replay-test
	bash tests/cli_test.sh $(BUILD)/xordepth
	bash tests/independence_test.sh

clean:
	rm -rf $(BUILD)
