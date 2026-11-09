CXX ?= c++
INCLUDE ?= include
BUILD ?= build

CPPFLAGS += -I$(INCLUDE)
CXXFLAGS ?= -O3 -DNDEBUG
override CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic -pthread

CORE = src/xd_json.cpp src/xd_core_mask.cpp src/xd_core_circuit.cpp src/xd_core_network.cpp
LOCAL = src/xd_depth.cpp src/xd_greedy.cpp src/xd_rebalance.cpp src/xd_rewrite.cpp src/xd_frontier.cpp src/xd_search.cpp src/xd_xsets.cpp

.PHONY: all generate test test-full clean
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

$(BUILD)/matrix-suite-test: tests/matrix_suite_test.cpp $(LOCAL) $(CORE) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/xordepth-matrix-gen: apps/xordepth_matrix_gen.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

generate: $(BUILD)/xordepth-matrix-gen
	$(BUILD)/xordepth-matrix-gen data/matrices

test: all $(BUILD)/xd-test $(BUILD)/xd-search-test $(BUILD)/trace-replay-test $(BUILD)/matrix-suite-test generate
	$(BUILD)/xd-test
	$(BUILD)/xd-search-test
	$(BUILD)/trace-replay-test
	$(BUILD)/matrix-suite-test data/matrices/sample.manifest
	bash tests/cli_test.sh $(BUILD)/xordepth
	bash tests/independence_test.sh

test-full: test
	$(BUILD)/matrix-suite-test data/matrices/full.manifest

clean:
	rm -rf $(BUILD)
