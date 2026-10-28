CXX ?= c++
BUILD ?= build
CPPFLAGS += -Iinclude -DRINGLAB_ROOT='"$(CURDIR)"'
XSETS_HASH := $(shell if command -v sha256sum >/dev/null 2>&1; then sha256sum third_party/xsets/X-Sets.cpp; else shasum -a 256 third_party/xsets/X-Sets.cpp; fi | cut -d" " -f1)
HASH_VALID := $(shell printf "%s" "$(XSETS_HASH)" | LC_ALL=C grep -E "^[0-9a-f]{64}$$")
ifeq ($(HASH_VALID),)
$(error A SHA-256 tool is required: sha256sum or shasum)
endif
CPPFLAGS += -DRINGLAB_XSETS_SHA256='"$(XSETS_HASH)"'
CXXFLAGS ?= -O3 -DNDEBUG
override CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic
CORE = $(BUILD)/json.o $(BUILD)/bitvector.o $(BUILD)/field.o $(BUILD)/circuit.o $(BUILD)/cse.o
LAB = $(CORE) $(BUILD)/matrices.o $(BUILD)/verification.o $(BUILD)/search.o $(BUILD)/backend.o
.PHONY: all test clean
all: $(BUILD)/ringlab $(BUILD)/xsets $(BUILD)/uber
$(BUILD):
	mkdir -p $@
$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@
$(BUILD)/backend.o: third_party/xsets/X-Sets.cpp
$(BUILD)/experiment-test: tests/experiment_test.cpp $(LAB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@
$(BUILD)/core-test: tests/core_test.cpp $(CORE)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@
$(BUILD)/ringlab: apps/ringlab.cpp $(LAB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@
$(BUILD)/xsets: third_party/xsets/X-Sets.cpp | $(BUILD)
	$(CXX) -std=c++17 -O3 $< -o $@
$(BUILD)/uber: third_party/uber/Uber.cpp | $(BUILD)
	$(CXX) -std=c++17 -O3 $< -o $@
test: all $(BUILD)/core-test $(BUILD)/experiment-test
	$(BUILD)/core-test data/certificates
	$(BUILD)/experiment-test data/certificates
	$(BUILD)/ringlab verify --exhaustive
	bash tests/cli_test.sh $(BUILD)/ringlab
	bash tests/xsets_quiet.sh $(BUILD)/xsets
clean:
	rm -rf $(BUILD)
-include $(wildcard $(BUILD)/*.d)
