# TrekIndex: index, baselines and benchmark driver (C++17, header-only sources, no dependencies).
# The paper's experiments were compiled with GCC 13 (g++ -O3) and run single-threaded.
#   make            build every binary into build/
#   make test       build and run the invariant tests (C++) and the brute-force reduction check (Python)
#   make example    build the index on a small synthetic DAG and compare its answers with the baselines
CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic
SANFLAGS = -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer
HEADERS = $(wildcard src/*.hpp)
BIN = build/zskeleton build/paper-baselines-native build/causenet_dag build/markov_workload build/dump_running_example
TESTS = build/zskeleton-invariants build/coretree-invariants build/causal-invariants build/reachability-gadget

all: $(BIN)

# Benchmark driver: TrekIndex (method "ctskel"), its variants and every baseline; subcommands
# bench (d-separation), causal (criteria), screen, replay, pcgen, stats.
build/zskeleton: src/zskeleton.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) src/zskeleton.cpp -o $@

# Standalone RP-separation driver for saved workloads (scripts/replay_paper_baselines.py; Exp-1 RP column).
build/paper-baselines-native: src/paper_baselines.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) src/paper_baselines.cpp -o $@

# Dataset and workload tools.
build/causenet_dag: tools/causenet_dag.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) tools/causenet_dag.cpp -o $@
build/markov_workload: tools/markov_workload.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) tools/markov_workload.cpp -o $@
# The index of the running example of the paper (Figure 1): what TrekIndex stores and how it answers Q1.
build/dump_running_example: tools/dump_running_example.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -Isrc tools/dump_running_example.cpp -o $@

build/zskeleton-invariants: tests/zskeleton_invariants.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) tests/zskeleton_invariants.cpp -o $@
build/coretree-invariants: tests/coretree_invariants.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) tests/coretree_invariants.cpp -o $@
build/causal-invariants: tests/causal_invariants.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) tests/causal_invariants.cpp -o $@
build/reachability-gadget: tests/reachability_gadget.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) tests/reachability_gadget.cpp -o $@

test: $(TESTS)
	./build/zskeleton-invariants
	./build/coretree-invariants
	./build/causal-invariants
	./build/reachability-gadget
	python3 tests/verify_zskeleton_reduction.py

# The invariant tests under AddressSanitizer and UndefinedBehaviorSanitizer.
test-sanitize:
	@mkdir -p build
	for t in zskeleton_invariants coretree_invariants causal_invariants reachability_gadget; do \
	  $(CXX) $(SANFLAGS) tests/$$t.cpp -o build/$$t-sanitize && ./build/$$t-sanitize || exit 1; done

example: build/zskeleton build/paper-baselines-native
	bash scripts/example.sh

clean:
	rm -rf build

.PHONY: all test test-sanitize example clean
