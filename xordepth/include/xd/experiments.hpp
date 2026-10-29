#pragma once
#include "xd/depth.hpp"
#include <string>
#include <vector>

namespace xd {
using ringlab::Json;
using ringlab::Rows;

struct Options {
    // Share thresholds to sweep; 0 or 1 means "no sharing at all".
    std::vector<int> shares{2, 3, 4, 0};
    bool depth_aware = false;
    // Stored result or certificate files whose circuits are offered as extra
    // candidates when they compute the same transform as the matrix.
    std::vector<std::string> also;
    std::string output;
};

// Count/depth Pareto set for one fixed matrix. Rows are the fixed-input targets.
Json frontier(const Rows &rows, int inputs, const Json &source, const Options &options);

// Run the hard-depth bounded search and serialize one complete result record.
Json bounded_result(const Rows &rows, int inputs, const Json &source,
                    const BoundedScheduleOptions &options);

// Depth audit of every circuit stored in the given files. Metrics only: the
// files themselves are never rewritten.
Json audit(const std::vector<std::string> &files);

// Re-check a saved frontier result: symbolic execution of every candidate
// against its stored target masks, plus a check of the declared cost and depth.
Json verify_result(const std::string &path);

// Re-check one hard-depth bounded result record.
Json verify_bounded_result(const std::string &path);

// Parse a 01 matrix: one output per line, one input per column.
Rows read_matrix(const std::string &path, int &inputs);
} // namespace xd
