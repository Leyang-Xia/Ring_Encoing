#pragma once
#include "xd/core.hpp"
#include "xd/json.hpp"
#include <iosfwd>
#include <string>

namespace xd {

// One parsed Plank X-Sets run. The "Schedule:" trailer of the program's
// output is an explicit gate list (sK = xJ aliases, sK = sA + sB XOR gates,
// bT = sK target bindings); it imports losslessly into a Circuit.
struct XsetsSchedule {
    Circuit circuit;
    // Value of the "Total XOR operations: N" trailer; -1 when absent. Always
    // equal to circuit.cost() for an accepted parse.
    int declared_xor = -1;
    // Value of the "Schedule Depth: N" trailer (depth over the bound
    // targets); -1 when absent. Cross-checked against the imported circuit.
    int declared_depth = -1;
};

// Parse an X-Sets output file. `inputs` is the number of primary inputs
// (x0 .. x{inputs-1}); x indices outside that range are rejected. The
// bindings b0..b{targets-1} must be present exactly once each and become the
// circuit outputs in that order, regardless of the order the b lines appear
// in. Everything before the Schedule: header (progress echo, bitmask dump)
// is skipped. Throws on any unrecognized or ill-formed line,
// use-before-definition, redefinition, or trailer/count mismatch.
XsetsSchedule parse_xsets_schedule(const std::string &path, int inputs);

// Same parser over an open stream, e.g. standard input piped straight from
// the X-Sets run.
XsetsSchedule parse_xsets_schedule(std::istream &in, int inputs);

// Import an X-Sets schedule, verify it gate-for-gate against a 01 matrix,
// re-associate it with rebalanced(), and return a complete result record:
// original and rebalanced (cost, depth, max_output_depth, depth_histogram),
// the rebalanced circuit, and the target masks. Throws when the schedule
// does not reproduce the matrix or the target counts differ.
Json rebalance_xsets(const std::string &schedule_path, const std::string &matrix_path);

// Same command with the schedule read from a stream; `schedule_label` names
// the source in the result record ("-" or "<stdin>").
Json rebalance_xsets(std::istream &schedule, const std::string &schedule_label,
                     const std::string &matrix_path);
} // namespace xd
