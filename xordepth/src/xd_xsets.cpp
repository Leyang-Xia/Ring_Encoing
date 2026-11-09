#include "xd/xsets.hpp"
#include "xd/experiments.hpp"
#include <algorithm>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <unordered_map>

namespace xd {
namespace {

std::string trim(const std::string &text) {
    std::size_t begin = text.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return {};
    std::size_t end = text.find_last_not_of(" \t\r\n");
    return text.substr(begin, end - begin + 1);
}

int parse_index(const std::string &kind, const std::string &digits) {
    std::size_t count = 0;
    long long value = 0;
    try {
        value = std::stoll(digits, &count);
    } catch (const std::exception &) {
        require(false, kind + " index is not a number: " + digits);
    }
    require(count == digits.size() && value >= 0 && value < 100000000,
            kind + " index out of range: " + digits);
    return static_cast<int>(value);
}

// Exact decimal-string mask list, matching the storage used by the other
// result records so merge/audit can read the same fields.
Json masks_to_json(const BitRows &rows) {
    auto out = Json::array();
    for (const auto &row : rows)
        out.arr.emplace_back(row.decimal());
    return out;
}

int max_output_depth(const DepthProfile &profile) {
    int deepest = 0;
    for (int d : profile.output_depth) deepest = std::max(deepest, d);
    return deepest;
}

// Distribution of per-output depths; shows how many targets sit on each
// level, which a single max alone hides.
Json depth_histogram(const DepthProfile &profile) {
    std::map<int, int> counts;
    for (int d : profile.output_depth) ++counts[d];
    auto out = Json::object();
    for (const auto &[depth, count] : counts)
        out.obj.emplace(std::to_string(depth), Json(count));
    return out;
}

} // namespace

XsetsSchedule parse_xsets_schedule(const std::string &path, int inputs) {
    std::ifstream in(path);
    require(in.good(), "cannot open xsets schedule: " + path);
    return parse_xsets_schedule(in, inputs);
}

XsetsSchedule parse_xsets_schedule(std::istream &in, int inputs) {
    require(inputs > 0, "xsets schedule needs a positive input count");
    require(static_cast<bool>(in), "xsets schedule stream is not readable");

    std::string line;
    bool in_schedule = false;
    XsetsSchedule out;
    // Signal name -> circuit node id. xJ aliases live at node J.
    std::unordered_map<std::string, int> defined;
    // Binding lines appear in arbitrary file order; targets keep b order.
    std::vector<int> binding_node;

    const std::regex alias_re(R"(^s(\d+) = x(\d+)$)");
    const std::regex gate_re(R"(^s(\d+) = s(\d+) ([+^]) s(\d+)$)");
    const std::regex bind_re(R"(^b(\d+) = (s\d+|x\d+)$)");
    const std::regex total_re(R"(^Total XOR operations: (\d+)$)");
    const std::regex depth_re(R"(^Schedule Depth: (\d+)$)");

    auto node_of = [&](const std::string &name) {
        auto it = defined.find(name);
        require(it != defined.end(), "schedule uses undefined signal: " + name);
        return it->second;
    };

    while (std::getline(in, line)) {
        std::string text = trim(line);
        if (text.empty()) continue;
        if (!in_schedule) {
            in_schedule = text == "Schedule:";
            continue;
        }
        std::smatch m;
        if (std::regex_match(text, m, alias_re)) {
            const std::string name = "s" + m.str(1);
            require(!defined.count(name), "schedule redefines signal: " + name);
            const int input = parse_index("input", m.str(2));
            require(input < inputs, "schedule input index exceeds the matrix: x" + m.str(2));
            defined.emplace(name, input);
        } else if (std::regex_match(text, m, gate_re)) {
            const std::string name = "s" + m.str(1);
            require(!defined.count(name), "schedule redefines signal: " + name);
            // Gate destinations are fresh node ids appended in file order, so
            // the schedule stays topological exactly like the source run.
            const int dst = inputs + static_cast<int>(out.circuit.ops.size());
            const int a = node_of("s" + m.str(2));
            const int b = node_of("s" + m.str(4));
            out.circuit.ops.push_back({dst, a, b});
            defined.emplace(name, dst);
        } else if (std::regex_match(text, m, bind_re)) {
            const int target = parse_index("target", m.str(1));
            if (target >= static_cast<int>(binding_node.size()))
                binding_node.resize(target + 1, -1);
            require(binding_node[target] < 0, "schedule binds target twice: b" + m.str(1));
            binding_node[target] = node_of(m.str(2));
        } else if (std::regex_match(text, m, total_re)) {
            require(out.declared_xor < 0, "duplicate xsets trailer");
            out.declared_xor = parse_index("trailer", m.str(1));
        } else if (std::regex_match(text, m, depth_re)) {
            require(out.declared_depth < 0, "duplicate xsets depth trailer");
            out.declared_depth = parse_index("trailer", m.str(1));
        } else if (text.rfind("# ", 0) == 0) {
            // Structured summary comment; the machine fields it carries are
            // re-derived here, so it is only informational.
        } else {
            require(false, "unrecognized xsets schedule line: " + text);
        }
    }
    require(in_schedule, "xsets schedule has no Schedule: section");
    for (std::size_t t = 0; t < binding_node.size(); ++t)
        require(binding_node[t] >= 0,
                "xsets schedule never binds target b" + std::to_string(t));
    out.circuit.outputs = binding_node;
    require(out.declared_xor < 0 ||
                out.declared_xor == static_cast<int>(out.circuit.ops.size()),
            "xsets trailer disagrees with the parsed gate count");
    out.circuit.inputs = inputs;
    out.circuit.has_zero = false;
    return out;
}

Json rebalance_xsets(const std::string &schedule_path, const std::string &matrix_path) {
    std::ifstream in(schedule_path);
    require(in.good(), "cannot open xsets schedule: " + schedule_path);
    return rebalance_xsets(in, schedule_path, matrix_path);
}

Json rebalance_xsets(std::istream &schedule_stream, const std::string &schedule_label,
                     const std::string &matrix_path) {
    int inputs = 0;
    const BitRows rows = read_matrix(matrix_path, inputs);
    XsetsSchedule schedule = parse_xsets_schedule(schedule_stream, inputs);
    Circuit original = std::move(schedule.circuit);
    require(static_cast<int>(original.outputs.size()) == static_cast<int>(rows.size()),
            "xsets schedule binds " + std::to_string(original.outputs.size()) +
                " targets but the matrix has " + std::to_string(rows.size()));

    // Gate-for-gate check of the import: this throws unless every output
    // mask reproduces its matrix row exactly.
    original.verify(rows);
    const DepthProfile before = depth_profile(original);
    require(schedule.declared_depth < 0 ||
                schedule.declared_depth == max_output_depth(before),
            "xsets Schedule Depth trailer disagrees with the imported circuit");

    Circuit rebalanced_circuit = rebalanced(original);
    rebalanced_circuit.verify(rows);
    const DepthProfile after = depth_profile(rebalanced_circuit);

    auto stage = [](const Circuit &circuit, const DepthProfile &profile) {
        return Json::object({{"cost", circuit.cost()},
                             {"depth", profile.depth},
                             {"max_output_depth", max_output_depth(profile)},
                             {"depth_histogram", depth_histogram(profile)}});
    };
    Json result = Json::object(
        {{"method", "XSETS_REBALANCE"},
         {"schedule", schedule_label},
         {"matrix", matrix_path},
         {"inputs", inputs},
         {"targets", static_cast<int>(rows.size())},
         {"declared_xor", schedule.declared_xor},
         {"declared_depth", schedule.declared_depth},
         {"original", stage(original, before)},
         {"rebalanced", stage(rebalanced_circuit, after)},
         {"cost", rebalanced_circuit.cost()},
         {"depth", after.depth},
         {"circuit", rebalanced_circuit.json()},
         {"target_masks", masks_to_json(rows)},
         {"note", "Imported from an X-Sets schedule, verified gate-for-gate "
                  "against the matrix, then re-associated with rebalanced(); "
                  "the gate count never grows."}});
    return result;
}
} // namespace xd
