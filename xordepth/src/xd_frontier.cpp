#include "xd/experiments.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace xd {
namespace {
using ringlab::require;

struct Found {
    std::string path;
    Circuit circuit;
    bool same(const Found &other) const {
        return circuit.inputs == other.circuit.inputs && circuit.ops == other.circuit.ops &&
               circuit.outputs == other.circuit.outputs;
    }
};

void collect(const Json &value, const std::string &path, std::vector<Found> &out) {
    if (value.type == Json::Type::Object) {
        if (value.has("ops") && value.has("outputs") && value.has("cost")) {
            out.push_back({path, Circuit::from_json(value)});
            return;
        }
        for (const auto &[key, child] : value.obj)
            collect(child, path + "/" + key, out);
    } else if (value.type == Json::Type::Array) {
        for (std::size_t i = 0; i < value.arr.size(); ++i)
            collect(value.arr[i], path + "/" + std::to_string(i), out);
    }
}

// Masks are stored as exact decimal strings: a binary linear form can exceed
// any fixed-width integer, so the JSON never carries them as numbers.
Json masks_to_json(const Rows &rows) {
    auto out = Json::array();
    for (const auto &row : rows)
        out.arr.emplace_back(row.decimal());
    return out;
}

Rows json_to_masks(const Json &value) {
    require(value.type == Json::Type::Array, "mask list must be an array");
    Rows rows;
    for (const auto &item : value.arr)
        rows.push_back(Mask::decimal(item.string()));
    return rows;
}

int heaviest_row(const Rows &rows) {
    int weight = 0;
    for (const auto &row : rows)
        weight = std::max(weight, row.count());
    return weight;
}

Json candidate(const std::string &label, const Circuit &c, const std::string &method, int share,
               bool depth_aware) {
    auto profile = depth_profile(c);
    auto depths = Json::array();
    for (auto depth : profile.output_depth)
        depths.arr.emplace_back(depth);
    return Json::object({{"label", label},
                         {"method", method},
                         {"share_threshold", share < 0 ? Json() : Json(share)},
                         {"depth_aware", depth_aware},
                         {"status", "verified"},
                         {"cost", c.cost()},
                         {"depth", profile.depth},
                         {"max_fanout", profile.max_fanout},
                         {"output_depth", depths},
                         {"circuit", c.json()}});
}

void mark_frontier(Json &result) {
    auto &items = result["candidates"].arr;
    auto front = Json::array();
    for (std::size_t i = 0; i < items.size(); ++i) {
        bool verified = items[i]["status"].string() == "verified";
        bool dominated = false;
        for (std::size_t j = 0; verified && j < items.size(); ++j) {
            if (i == j || items[j]["status"].string() != "verified")
                continue;
            int ci = items[i]["cost"].int32(), di = items[i]["depth"].int32();
            int cj = items[j]["cost"].int32(), dj = items[j]["depth"].int32();
            dominated = dominated || (cj <= ci && dj <= di && (cj < ci || dj < di));
        }
        items[i]["pareto"] = verified && !dominated;
        if (verified && !dominated)
            front.arr.push_back(items[i]["label"]);
    }
    result["pareto_front"] = front;
}

std::string share_label(int share, bool depth_aware) {
    std::string label = share > 1 ? "CSE-share" + std::to_string(share) : "CSE-noshare";
    if (depth_aware)
        label += "-depth";
    return label;
}
} // namespace

Rows read_matrix(const std::string &path, int &inputs) {
    std::ifstream file(path);
    require(bool(file), "cannot read matrix: " + path);
    std::string line;
    Rows rows;
    inputs = 0;
    while (std::getline(file, line)) {
        Mask row;
        int bit = 0;
        for (char c : line) {
            if (c == ' ' || c == '\t' || c == '\r')
                continue;
            require(c == '0' || c == '1', "matrix must contain only 0/1 and whitespace");
            if (c == '1')
                row.toggle(bit);
            ++bit;
        }
        if (!bit)
            continue;
        if (!inputs)
            inputs = bit;
        require(bit == inputs, "ragged binary matrix");
        rows.push_back(row);
    }
    require(!rows.empty(), "empty binary matrix");
    return rows;
}

Json frontier(const Rows &rows, int inputs, const Json &source, const Options &o) {
    require(inputs > 0, "matrix needs at least one input column");
    auto result = Json::object({{"status", "verified"},
                                {"scope", "fixed-input-to-fixed-output XOR network"},
                                {"source", source},
                                {"inputs", inputs},
                                {"outputs", int(rows.size())},
                                {"max_row_weight", heaviest_row(rows)},
                                {"depth_lower_bound", depth_lower_bound(rows)},
                                {"target_masks", masks_to_json(rows)},
                                {"candidates", Json::array()},
                                {"note", "Costs count live XOR gates. depth counts two-input XOR "
                                         "gates from a primary input to an output. Every entry is "
                                         "a heuristic upper bound; the frontier is not an "
                                         "optimality claim."}});
    auto chained = ringlab::synthesize(rows, inputs);
    result["candidates"].arr.push_back(candidate("CSE-chained", chained, "CSE_CHAINED", 2, false));
    result["candidates"].arr.push_back(
        candidate("CSE-chained-rebalanced", rebalanced(chained), "REBALANCED", 2, false));
    for (int share : o.shares) {
        ScheduleOptions options;
        options.share_threshold = share;
        options.depth_aware = o.depth_aware;
        auto c = synthesize_schedule(rows, inputs, options);
        result["candidates"].arr.push_back(
            candidate(share_label(share, o.depth_aware), c, "CSE", share, o.depth_aware));
    }
    Json also = Json::array();
    for (const auto &path : o.also) {
        auto value = Json::read(path);
        std::vector<Found> found;
        collect(value, "", found);
        int used = 0, skipped = 0;
        std::vector<bool> drop(found.size(), false);
        for (std::size_t i = 0; i < found.size(); ++i)
            for (std::size_t j = 0; j < i; ++j)
                if (found[i].same(found[j]))
                    drop[i] = true;
        for (std::size_t i = 0; i < found.size(); ++i) {
            const auto &entry = found[i];
            if (drop[i]) {
                ++skipped;
                continue;
            }
            try {
                entry.circuit.verify(rows);
            } catch (const std::exception &) {
                ++skipped;
                continue;
            }
            ++used;
            auto label = std::filesystem::path(path).filename().string() + ":" + entry.path;
            result["candidates"].arr.push_back(
                candidate(label, entry.circuit, "STORED", -1, false));
            result["candidates"].arr.push_back(
                candidate(label + "-rebalanced", rebalanced(entry.circuit), "REBALANCED", -1,
                          false));
        }
        also.arr.push_back(Json::object({{"file", path}, {"used", used}, {"skipped", skipped}}));
    }
    if (!o.also.empty())
        result["stored_candidates"] = also;
    mark_frontier(result);
    return result;
}

Json bounded_result(const Rows &rows, int inputs, const Json &source,
                    const BoundedScheduleOptions &options) {
    auto search = synthesize_bounded(rows, inputs, options);
    auto result = Json::object({{"status", search.status},
                                {"method", "DEPTH_BOUNDED_XSET_BEAM"},
                                {"scope", "fixed-input-to-fixed-output XOR network"},
                                {"source", source},
                                {"inputs", inputs},
                                {"outputs", int(rows.size())},
                                {"target_masks", masks_to_json(rows)},
                                {"depth_limit", options.max_depth},
                                {"depth_lower_bound", search.depth_lower_bound},
                                {"beam_width", options.beam_width},
                                {"branch_width", options.branch_width},
                                {"max_representations", options.max_representations},
                                {"restarts", options.restarts},
                                {"seed", std::uint64_t(options.seed)},
                                {"timeout_seconds", options.timeout_seconds},
                                {"explored_states", std::uint64_t(search.explored_states)},
                                {"pruned_depth", std::uint64_t(search.pruned_depth)},
                                {"pruned_bound", std::uint64_t(search.pruned_bound)},
                                {"cost", Json()},
                                {"note", "Depth is a hard output constraint; verified results are "
                                         "heuristic upper bounds."}});
    if (search.status == "verified") {
        auto profile = depth_profile(search.circuit);
        require(profile.depth <= options.max_depth, "bounded result exceeds depth limit");
        search.circuit.verify(rows);
        result["cost"] = search.circuit.cost();
        result["depth"] = profile.depth;
        result["max_fanout"] = profile.max_fanout;
        result["output_depth"] = Json::array();
        for (auto value : profile.output_depth) result["output_depth"].arr.emplace_back(value);
        result["circuit"] = search.circuit.json();
    }
    return result;
}

Json audit(const std::vector<std::string> &files) {
    require(!files.empty(), "audit requires at least one file");
    auto reports = Json::array();
    int total = 0, improved = 0, removed = 0;
    for (const auto &path : files) {
        std::vector<Found> found;
        collect(Json::read(path), "", found);
        auto entries = Json::array();
        for (const auto &entry : found) {
            auto before = depth_profile(entry.circuit);
            auto fixed = rebalanced(entry.circuit);
            auto after = depth_profile(fixed);
            require(fixed.cost() <= entry.circuit.cost(), "audit cost regression");
            require(after.depth <= before.depth, "audit depth regression");
            ++total;
            if (after.depth < before.depth)
                ++improved;
            removed += before.gates - after.gates;
            entries.arr.push_back(Json::object({{"path", entry.path},
                                                {"cost", before.gates},
                                                {"depth", before.depth},
                                                {"max_fanout", before.max_fanout},
                                                {"rebalanced_cost", after.gates},
                                                {"rebalanced_depth", after.depth}}));
        }
        reports.arr.push_back(Json::object(
            {{"file", path}, {"circuits", int(found.size())}, {"entries", entries}}));
    }
    return Json::object({{"status", "verified"},
                         {"scope", "depth audit of stored circuits; no file is rewritten"},
                         {"files", reports},
                         {"circuits", total},
                         {"depth_improved", improved},
                         {"gates_removed", removed},
                         {"note", "rebalanced_depth is the exact result of rebalanced(): every "
                                  "output mask is compared and the gate count never increases."}});
}

Json verify_result(const std::string &path) {
    auto data = Json::read(path);
    require(data.has("candidates") && data.has("target_masks"),
            "not a frontier result: need candidates and target_masks");
    auto rows = json_to_masks(data["target_masks"]);
    auto reports = Json::array();
    int checked = 0;
    for (const auto &entry : data["candidates"].arr) {
        const std::string label = entry["label"].string();
        const std::string status = entry["status"].string();
        if (status != "verified") {
            require(status == "timeout" || status == "cancelled",
                    "result contains an execution error");
            reports.arr.push_back(Json::object({{"label", label}, {"status", status}}));
            continue;
        }
        auto circuit = Circuit::from_json(entry["circuit"]);
        circuit.verify(rows);
        require(circuit.cost() == entry["cost"].integer(), "declared cost mismatch");
        auto profile = depth_profile(circuit);
        require(profile.depth == entry["depth"].integer(), "declared depth mismatch");
        ++checked;
        reports.arr.push_back(Json::object({{"label", label},
                                            {"status", "verified"},
                                            {"cost", circuit.cost()},
                                            {"depth", profile.depth}}));
    }
    return Json::object({{"status", "verified"},
                         {"file", path},
                         {"candidates_verified", checked},
                         {"candidates", reports}});
}

Json verify_bounded_result(const std::string &path) {
    auto data = Json::read(path);
    require(data.has("method") && data["method"].string() == "DEPTH_BOUNDED_XSET_BEAM",
            "not a bounded result");
    require(data.has("target_masks") && data.has("depth_limit"),
            "bounded result is missing target masks or depth limit");
    auto rows = json_to_masks(data["target_masks"]);
    const std::string status = data["status"].string();
    require(status == "verified" || status == "timeout" || status == "search_failed" ||
                status == "infeasible_by_lower_bound",
            "invalid bounded result status");
    if (status != "verified") {
        require(data.has("cost") && data["cost"].null(),
                "unfinished bounded result must have null cost");
        require(!data.has("circuit"), "unfinished bounded result cannot contain a circuit");
        return Json::object({{"status", "verified"},
                             {"file", path},
                             {"result_status", status},
                             {"candidates_verified", 0}});
    }
    require(data.has("circuit") && data.has("cost") && data.has("depth"),
            "verified bounded result is missing circuit metrics");
    auto circuit = Circuit::from_json(data["circuit"]);
    circuit.verify(rows);
    require(circuit.cost() == data["cost"].integer(), "bounded cost mismatch");
    auto profile = depth_profile(circuit);
    require(profile.depth == data["depth"].integer(), "bounded depth mismatch");
    require(profile.depth <= data["depth_limit"].integer(), "bounded depth limit exceeded");
    return Json::object({{"status", "verified"},
                         {"file", path},
                         {"result_status", "verified"},
                         {"cost", circuit.cost()},
                         {"depth", profile.depth},
                         {"candidates_verified", 1}});
}
} // namespace xd
