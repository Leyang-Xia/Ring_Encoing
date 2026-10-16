#include "ringlab/core.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace ringlab {
Circuit Circuit::from_json(const Json &j) {
    Circuit c;
    c.inputs = j["inputs"].int32();
    require(c.inputs > 0 && c.inputs < 1000000, "invalid circuit input count");
    c.has_zero = !j["zero_node"].null();
    if (c.has_zero)
        c.zero_node = j["zero_node"].int32();
    require(j["ops"].type == Json::Type::Array && j["outputs"].type == Json::Type::Array,
            "circuit ops and outputs must be arrays");
    for (auto &op : j["ops"].arr) {
        require(op.type == Json::Type::Array && op.arr.size() == 3, "XOR gate arity");
        c.ops.push_back({op[0].int32(), op[1].int32(), op[2].int32()});
    }
    for (auto &x : j["outputs"].arr)
        c.outputs.push_back(x.int32());
    require(c.cost() == j["cost"].integer(), "declared circuit cost mismatch");
    return c;
}
Json Circuit::json() const {
    auto j = Json::object({{"inputs", inputs},
                           {"zero_node", has_zero ? Json(zero_node) : Json()},
                           {"ops", Json::array()},
                           {"outputs", Json::array()},
                           {"cost", cost()}});
    for (auto &op : ops)
        j["ops"].arr.push_back(Json::array({op[0], op[1], op[2]}));
    for (auto x : outputs)
        j["outputs"].arr.emplace_back(x);
    return j;
}
Circuit Circuit::pruned() const {
    std::unordered_map<int, std::array<int, 2>> parents;
    for (auto &op : ops)
        require(parents.emplace(op[0], std::array<int, 2>{op[1], op[2]}).second,
                "duplicate gate node");
    std::set<int> live;
    std::vector<int> stack = outputs;
    while (!stack.empty()) {
        int node = stack.back();
        stack.pop_back();
        auto it = parents.find(node);
        if (it == parents.end() || !live.insert(node).second)
            continue;
        stack.push_back(it->second[0]);
        stack.push_back(it->second[1]);
    }
    Circuit c = *this;
    c.ops.clear();
    for (auto &op : ops)
        if (live.count(op[0]))
            c.ops.push_back(op);
    return c;
}
Values Circuit::verify(const Rows &targets) const {
    require(inputs > 0 && inputs < 1000000, "invalid circuit input count");
    Values values;
    for (int i = 0; i < inputs; ++i)
        values.emplace(i, Mask::unit(i));
    if (has_zero)
        require(values.emplace(zero_node, Mask()).second, "zero node aliases input");
    for (auto &op : ops) {
        require(!values.count(op[0]) && values.count(op[1]) && values.count(op[2]),
                "duplicate or unavailable gate node");
        values.emplace(op[0], values.at(op[1]) ^ values.at(op[2]));
    }
    require(outputs.size() == targets.size(), "circuit output count mismatch");
    for (std::size_t i = 0; i < outputs.size(); ++i)
        require(values.count(outputs[i]) && values.at(outputs[i]) == targets[i],
                "circuit output mask mismatch");
    require(pruned().cost() == cost(), "counted circuit contains dead gates");
    return values;
}
std::vector<int> Circuit::evaluate(const Mask &word) const {
    require(inputs > 0, "invalid circuit input count");
    // Certificates normally use compact node IDs. Use contiguous storage for
    // repeated bit-level execution, retaining a sparse fallback for arbitrary IDs.
    std::size_t bound = 4 * (std::size_t(inputs) + ops.size() + 1);
    int largest = inputs - 1;
    bool dense = true;
    auto inspect = [&](int node) {
        if (has_zero && node == zero_node && node < 0)
            return;
        if (node < 0 || std::size_t(node) >= bound)
            dense = false;
        else
            largest = std::max(largest, node);
    };
    if (has_zero)
        inspect(zero_node);
    for (const auto &op : ops)
        for (auto node : op)
            inspect(node);
    for (auto node : outputs)
        inspect(node);
    if (dense) {
        std::vector<int> values(std::size_t(largest) + 2, -1);
        auto index = [&](int node) { return node < 0 ? largest + 1 : node; };
        for (int i = 0; i < inputs; ++i)
            values[i] = word.bit(i);
        if (has_zero) {
            require(values[index(zero_node)] == -1, "zero node aliases input");
            values[index(zero_node)] = 0;
        }
        for (const auto &op : ops) {
            int node = index(op[0]), a = index(op[1]), b = index(op[2]);
            require(values[node] == -1 && values[a] >= 0 && values[b] >= 0,
                    "duplicate or unavailable evaluation node");
            values[node] = values[a] ^ values[b];
        }
        std::vector<int> out;
        out.reserve(outputs.size());
        for (auto node : outputs) {
            require(values[index(node)] >= 0, "unavailable evaluation output");
            out.push_back(values[index(node)]);
        }
        return out;
    }
    std::unordered_map<int, int> v;
    for (int i = 0; i < inputs; ++i)
        v[i] = word.bit(i);
    if (has_zero)
        v[zero_node] = 0;
    for (auto &op : ops)
        v[op[0]] = v.at(op[1]) ^ v.at(op[2]);
    std::vector<int> out;
    for (auto node : outputs)
        out.push_back(v.at(node));
    return out;
}
Composed compose(const Circuit &main, const Circuit &projection, const Rows &virtuals,
                 const Rows &middle, const Rows &targets) {
    require(main.outputs.size() == middle.size() && projection.inputs == int(virtuals.size()),
            "composition dimensions");
    GreedyCircuit c(main.inputs);
    std::unordered_map<int, int> map;
    for (int i = 0; i < main.inputs; ++i)
        map[i] = i;
    if (main.has_zero)
        map[main.zero_node] = c.inputs;
    for (auto &op : main.ops)
        map[op[0]] = c.xor_nodes(map.at(op[1]), map.at(op[2]));
    std::unordered_map<Mask, int, MaskHash> bymask;
    for (std::size_t i = 0; i < middle.size(); ++i)
        if (!middle[i].empty())
            bymask[middle[i]] = map.at(main.outputs[i]);
    std::set<int> mainnodes;
    for (auto &op : c.ops)
        mainnodes.insert(op[0]);
    map.clear();
    for (std::size_t i = 0; i < virtuals.size(); ++i)
        map[i] = bymask.at(virtuals[i]);
    if (projection.has_zero)
        map[projection.zero_node] = c.inputs;
    for (auto &op : projection.ops)
        map[op[0]] = c.xor_nodes(map.at(op[1]), map.at(op[2]));
    std::vector<int> out;
    for (auto node : projection.outputs)
        out.push_back(map.at(node));
    Composed result;
    result.circuit = c.export_circuit(out);
    result.circuit.verify(targets);
    for (auto &op : result.circuit.ops)
        if (mainnodes.count(op[0]))
            ++result.main_cost;
        else
            ++result.projection_cost;
    return result;
}
Circuit compose_layers(const std::vector<Circuit> &stages, const Rows &targets) {
    require(!stages.empty(), "empty composition");
    GreedyCircuit c(stages[0].inputs);
    std::vector<int> values(c.inputs);
    std::iota(values.begin(), values.end(), 0);
    for (auto &stage : stages) {
        require(stage.inputs == int(values.size()), "stage composition dimensions");
        std::unordered_map<int, int> map;
        for (int i = 0; i < stage.inputs; ++i)
            map[i] = values[i];
        if (stage.has_zero)
            map[stage.zero_node] = c.inputs;
        for (auto &op : stage.ops)
            map[op[0]] = c.xor_nodes(map.at(op[1]), map.at(op[2]));
        values.clear();
        for (auto node : stage.outputs)
            values.push_back(map.at(node));
    }
    auto result = c.export_circuit(values);
    result.verify(targets);
    return result;
}

} // namespace ringlab
