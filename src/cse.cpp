#include "ringlab/core.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace ringlab {
SeedRandom::SeedRandom(std::uint64_t seed) {
    state_[0] = 19650218U;
    for (unsigned i = 1; i < 624; ++i)
        state_[i] = 1812433253U * (state_[i - 1] ^ (state_[i - 1] >> 30)) + i;
    std::vector<std::uint32_t> key{std::uint32_t(seed)};
    if (seed >> 32)
        key.push_back(std::uint32_t(seed >> 32));
    unsigned i = 1, j = 0;
    for (unsigned k = 0; k < std::max<std::size_t>(624, key.size()); ++k) {
        state_[i] = (state_[i] ^ ((state_[i - 1] ^ (state_[i - 1] >> 30)) * 1664525U)) + key[j] + j;
        if (++i >= 624) {
            state_[0] = state_[623];
            i = 1;
        }
        if (++j >= key.size())
            j = 0;
    }
    for (unsigned k = 0; k < 623; ++k) {
        state_[i] = (state_[i] ^ ((state_[i - 1] ^ (state_[i - 1] >> 30)) * 1566083941U)) - i;
        if (++i >= 624) {
            state_[0] = state_[623];
            i = 1;
        }
    }
    state_[0] = 0x80000000U;
}
std::uint32_t SeedRandom::next() {
    if (index_ == 624) {
        for (int i = 0; i < 624; ++i) {
            auto y = (state_[i] & 0x80000000U) | (state_[(i + 1) % 624] & 0x7fffffffU);
            state_[i] = state_[(i + 397) % 624] ^ (y >> 1) ^ ((y & 1) ? 0x9908b0dfU : 0);
        }
        index_ = 0;
    }
    auto y = state_[index_++];
    y ^= y >> 11;
    y ^= (y << 7) & 0x9d2c5680U;
    y ^= (y << 15) & 0xefc60000U;
    y ^= y >> 18;
    return y;
}
std::size_t SeedRandom::below(std::size_t n) {
    require(n > 0 && n < 0xffffffffULL, "invalid random bound");
    int bits = 64 - __builtin_clzll(n);
    std::uint32_t v;
    do {
        v = next() >> (32 - bits);
    } while (v >= n);
    return v;
}
void SeedRandom::shuffle(std::vector<int> &v) {
    for (std::size_t i = v.size(); i > 1; --i)
        std::swap(v[i - 1], v[below(i)]);
}
GreedyCircuit::GreedyCircuit(int n) : inputs(n) {
    require(n > 0 && n < 1000000, "invalid CSE input count");
    for (int i = 0; i < n; ++i) {
        masks.push_back(Mask::unit(i));
        by_mask_.emplace(masks.back(), i);
    }
    masks.emplace_back();
    by_mask_.emplace(Mask(), n);
}
int GreedyCircuit::xor_nodes(int a, int b) {
    auto mask = masks.at(a) ^ masks.at(b);
    auto it = by_mask_.find(mask);
    if (it != by_mask_.end())
        return it->second;
    int node = masks.size();
    masks.push_back(std::move(mask));
    by_mask_.emplace(masks.back(), node);
    ops.push_back({node, a, b});
    return node;
}
std::vector<int> GreedyCircuit::synthesize(std::vector<std::vector<int>> expressions,
                                           std::uint64_t seed) {
    for (auto &row : expressions) {
        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end()), row.end());
        row.erase(std::remove(row.begin(), row.end(), inputs), row.end());
    }
    SeedRandom random(seed);
    // Retain counts between steps; only changed rows need pair updates. Rebuilding
    // the entire hash table at every gate dominated the historical implementation.
    std::unordered_map<std::uint64_t, int> freq;
    auto update_pairs = [&](const std::vector<int> &row, int delta) {
        for (std::size_t i = 0; i < row.size(); ++i)
            for (std::size_t j = i + 1; j < row.size(); ++j) {
                auto key = (std::uint64_t(row[i]) << 32) | unsigned(row[j]);
                if (delta > 0) {
                    ++freq[key];
                } else {
                    auto it = freq.find(key);
                    require(it != freq.end() && it->second > 0, "invalid CSE pair frequency");
                    if (--it->second == 0)
                        freq.erase(it);
                }
            }
    };
    for (const auto &row : expressions)
        update_pairs(row, 1);
    while (true) {
        int best = 0;
        std::vector<std::uint64_t> choices;
        for (auto [pair, count] : freq) {
            if (count > best) {
                best = count;
                choices.clear();
            }
            if (count == best)
                choices.push_back(pair);
        }
        if (best < 2)
            break;
        std::sort(choices.begin(), choices.end());
        auto pair = choices[seed ? random.below(choices.size()) : 0];
        int a = pair >> 32, b = std::uint32_t(pair), node = xor_nodes(a, b);
        for (auto &row : expressions)
            if (std::binary_search(row.begin(), row.end(), a) &&
                std::binary_search(row.begin(), row.end(), b)) {
                update_pairs(row, -1);
                row.erase(std::lower_bound(row.begin(), row.end(), a));
                row.erase(std::lower_bound(row.begin(), row.end(), b));
                if (node != inputs) {
                    auto it = std::lower_bound(row.begin(), row.end(), node);
                    if (it != row.end() && *it == node)
                        row.erase(it);
                    else
                        row.insert(it, node);
                }
                update_pairs(row, 1);
            }
    }
    std::vector<int> outputs;
    for (auto row : expressions) {
        Mask mask;
        for (auto node : row)
            mask ^= masks.at(node);
        auto it = by_mask_.find(mask);
        if (it != by_mask_.end()) {
            outputs.push_back(it->second);
            continue;
        }
        if (seed)
            random.shuffle(row);
        int node = inputs;
        for (auto next : row)
            node = xor_nodes(node, next);
        outputs.push_back(node);
    }
    return outputs;
}
std::vector<int> GreedyCircuit::synthesize_masks(const Rows &targets, std::uint64_t seed) {
    std::vector<std::vector<int>> expressions;
    for (auto &mask : targets) {
        require(mask.bits() <= std::size_t(inputs), "mask exceeds CSE input count");
        std::vector<int> row;
        for (int bit = 0; bit < inputs; ++bit)
            if (mask.bit(bit))
                row.push_back(bit);
        expressions.push_back(std::move(row));
    }
    return synthesize(std::move(expressions), seed);
}
Circuit GreedyCircuit::export_circuit(const std::vector<int> &out) const {
    Circuit c;
    c.inputs = inputs;
    c.zero_node = inputs;
    c.ops = ops;
    c.outputs = out;
    return c.pruned();
}
Circuit synthesize(const Rows &rows, int inputs, std::uint64_t seed) {
    GreedyCircuit c(inputs);
    auto result = c.export_circuit(c.synthesize_masks(rows, seed));
    result.verify(rows);
    return result;
}

} // namespace ringlab
