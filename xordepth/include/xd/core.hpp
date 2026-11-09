#pragma once
#include "xd/json.hpp"
#include <array>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>
namespace xd {
struct BitMask {
    // Exact multi-word bit set. Up to kInlineWords (128 bits) live inside the
    // object so the synthesis hot paths (pair xor, hash lookups, state copies)
    // never touch the allocator; wider masks spill to the heap with identical
    // value semantics. Still an exact multi-word linear-form mask: never a
    // single 64-bit integer.
    static constexpr std::uint32_t kInlineWords = 2;
    std::uint64_t inline_[kInlineWords] = {0, 0};
    std::uint64_t *words_ = inline_;
    std::uint32_t size_ = 0;
    std::uint32_t capacity_ = kInlineWords;

    BitMask() = default;
    explicit BitMask(std::uint64_t value);
    BitMask(const BitMask &other);
    BitMask(BitMask &&other) noexcept;
    BitMask &operator=(const BitMask &other);
    BitMask &operator=(BitMask &&other) noexcept;
    ~BitMask();

    static BitMask unit(std::size_t bit);
    static BitMask decimal(const std::string &value);
    std::string decimal() const;
    void toggle(std::size_t bit);
    bool bit(std::size_t bit) const;
    std::size_t bits() const;
    int count() const;
    bool empty() const { return size_ == 0; }
    BitMask &operator^=(const BitMask &other);
    friend BitMask operator^(BitMask a, const BitMask &b) { return a ^= b; }
    bool operator==(const BitMask &b) const;
    bool operator!=(const BitMask &b) const { return !(*this == b); }
    bool operator<(const BitMask &b) const;
    void trim();

    // Implementation accessors; the word layout is an internal detail.
    const std::uint64_t *word_data() const { return words_; }
    std::uint32_t word_size() const { return size_; }

  private:
    void reserve_words(std::uint32_t words);
};
struct BitMaskHash { std::size_t operator()(const BitMask &mask) const; };
using BitRows = std::vector<BitMask>;
using BitValues = std::unordered_map<int, BitMask>;
struct Circuit {
    int inputs = 0;
    bool has_zero = true;
    int zero_node = -1;
    std::vector<std::array<int, 3>> ops;
    std::vector<int> outputs;
    int cost() const { return static_cast<int>(ops.size()); }
    static Circuit from_json(const Json &value);
    Json json() const;
    Circuit pruned() const;
    BitValues verify(const BitRows &targets) const;
    std::vector<int> evaluate(const BitMask &data) const;
};
class SeedRandom {
    std::array<std::uint32_t, 624> state_{};
    int index_ = 624;
    std::uint32_t next();
  public:
    explicit SeedRandom(std::uint64_t seed);
    std::size_t below(std::size_t n);
    void shuffle(std::vector<int> &values);
};
class GreedyNetwork {
    std::unordered_map<BitMask, int, BitMaskHash> by_mask_;
  public:
    int inputs;
    std::vector<BitMask> masks;
    std::vector<std::array<int, 3>> ops;
    explicit GreedyNetwork(int input_count);
    int xor_nodes(int a, int b);
    std::vector<int> synthesize(std::vector<std::vector<int>> expressions, std::uint64_t seed = 0);
    std::vector<int> synthesize_masks(const BitRows &targets, std::uint64_t seed = 0);
    Circuit export_circuit(const std::vector<int> &outputs) const;
};
Circuit synthesize(const BitRows &targets, int inputs, std::uint64_t seed = 0);
} // namespace xd
