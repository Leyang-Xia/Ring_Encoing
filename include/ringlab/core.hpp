#pragma once
#include "ringlab/json.hpp"
#include <array>
#include <functional>
#include <unordered_map>

namespace ringlab {
struct Mask {
    std::vector<std::uint64_t> words;
    Mask() = default;
    explicit Mask(std::uint64_t value);
    static Mask unit(std::size_t bit);
    static Mask decimal(const std::string &value);
    std::string decimal() const;
    void toggle(std::size_t bit);
    bool bit(std::size_t bit) const;
    std::size_t bits() const;
    int count() const;
    bool empty() const {
        return words.empty();
    }
    Mask &operator^=(const Mask &other);
    friend Mask operator^(Mask a, const Mask &b) {
        return a ^= b;
    }
    bool operator==(const Mask &b) const {
        return words == b.words;
    }
    bool operator!=(const Mask &b) const {
        return !(*this == b);
    }
    bool operator<(const Mask &b) const;
    void trim();
};
struct MaskHash {
    std::size_t operator()(const Mask &mask) const;
};
using Rows = std::vector<Mask>;
using Matrix = std::vector<std::vector<std::uint64_t>>;
using Values = std::unordered_map<int, Mask>;
std::uint64_t polynomial_product(std::uint64_t a, std::uint64_t b);
std::uint64_t cyclic_multiply(std::uint64_t a, std::uint64_t b, int width);
std::uint64_t remainder(std::uint64_t a, std::uint64_t modulus);
std::uint64_t multiply(std::uint64_t a, std::uint64_t b, int width, std::uint64_t modulus);
std::uint64_t power(std::uint64_t a, std::uint64_t exponent, int width, std::uint64_t modulus);
bool irreducible(std::uint64_t polynomial, int width);
Rows expand(const Matrix &matrix, int input_width, int output_width, std::uint64_t modulus);
Matrix matrix_product(const Matrix &a, const Matrix &b, int w, std::uint64_t p);
Matrix inverse_matrix(const Matrix &a, int w, std::uint64_t p);
bool full_rank(Matrix a, int w, std::uint64_t p);
int check_mds(const Matrix &a, int w, std::uint64_t p);
Matrix matrix_json(const Json &value);
Json to_json(const Matrix &matrix);
Rows rows_json(const Json &value);
Json to_json(const Rows &rows);

struct Circuit {
    int inputs = 0;
    bool has_zero = true;
    int zero_node = -1;
    std::vector<std::array<int, 3>> ops;
    std::vector<int> outputs;
    int cost() const {
        return static_cast<int>(ops.size());
    }
    static Circuit from_json(const Json &value);
    Json json() const;
    Circuit pruned() const;
    Values verify(const Rows &targets) const;
    std::vector<int> evaluate(const Mask &data) const;
    bool same(const Circuit &b) const {
        return ops == b.ops && outputs == b.outputs && inputs == b.inputs;
    }
};

// MT19937's integer seeding, getrandbits and randbelow match the historical seeds.
class SeedRandom {
    std::array<std::uint32_t, 624> state_{};
    int index_ = 624;
    std::uint32_t next();

  public:
    explicit SeedRandom(std::uint64_t seed);
    std::size_t below(std::size_t n);
    void shuffle(std::vector<int> &values);
};
class GreedyCircuit {
    std::unordered_map<Mask, int, MaskHash> by_mask_;

  public:
    int inputs;
    std::vector<Mask> masks;
    std::vector<std::array<int, 3>> ops;
    explicit GreedyCircuit(int input_count);
    int xor_nodes(int a, int b);
    std::vector<int> synthesize(std::vector<std::vector<int>> expressions, std::uint64_t seed = 0);
    std::vector<int> synthesize_masks(const Rows &targets, std::uint64_t seed = 0);
    Circuit export_circuit(const std::vector<int> &outputs) const;
};
Circuit synthesize(const Rows &targets, int inputs, std::uint64_t seed = 0);
struct Composed {
    Circuit circuit;
    int main_cost = 0;
    int projection_cost = 0;
};
Composed compose(const Circuit &main, const Circuit &projection, const Rows &virtual_inputs,
                 const Rows &middle, const Rows &targets);
Circuit compose_layers(const std::vector<Circuit> &stages, const Rows &targets);
} // namespace ringlab
