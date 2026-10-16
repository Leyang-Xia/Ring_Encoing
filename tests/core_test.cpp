#include "ringlab/core.hpp"
#include <iostream>
using namespace ringlab;
template <class F> void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const std::exception &) {
        rejected = true;
    }
    require(rejected, "invalid input accepted");
}
int main(int argc, char **argv) {
    try {
        require(argc == 2, "certificate directory required");
        const auto big = std::string("5192296858534827628530496329220097");
        require(Json::parse(big).dump() == big, "JSON integer precision");
        require(Mask::decimal(big).decimal() == big, "mask decimal precision");
        require(Mask::decimal(big) == (Mask::unit(112) ^ Mask(1)), "112-bit mask");
        require(Json::parse("{\"s\":\"\\u4e2d\\u6587\",\"a\":[true,null,-1,1.25]}")["s"].string() ==
                    "中文",
                "unicode JSON");
        rejects([] { Json::parse("[1,]"); });
        rejects([] { Json::parse("{\"a\":1,\"a\":2}"); });
        rejects([] { Mask::decimal("-1"); });
        require(polynomial_product(313, 633) == 131073, "ring modulus product");
        require(irreducible(313, 8) && !irreducible(257, 8), "irreducibility");
        require(remainder(256, 31) == 8, "polynomial remainder");
        for (std::uint64_t a = 0; a < 32; ++a)
            for (std::uint64_t b = 0; b < 32; ++b)
                require(cyclic_multiply(a, b, 5) == multiply(a, b, 5, 33),
                        "independent cyclic multiplication");
        auto c = synthesize({Mask(), Mask(1), Mask(3), Mask(3)}, 3);
        require(c.cost() == 1, "free aliases");
        c.verify({Mask(), Mask(1), Mask(3), Mask(3)});
        auto relabeled = c;
        relabeled.zero_node = -7;
        for (auto &node : relabeled.outputs)
            if (node == c.zero_node)
                node = relabeled.zero_node;
        relabeled.verify({Mask(), Mask(1), Mask(3), Mask(3)});
        require(relabeled.evaluate(Mask(1)) == std::vector<int>({0, 1, 1, 1}),
                "negative zero label evaluation");
        auto old_node = relabeled.ops[0][0];
        relabeled.ops[0][0] = 1000007;
        for (auto &node : relabeled.outputs)
            if (node == old_node)
                node = 1000007;
        relabeled.verify({Mask(), Mask(1), Mask(3), Mask(3)});
        require(relabeled.evaluate(Mask(1)) == c.evaluate(Mask(1)),
                "sparse labels preserve the transform");
        auto invalid = c;
        invalid.ops[0][1] = 9999;
        rejects([&] { invalid.verify({Mask(), Mask(1), Mask(3), Mask(3)}); });
        invalid = c;
        invalid.ops.push_back({100, 0, 2});
        rejects([&] { invalid.verify({Mask(), Mask(1), Mask(3), Mask(3)}); });
        rejects([&] { c.verify({Mask(7)}); });
        auto encoded = c.json();
        encoded["inputs"] = Json::numeric("4294967299");
        rejects([&] { Circuit::from_json(encoded); });
        encoded = c.json();
        encoded["ops"][0][1] = Json::numeric("4294967296");
        rejects([&] { Circuit::from_json(encoded); });
        encoded = c.json();
        encoded["ops"] = Json::object();
        rejects([&] { Circuit::from_json(encoded); });
        const auto data = Json::read(std::string(argv[1]) + "/ring_lift_example.json");
        auto targets = expand(matrix_json(data["A"]), 4, 4, 31);
        std::map<int, int> histogram;
        for (int seed = 0; seed < 100; ++seed) {
            auto generated = synthesize(targets, 16, seed);
            ++histogram[generated.cost()];
            if (seed == data["direct_seed"].integer())
                require(generated.same(Circuit::from_json(data["direct"])),
                        "historical CSE circuit");
        }
        for (auto &entry : histogram)
            require(entry.second == data["direct_histogram"][std::to_string(entry.first)].integer(),
                    "historical histogram");
        auto inverse = inverse_matrix({{1, 1}, {1, 2}}, 8, 285);
        require(matrix_product({{1, 1}, {1, 2}}, inverse, 8, 285) == Matrix({{1, 0}, {0, 1}}),
                "matrix inverse");
        std::cout << "core: exact integers, invalid circuits, field arithmetic, 100 historical CSE "
                     "seeds passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
