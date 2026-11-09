#include "xd/experiments.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;
using namespace xd;

namespace {

struct MatrixCase {
    std::string name;
    fs::path path;
};

std::vector<MatrixCase> read_manifest(const fs::path &manifest) {
    std::ifstream input(manifest);
    require(bool(input), "cannot read matrix manifest: " + manifest.string());
    std::vector<MatrixCase> cases;
    std::set<std::string> seen;
    std::string line;
    while (std::getline(input, line)) {
        const auto comment = line.find('#');
        if (comment != std::string::npos)
            line.erase(comment);
        std::istringstream tokens(line);
        std::string name;
        if (!(tokens >> name))
            continue;
        require(!(tokens >> line), "matrix manifest entry contains extra tokens");
        const fs::path relative(name);
        require(!relative.is_absolute(), "matrix manifest entries must be relative");
        require(seen.insert(name).second, "duplicate matrix manifest entry: " + name);
        cases.push_back({name, manifest.parent_path() / relative});
    }
    return cases;
}

std::set<std::string> names(const std::vector<MatrixCase> &cases, bool generated) {
    std::set<std::string> result;
    for (const auto &item : cases) {
        if ((item.path.parent_path().filename() == "generated") == generated)
            result.insert(item.name);
    }
    return result;
}

void require_names(const std::vector<MatrixCase> &cases, bool generated,
                   const std::set<std::string> &expected, const std::string &label) {
    const auto actual = names(cases, generated);
    require(actual == expected, label + " matrix manifest entries do not match the required set");
}

using GeneratedKey = std::tuple<std::string, int, int, int>;

GeneratedKey generated_key(const MatrixCase &item) {
    static const std::regex pattern(
        R"(^(dense|sparse|blocks|modular)_([0-9]+)x([0-9]+)_seed([0-9]+)\.txt$)");
    std::smatch match;
    const auto filename = item.path.filename().string();
    require(std::regex_match(filename, match, pattern),
            "invalid generated matrix filename: " + item.name);
    return {match[1], std::stoi(match[2]), std::stoi(match[3]), std::stoi(match[4])};
}

void check_suite_shape(const std::vector<MatrixCase> &cases, const std::string &suite) {
    if (suite == "sample") {
        require(cases.size() == 12, "sample manifest must contain 12 matrices");
        require_names(cases, false,
                      {"cauchy_7_4_direct.txt", "systematic_rs_7_4_direct.txt",
                       "vandermonde_8to17_direct.txt", "gf8_10x10_expanded.txt"},
                      "sample representative");
        require_names(cases, true,
                      {"generated/dense_10x10_seed123.txt",
                       "generated/dense_40x64_seed123.txt",
                       "generated/sparse_12x20_seed123.txt",
                       "generated/sparse_100x100_seed123.txt",
                       "generated/blocks_24x40_seed123.txt",
                       "generated/blocks_100x100_seed456.txt",
                       "generated/modular_10x10_seed456.txt",
                       "generated/modular_64x96_seed456.txt"},
                      "sample generated");
    } else if (suite == "full") {
        require(cases.size() == 57, "full manifest must contain 57 matrices");
        require_names(cases, false,
                      {"tiny_hand_2x5.txt", "dense_8x16_seed123.txt",
                       "dense_16x32_seed123.txt", "blocks_32x64_seed123.txt",
                       "cauchy_7_4_direct.txt", "systematic_rs_7_4_direct.txt",
                       "vandermonde_8to17_direct.txt",
                       "vandermonde_8to17_direct_canonical.txt", "gf8_10x10_expanded.txt"},
                      "full existing");
        std::set<GeneratedKey> expected;
        const std::vector<std::pair<int, int>> sizes{
            {10, 10}, {12, 20}, {24, 40}, {40, 64}, {64, 96}, {100, 100}};
        for (const std::string kind : {"dense", "sparse", "blocks", "modular"})
            for (const auto &[outputs, inputs] : sizes)
                for (int seed : {123, 456})
                    expected.emplace(kind, outputs, inputs, seed);
        std::set<GeneratedKey> actual;
        for (const auto &item : cases)
            if (item.path.parent_path().filename() == "generated")
                actual.insert(generated_key(item));
        require(actual == expected, "full generated matrix coverage is incomplete");
    }
}

} // namespace

int main(int argc, char **argv) {
    try {
        const fs::path manifest = argc > 1 ? argv[1] : "data/matrices/sample.manifest";
        const auto cases = read_manifest(manifest);
        check_suite_shape(cases, manifest.stem().string());
        ScheduleOptions options;
        options.share_threshold = 0;
        for (const auto &item : cases) {
            int inputs = 0;
            auto rows = read_matrix(item.path.string(), inputs);
            require(!rows.empty() && inputs > 0, "matrix suite contains an empty case");
            if (item.path.parent_path().filename() == "generated") {
                const auto &[kind, outputs, expected_inputs, seed] = generated_key(item);
                (void)kind;
                (void)seed;
                require(outputs >= 10 && outputs <= 100 && expected_inputs >= 10 &&
                            expected_inputs <= 100,
                        "generated matrix dimensions are outside 10x10 through 100x100");
                require(static_cast<int>(rows.size()) == outputs && inputs == expected_inputs,
                        "generated matrix dimensions do not match its filename: " + item.name);
            }
            auto lower = depth_lower_bound(rows);
            auto circuit = synthesize_schedule(rows, inputs, options);
            circuit.verify(rows);
            auto profile = depth_profile(circuit);
            require(profile.depth == lower, "unshared schedule missed row-weight depth bound");
            auto rewritten = rewrite_depth_gates(circuit, rows, profile.depth, 64);
            rewritten.verify(rows);
            require(rewritten.cost() <= circuit.cost(), "matrix suite rewrite increased cost");
            require(depth_profile(rewritten).depth <= profile.depth,
                    "matrix suite rewrite increased depth");
            auto completable = synthesize_completable(rows, inputs, lower, 2, 1);
            completable.verify(rows);
            require(depth_profile(completable).depth <= lower,
                    "matrix suite completable greedy exceeded the row-weight depth bound: " +
                        item.name);
        }
        std::cout << "matrix suite (" << manifest.stem().string() << "): " << cases.size()
                  << " matrices passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
