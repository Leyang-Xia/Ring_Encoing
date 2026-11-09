#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

using Matrix = std::vector<std::string>;

void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}

void ensure_nonzero(std::string &row, std::mt19937_64 &random) {
    if (row.find('1') == std::string::npos)
        row[static_cast<std::size_t>(random() % row.size())] = '1';
}

Matrix dense_matrix(int outputs, int inputs, std::uint64_t seed) {
    std::mt19937_64 random(seed);
    Matrix matrix(outputs, std::string(static_cast<std::size_t>(inputs), '0'));
    for (auto &row : matrix) {
        for (char &bit : row)
            bit = (random() & 1U) != 0 ? '1' : '0';
        ensure_nonzero(row, random);
    }
    return matrix;
}

Matrix sparse_matrix(int outputs, int inputs, std::uint64_t seed) {
    std::mt19937_64 random(seed);
    Matrix matrix(outputs, std::string(static_cast<std::size_t>(inputs), '0'));
    std::vector<int> columns(static_cast<std::size_t>(inputs));
    std::iota(columns.begin(), columns.end(), 0);
    const int weight = std::min(3, inputs);
    for (auto &row : matrix) {
        std::iota(columns.begin(), columns.end(), 0);
        for (int position = 0; position < weight; ++position) {
            const int selected = position +
                                 static_cast<int>(random() % static_cast<std::uint64_t>(inputs - position));
            std::swap(columns[static_cast<std::size_t>(position)],
                      columns[static_cast<std::size_t>(selected)]);
            row[static_cast<std::size_t>(columns[static_cast<std::size_t>(position)])] = '1';
        }
    }
    return matrix;
}

Matrix blocks_matrix(int outputs, int inputs, std::uint64_t seed) {
    std::mt19937_64 random(seed);
    Matrix matrix(outputs, std::string(static_cast<std::size_t>(inputs), '0'));
    const int blocks = std::min(4, outputs);
    for (int output = 0; output < outputs; ++output) {
        const int block = std::min(blocks - 1, output * blocks / outputs);
        const int begin = block * inputs / blocks;
        const int end = (block + 1) * inputs / blocks;
        auto &row = matrix[static_cast<std::size_t>(output)];
        for (int input = begin; input < end; ++input)
            row[static_cast<std::size_t>(input)] = (random() & 1U) != 0 ? '1' : '0';
        if (row.find('1') == std::string::npos)
            row[static_cast<std::size_t>(begin + static_cast<int>(random() %
                                                                  static_cast<std::uint64_t>(end - begin)))] = '1';
    }
    return matrix;
}

Matrix modular_matrix(int outputs, int inputs, std::uint64_t seed) {
    std::mt19937_64 random(seed);
    Matrix matrix(outputs, std::string(static_cast<std::size_t>(inputs), '0'));
    for (int output = 0; output < outputs; ++output) {
        const std::uint64_t multiplier = random() % 17 + 1;
        const std::uint64_t offset = random() % 19;
        auto &row = matrix[static_cast<std::size_t>(output)];
        for (int input = 0; input < inputs; ++input) {
            const auto value = (multiplier * static_cast<std::uint64_t>(input + 1) + offset +
                                static_cast<std::uint64_t>(output) * seed) %
                               23;
            row[static_cast<std::size_t>(input)] = value < 11 ? '1' : '0';
        }
        ensure_nonzero(row, random);
    }
    return matrix;
}

Matrix make_matrix(const std::string &kind, int outputs, int inputs, std::uint64_t seed) {
    if (kind == "dense")
        return dense_matrix(outputs, inputs, seed);
    if (kind == "sparse")
        return sparse_matrix(outputs, inputs, seed);
    if (kind == "blocks")
        return blocks_matrix(outputs, inputs, seed);
    if (kind == "modular")
        return modular_matrix(outputs, inputs, seed);
    throw std::runtime_error("unknown generated matrix kind: " + kind);
}

std::string generated_name(const std::string &kind, int outputs, int inputs, int seed) {
    return "generated/" + kind + "_" + std::to_string(outputs) + "x" +
           std::to_string(inputs) + "_seed" + std::to_string(seed) + ".txt";
}

void write_matrix(const fs::path &path, const Matrix &matrix) {
    std::ofstream output(path);
    require(bool(output), "cannot write generated matrix: " + path.string());
    for (const auto &row : matrix)
        output << row << '\n';
    require(bool(output), "failed while writing generated matrix: " + path.string());
}

void write_manifest(const fs::path &path, const std::string &title,
                    const std::vector<std::string> &entries) {
    std::ofstream output(path);
    require(bool(output), "cannot write matrix manifest: " + path.string());
    output << "# " << title << "; paths are relative to this manifest.\n";
    for (const auto &entry : entries)
        output << entry << '\n';
    require(bool(output), "failed while writing matrix manifest: " + path.string());
}

} // namespace

int main(int argc, char **argv) {
    try {
        require(argc <= 2, "usage: xordepth-matrix-gen [matrix-output-directory]");
        const fs::path output_directory = argc == 2 ? argv[1] : "data/matrices";
        fs::create_directories(output_directory / "generated");
        constexpr std::array<std::pair<int, int>, 6> sizes{
            std::pair{10, 10}, std::pair{12, 20}, std::pair{24, 40},
            std::pair{40, 64}, std::pair{64, 96}, std::pair{100, 100}};
        constexpr std::array<int, 2> seeds{123, 456};
        const std::array<std::string, 4> kinds{"dense", "sparse", "blocks", "modular"};

        std::vector<std::string> generated;
        for (const auto &kind : kinds) {
            for (const auto &[outputs, inputs] : sizes) {
                for (int seed : seeds) {
                    const auto name = generated_name(kind, outputs, inputs, seed);
                    write_matrix(output_directory / name,
                                 make_matrix(kind, outputs, inputs, static_cast<std::uint64_t>(seed)));
                    generated.push_back(name);
                }
            }
        }

        const std::vector<std::string> existing{
            "tiny_hand_2x5.txt", "dense_8x16_seed123.txt", "dense_16x32_seed123.txt",
            "blocks_32x64_seed123.txt", "cauchy_7_4_direct.txt",
            "systematic_rs_7_4_direct.txt", "vandermonde_8to17_direct.txt",
            "vandermonde_8to17_direct_canonical.txt", "gf8_10x10_expanded.txt"};
        const std::vector<std::string> sample_generated{
            "generated/dense_10x10_seed123.txt", "generated/dense_40x64_seed123.txt",
            "generated/sparse_12x20_seed123.txt", "generated/sparse_100x100_seed123.txt",
            "generated/blocks_24x40_seed123.txt", "generated/blocks_100x100_seed456.txt",
            "generated/modular_10x10_seed456.txt", "generated/modular_64x96_seed456.txt"};
        const std::vector<std::string> sample_existing{
            "cauchy_7_4_direct.txt", "systematic_rs_7_4_direct.txt",
            "vandermonde_8to17_direct.txt", "gf8_10x10_expanded.txt"};

        std::vector<std::string> sample = sample_generated;
        sample.insert(sample.end(), sample_existing.begin(), sample_existing.end());
        std::vector<std::string> full = generated;
        full.insert(full.end(), existing.begin(), existing.end());
        require(generated.size() == 48 && sample.size() == 12 && full.size() == 57,
                "internal matrix manifest count mismatch");
        write_manifest(output_directory / "sample.manifest", "12-matrix correctness sample", sample);
        write_manifest(output_directory / "full.manifest", "57-matrix full correctness suite", full);
        std::cout << "Generated " << generated.size() << " matrices; sample=" << sample.size()
                  << ", full=" << full.size() << "\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
