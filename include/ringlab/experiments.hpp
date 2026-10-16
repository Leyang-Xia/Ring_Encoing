#pragma once
#include "ringlab/core.hpp"
#include <array>

namespace ringlab {
struct Experiment {
    int w = 0, n = 0, k = 0, m = 0;
    std::uint64_t p = 0, modulus = 0;
    Matrix a, lifted, post;
    Rows direct, middle, virtuals, projection;
};
Experiment build_case(const Json &definition);
bool is_vandermonde(const Matrix &a, int w, std::uint64_t p);
Json verify_all(const std::string &directory, const std::string &group, int samples,
                bool exhaustive);
Json verify_record(const Json &record, int samples, const Json *original = nullptr);
Json verify_file(const Json &data, int samples, const Json *original = nullptr);
struct ProcessResult {
    int exit_code = -1;
    int cancel_signal = 0;
    bool timeout = false;
    std::string output, error;
    double seconds = 0;
};
ProcessResult execute(const std::vector<std::string> &arguments, const std::string &input,
                      double timeout);
Json run_backend(const std::string &binary, const Rows &rows, int inputs, const std::string &method,
                 const std::array<int, 3> &params, double timeout, const std::string &path);
struct RunOptions {
    std::string binary, method = "MW_SQ", path = "direct", family = "vandermonde", output;
    std::array<int, 3> params{0, 3, 0};
    double timeout = 180;
    int limit = 10;
};
Json profile(const Json &definition, const RunOptions &options);
Json search(const Json &data, const RunOptions &options);
Json compare(const Json &small, const RunOptions &options);
std::string backend_source_hash();
} // namespace ringlab
