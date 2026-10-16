#include "ringlab/experiments.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <unistd.h>

namespace fs = std::filesystem;
using namespace ringlab;
namespace {
const char *help = R"(ringlab: finite-field / quotient-ring XOR experiments (C++17)

  ringlab verify [--group all|small|original|family|practical] [--exhaustive]
                 [--samples N] [--input FILE] [--output FILE]
  ringlab profile NAME [--path direct|main|projection] --output FILE
  ringlab profile --input DEFINITION.json [optimizer options] --output FILE
  ringlab search [--family vandermonde|all] [--limit N] --output FILE
  ringlab compare --output FILE
  ringlab analyze MATRIX.txt [--method CSE|MW_SQ|UBER_XSET] [--seeds N]

Optimizer options: --method MW_SQ|UBER_XSET --params THRESH NSTART GIVEUP
                   --timeout SECONDS --backend PATH
Common options:    --data-dir PATH --output FILE

NAME: original, small, or a named family / systematic RS / Cauchy certificate.
All costs count live XOR gates. Search results are heuristic upper bounds.
Timeouts are saved as inconclusive; execution errors produce nonzero exit status.
)";
int parse_int(const std::string &text) {
    std::size_t count = 0;
    auto n = std::stoll(text, &count);
    require(count == text.size() && n >= 0 && n <= 1000000,
            "expected integer in [0,1000000]: " + text);
    return int(n);
}
Json find_definition(const std::string &dir, const std::string &name) {
    if (name == "small")
        return Json::read(dir + "/ring_lift_example.json");
    auto family = Json::read(dir + "/structured_ring_family.json");
    if (name == "original")
        return family["primary_case"];
    if (family["examples"].has(name))
        return family["examples"][name];
    auto practical = Json::read(dir + "/practical_mds_examples.json");
    require(practical["examples"].has(name), "unknown example: " + name);
    return practical["examples"][name];
}
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
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "help") {
            std::cout << help;
            return 0;
        }
        std::string command = argv[1], name, group = "all", input;
        RunOptions o;
        std::string dir = std::string(RINGLAB_ROOT) + "/data/certificates";
        int samples = 1000, seeds = 1;
        bool exhaustive = false;
        fs::path executable = fs::absolute(argv[0]);
        char resolved[4096];
        auto len = readlink("/proc/self/exe", resolved, sizeof(resolved) - 1);
        if (len > 0) {
            resolved[len] = '\0';
            executable = resolved;
        }
        o.binary = (executable.parent_path() / "xsets").string();
        for (int i = 2; i < argc; ++i) {
            std::string key = argv[i];
            auto value = [&]() {
                require(i + 1 < argc, "missing value for " + key);
                return std::string(argv[++i]);
            };
            if (key == "--help") {
                std::cout << help;
                return 0;
            } else if (key == "--group")
                group = value();
            else if (key == "--samples")
                samples = parse_int(value());
            else if (key == "--exhaustive")
                exhaustive = true;
            else if (key == "--data-dir")
                dir = value();
            else if (key == "--input")
                input = value();
            else if (key == "--output")
                o.output = value();
            else if (key == "--path")
                o.path = value();
            else if (key == "--method")
                o.method = value();
            else if (key == "--family")
                o.family = value();
            else if (key == "--limit")
                o.limit = parse_int(value());
            else if (key == "--seeds")
                seeds = parse_int(value());
            else if (key == "--timeout") {
                auto s = value();
                std::size_t count = 0;
                o.timeout = std::stod(s, &count);
                require(count == s.size(), "invalid timeout");
            } else if (key == "--backend")
                o.binary = value();
            else if (key == "--params") {
                for (int p = 0; p < 3; ++p)
                    o.params[p] = parse_int(value());
            } else {
                require(key.rfind("--", 0) != 0 && name.empty(), "unexpected argument: " + key);
                name = key;
            }
        }
        require(o.timeout > 0 && std::isfinite(o.timeout) && o.limit > 0 && seeds > 0,
                "timeout, limit and seeds must be positive");
        require(o.family == "all" || o.family == "vandermonde",
                "family must be vandermonde or all");
        require(o.path == "main" || o.path == "projection" || o.path == "direct", "invalid path");
        require(o.method == "MW_SQ" || o.method == "UBER_XSET" ||
                    (command == "analyze" && o.method == "CSE"),
                "invalid optimizer method");
        Json result;
        if (command == "verify") {
            require(name.empty(), "verify accepts --input, not a positional argument");
            if (input.empty()) {
                result = verify_all(dir, group, samples, exhaustive);
            } else {
                std::optional<Json> original;
                auto path = dir + "/structured_ring_example.json";
                if (fs::exists(path))
                    original = Json::read(path);
                result = verify_file(Json::read(input), samples, original ? &*original : nullptr);
            }
        } else if (command == "profile") {
            require((!name.empty() != !input.empty()) && !o.output.empty(),
                    "profile requires either NAME or --input, and --output");
            result = profile(input.empty() ? find_definition(dir, name) : Json::read(input), o);
        } else if (command == "search") {
            require(name.empty() && !o.output.empty(), "search requires --output");
            search(Json::read(dir + "/structured_ring_family.json"), o);
            return 0;
        } else if (command == "compare") {
            require(name.empty() && !o.output.empty(), "compare requires --output");
            compare(Json::read(dir + "/ring_lift_example.json"), o);
            return 0;
        } else if (command == "analyze") {
            require(!name.empty(), "analyze requires a binary matrix file");
            int inputs;
            auto rows = read_matrix(name, inputs);
            if (o.method == "CSE") {
                auto c = synthesize(rows, inputs);
                for (int seed = 1; seed < seeds; ++seed) {
                    auto candidate = synthesize(rows, inputs, seed);
                    if (candidate.cost() < c.cost())
                        c = candidate;
                }
                result = Json::object({{"status", "verified"},
                                       {"method", "CSE"},
                                       {"seeds", seeds},
                                       {"cost", c.cost()},
                                       {"circuit", c.json()}});
            } else
                result =
                    run_backend(o.binary, rows, inputs, o.method, o.params, o.timeout, "direct");
            result["target_masks"] = to_json(rows);
        } else
            throw std::runtime_error("unknown command: " + command);
        if (!o.output.empty())
            result.write(o.output);
        auto summary = result;
        if (summary.has("circuit"))
            summary.obj.erase("circuit");
        if (summary.has("definition"))
            summary.obj.erase("definition");
        std::cout << summary.dump() << '\n';
        if (result.has("status") && result["status"].string() == "cancelled")
            return 128 + result["cancel_signal"].int32();
        return result.has("status") && result["status"].string() == "error" ? 1 : 0;
    } catch (const std::exception &e) {
        std::cerr << "ringlab: " << e.what() << '\n';
        return 1;
    }
}
