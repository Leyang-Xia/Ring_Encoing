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
  ringlab analyze MATRIX.txt [--method CSE|OUR_DCCSE|OUR_DCCSE_V2|MW_SQ|UBER_XSET] [--seeds N]

  ringlab generate MATRIX.txt --rows M --columns N --pattern sparse|dense|blocks
                   --seed N --output manifest.json

Optimizer options: --method MW_SQ|UBER_XSET --params THRESH NSTART GIVEUP
                   --timeout SECONDS --backend PATH
Depth-CSE options (analyze --method OUR_DCCSE):
                   --depth D --seeds N --seed N --lookahead --no-refine
V2 options:        --incumbent FILE (verified circuit or same-target frontier)
                   V2 runs both greedy variants: 2*N total restarts.
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
        DepthCseOptions depth_cse;
        std::string incumbent;
        int bench_rows=0,bench_columns=0;std::string bench_pattern="dense";
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
            else if (key == "--rows")
                bench_rows=parse_int(value());
            else if (key == "--columns")
                bench_columns=parse_int(value());
            else if (key == "--pattern")
                bench_pattern=value();
            else if (key == "--depth")
                depth_cse.max_depth = parse_int(value());
            else if (key == "--seed")
                depth_cse.seed = std::uint64_t(parse_int(value()));
            else if (key == "--incumbent")
                incumbent=value();
            else if (key == "--lookahead")
                depth_cse.lookahead = true;
            else if (key == "--no-refine")
                depth_cse.refine = false;
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
                    (command == "analyze" && (o.method == "CSE" || o.method == "OUR_DCCSE" || o.method == "OUR_DCCSE_V2")),
                "invalid optimizer method");
        require(incumbent.empty() || (command=="analyze" && o.method=="OUR_DCCSE_V2"),"--incumbent requires analyze --method OUR_DCCSE_V2");
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
        } else if (command == "generate") {
            require(!name.empty(),"generate requires a matrix output path");
            auto rows=benchmark_inputs(bench_rows,bench_columns,bench_pattern,depth_cse.seed);
            std::ofstream file(name);require(bool(file),"cannot write benchmark matrix");
            int ones=0;for(auto &row:rows){ones+=row.count();for(int j=0;j<bench_columns;++j)file<<(row.bit(j)?'1':'0');file<<'\n';}
            require(bool(file),"benchmark matrix write failed");
            result=Json::object({{"status","generated"},{"matrix",name},{"outputs",bench_rows},{"inputs",bench_columns},{"pattern",bench_pattern},{"seed",depth_cse.seed},{"ones",ones},{"depth_lower_bound",linear_depth_lower_bound(rows)},{"target_masks",to_json(rows)}});
        } else if (command == "analyze") {
            require(!name.empty(), "analyze requires a binary matrix file");
            int inputs;
            auto rows = read_matrix(name, inputs);
            if (o.method == "OUR_DCCSE" || o.method == "OUR_DCCSE_V2") {
                depth_cse.restarts=seeds;
                auto c=synthesize_depth_cse(rows,inputs,depth_cse);
                int limit=depth_cse.max_depth<0?linear_depth_lower_bound(rows):depth_cse.max_depth;
                int original_cost=c.cost(), standalone_cost=c.cost();
                std::string initial_source="own synthesis";
                if(o.method=="OUR_DCCSE_V2") {
                    auto alternate_options=depth_cse;
                    alternate_options.lookahead=!depth_cse.lookahead;
                    auto alternate=synthesize_depth_cse(rows,inputs,alternate_options);
                    if(alternate.cost()<c.cost() || (alternate.cost()==c.cost() && circuit_depth(alternate)<circuit_depth(c))) c=alternate;
                    c=optimize_depth_gates(c,rows,limit); standalone_cost=c.cost();
                    if(!incumbent.empty()) {
                        auto evidence=Json::read(incumbent);
                        std::vector<Json> entries;
                        if(evidence.has("candidates")) entries=evidence["candidates"].arr;
                        else entries.push_back(evidence);
                        for(auto &entry:entries) {
                            require(entry.has("circuit"),"incumbent requires circuit or frontier candidates");
                            auto candidate=Circuit::from_json(entry["circuit"]).pruned();
                            require(candidate.inputs==inputs,"incumbent input mismatch");
                            candidate.verify(rows);
                            if(circuit_depth(candidate)<=limit && candidate.cost()<c.cost()) {
                                c=candidate;initial_source=incumbent;
                                if(entry.has("label")) initial_source+=":"+entry["label"].string();
                            }
                        }
                        c=optimize_depth_gates(c,rows,limit);
                    }
                }
                result=Json::object({{"original_ours_cost",original_cost},{"standalone_cost",standalone_cost},{"initial_source",initial_source},{"status","verified"},{"method",o.method},{"seeds",seeds},{"seed",depth_cse.seed},{"lookahead",depth_cse.lookahead},{"refine",depth_cse.refine},{"depth_limit",limit},{"depth",circuit_depth(c)},{"cost",c.cost()},{"circuit",c.json()}});
                if(o.method=="OUR_DCCSE_V2") result["total_synthesis_restarts"]=2*seeds;
                else {
                    result.obj.erase("original_ours_cost");
                    result.obj.erase("standalone_cost");
                    result.obj.erase("initial_source");
                }
            } else if (o.method == "CSE") {
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
