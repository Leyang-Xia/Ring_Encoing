#include "xd/experiments.hpp"
#include <iostream>
#include <sstream>

using namespace xd;

namespace {
const char *help = R"(xordepth: XOR count / critical-path depth frontier (C++17)

  xordepth frontier MATRIX.txt [options] --output FILE
  xordepth bounded MATRIX.txt --depth D [options] --output FILE
  xordepth audit FILE.json... [--output FILE]
  xordepth verify FILE.json

Options:
  --share LIST       share thresholds to sweep, e.g. 2,3,4,0 (0 = no sharing)
  --depth-aware      break equal-frequency ties towards the earliest-ready pair
  --also FILE.json   offer stored circuits that compute the same transform
  --depth D          hard maximum output depth for bounded search
  --beam N           bounded Beam width
  --branch N         bounded candidates expanded per state
  --representations N  bounded representations retained per target
  --restarts N       bounded deterministic restarts
  --seed N           bounded candidate-order seed
  --timeout SEC      bounded search budget
  --output FILE      write the complete result

MATRIX.txt holds one output per line and one input per column, 0/1 and whitespace only.
cost counts live XOR gates; depth counts two-input XOR gates on the critical path.
Frontier results keep every candidate circuit and the target masks, so 'verify'
re-checks them symbolically. Heuristic upper bounds, not optimality claims.
)";

int parse_int(const std::string &text) {
    std::size_t count = 0;
    auto n = std::stoll(text, &count);
    ringlab::require(count == text.size() && n >= 0 && n <= 1000000,
                     "expected integer in [0,1000000]: " + text);
    return int(n);
}

// The full result carries every circuit; stdout only shows the measurements.
Json condensed(Json result) {
    if (result.has("target_masks"))
        result.obj.erase("target_masks");
    if (result.has("circuit"))
        result.obj.erase("circuit");
    if (result.has("candidates"))
        for (auto &item : result["candidates"].arr) {
            if (item.has("circuit"))
                item.obj.erase("circuit");
            if (item.has("output_depth"))
                item.obj.erase("output_depth");
        }
    if (result.has("files"))
        for (auto &file : result["files"].arr) {
            if (file.has("entries"))
                file.obj.erase("entries");
        }
    return result;
}
} // namespace

int main(int argc, char **argv) {
    try {
        if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
            std::cout << help;
            return 0;
        }
        std::string command = argv[1];
        Options o;
        BoundedScheduleOptions bounded_options;
        bool depth_set = false;
        std::vector<std::string> positional;
        for (int i = 2; i < argc; ++i) {
            std::string key = argv[i];
            auto value = [&]() {
                ringlab::require(i + 1 < argc, "missing value for " + key);
                return std::string(argv[++i]);
            };
            if (key == "--help") {
                std::cout << help;
                return 0;
            } else if (key == "--share") {
                std::istringstream list(value());
                std::string token;
                o.shares.clear();
                while (std::getline(list, token, ','))
                    o.shares.push_back(parse_int(token));
                ringlab::require(!o.shares.empty(), "--share requires at least one threshold");
            } else if (key == "--depth-aware") {
                o.depth_aware = true;
            } else if (key == "--also") {
                o.also.push_back(value());
            } else if (key == "--depth") {
                bounded_options.max_depth = parse_int(value());
                depth_set = true;
            } else if (key == "--beam") {
                bounded_options.beam_width = parse_int(value());
            } else if (key == "--branch") {
                bounded_options.branch_width = parse_int(value());
            } else if (key == "--representations") {
                bounded_options.max_representations = parse_int(value());
            } else if (key == "--restarts") {
                bounded_options.restarts = parse_int(value());
            } else if (key == "--seed") {
                bounded_options.seed = std::stoull(value());
            } else if (key == "--timeout") {
                auto text = value();
                std::size_t count = 0;
                bounded_options.timeout_seconds = std::stod(text, &count);
                ringlab::require(count == text.size(), "invalid timeout");
            } else if (key == "--output") {
                o.output = value();
            } else if (key.rfind("--", 0) == 0) {
                ringlab::require(false, "unknown option: " + key);
            } else {
                positional.push_back(key);
            }
        }
        Json result;
        if (command == "frontier") {
            ringlab::require(positional.size() == 1, "frontier takes exactly one matrix");
            ringlab::require(!o.output.empty(), "frontier requires --output");
            int inputs = 0;
            auto rows = read_matrix(positional[0], inputs);
            result = frontier(rows, inputs,
                              Json::object({{"matrix", positional[0]}}), o);
        } else if (command == "bounded") {
            ringlab::require(positional.size() == 1, "bounded takes exactly one matrix");
            ringlab::require(depth_set, "bounded requires --depth");
            ringlab::require(!o.output.empty(), "bounded requires --output");
            int inputs = 0;
            auto rows = read_matrix(positional[0], inputs);
            result = bounded_result(rows, inputs,
                                    Json::object({{"matrix", positional[0]}}), bounded_options);
        } else if (command == "audit") {
            ringlab::require(!positional.empty(), "audit needs at least one JSON file");
            result = audit(positional);
        } else if (command == "verify") {
            ringlab::require(positional.size() == 1, "verify takes exactly one result file");
            auto input = Json::read(positional[0]);
            if (input.has("method") && input["method"].string() == "DEPTH_BOUNDED_XSET_BEAM")
                result = verify_bounded_result(positional[0]);
            else
                result = verify_result(positional[0]);
        } else {
            throw std::runtime_error("unknown command: " + command);
        }
        if (!o.output.empty())
            result.write(o.output);
        std::cout << condensed(result).dump() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "xordepth: " << e.what() << '\n';
        return 1;
    }
}
