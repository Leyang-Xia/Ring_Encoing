#include "ringlab/experiments.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>

namespace ringlab {
Json profile(const Json &definition, const RunOptions &o) {
    auto e = build_case(definition);
    require(o.path == "direct" || o.path == "main" || o.path == "projection", "unknown path");
    auto r = run_backend(o.binary,
                         o.path == "direct" ? e.direct
                         : o.path == "main" ? e.middle
                                            : e.projection,
                         o.path == "projection" ? e.virtuals.size() : e.w * e.k, o.method, o.params,
                         o.timeout, o.path);
    r["definition"] = definition;
    return r;
}
Json search(const Json &data, const RunOptions &o) {
    require(o.family == "vandermonde" || o.family == "all", "unknown family");
    require(o.limit > 0, "limit must be positive");
    std::vector<Json> cases;
    for (auto &[phase, entries] : data["screening"].obj)
        for (auto &c : entries.arr) {
            auto e = build_case(c);
            if (o.family == "all" || is_vandermonde(e.a, e.w, e.p))
                cases.push_back(c);
        }
    std::stable_sort(cases.begin(), cases.end(), [](const Json &a, const Json &b) {
        return a["percent"].number() > b["percent"].number();
    });
    auto results = Json::array();
    for (int i = 0; i < std::min<int>(o.limit, cases.size()); ++i) {
        auto &c = cases[i];
        auto e = build_case(c);
        auto main = synthesize(e.middle, e.w * e.k),
             projection = synthesize(e.projection, e.virtuals.size()),
             direct = synthesize(e.direct, e.w * e.k);
        for (int seed = 1; seed < 10; ++seed) {
            auto p = synthesize(e.projection, e.virtuals.size(), seed);
            if (p.cost() < projection.cost())
                projection = std::move(p);
        }
        for (int seed = 1; seed < 20; ++seed) {
            auto d = synthesize(e.direct, e.w * e.k, seed);
            if (d.cost() < direct.cost())
                direct = std::move(d);
        }
        auto measurements = Json::array();
        int baseline = direct.cost();
        for (auto path : {"projection", "direct"})
            for (auto method : {"MW_SQ", "UBER_XSET"}) {
                bool p = std::string(path) == "projection";
                auto r = run_backend(o.binary, p ? e.projection : e.direct,
                                     p ? e.virtuals.size() : e.w * e.k, method, o.params, o.timeout,
                                     path);
                if (r["status"].string() == "cancelled") {
                    Json::object(
                        {{"completed", results}, {"case", c}, {"interrupted_measurement", r}})
                        .write(o.output);
                    throw std::runtime_error(
                        "search cancelled; completed cases and interrupted measurement saved");
                }
                if (r["status"].string() == "verified") {
                    auto circuit = Circuit::from_json(r["circuit"]);
                    if (p && circuit.cost() < projection.cost())
                        projection = circuit;
                    if (!p)
                        baseline = std::min(baseline, circuit.cost());
                }
                measurements.arr.push_back(r);
            }
        auto ring = compose(main, projection, e.virtuals, e.middle, e.direct);
        auto result = Json::object(
            {{"case", c},
             {"matrix_family", is_vandermonde(e.a, e.w, e.p) ? "vandermonde" : "scaled_extension"},
             {"A", to_json(e.a)},
             {"lifted", to_json(e.lifted)},
             {"main", main.json()},
             {"projection", projection.json()},
             {"projection_input_masks", to_json(e.virtuals)},
             {"ring", ring.circuit.json()},
             {"direct_cse", direct.json()},
             {"measurements", measurements},
             {"best_completed_direct", baseline},
             {"percent", 100. * (baseline - ring.circuit.cost()) / baseline},
             {"checked_nonzero_minors", check_mds(e.a, e.w, e.p)},
             {"xsets_source_sha256", backend_source_hash()}});
        results.arr.push_back(result);
        results.write(o.output);
        std::cout << Json::object({{"case", c},
                                   {"ring", ring.circuit.cost()},
                                   {"baseline", baseline},
                                   {"results_saved", o.output}})
                         .dump(-1)
                  << '\n'
                  << std::flush;
        for (auto &r : measurements.arr)
            require(r["status"].string() != "error",
                    "backend error (partial search results saved)");
    }
    return results;
}
Json compare(const Json &small, const RunOptions &o) {
    auto e = build_case(small);
    auto result = Json::object({{"source_sha256", backend_source_hash()},
                                {"definition", small},
                                {"projection_cost", 8},
                                {"records", Json::array()}});
    std::vector<std::pair<std::string, std::array<int, 3>>> profiles{{"default", {10, 10, 100}},
                                                                     {"enum2", {10, 2, 0}},
                                                                     {"enum3", {10, 3, 0}},
                                                                     {"enum3-thresh0", {0, 3, 0}}};
    for (auto &[name, params] : profiles)
        for (auto method : {"MW_SQ", "UBER_XSET"})
            for (auto path : {"direct", "ring"}) {
                bool ring = std::string(path) == "ring";
                auto r = run_backend(o.binary, ring ? e.middle : e.direct, e.w * e.k, method,
                                     params, o.timeout, path);
                r["profile"] = name;
                if (r["status"].string() == "verified") {
                    auto c = Circuit::from_json(r["circuit"]);
                    if (ring) {
                        int node = e.w * e.k;
                        for (auto &op : c.ops)
                            node = std::max(node, op[0] + 1);
                        if (c.has_zero && c.zero_node >= node)
                            node = c.zero_node + 1;
                        auto z = c.outputs;
                        c.outputs.clear();
                        for (int i = 0; i < 2; ++i)
                            for (int bit = 0; bit < 4; ++bit) {
                                c.ops.push_back({node, z[i * 5 + bit], z[i * 5 + 4]});
                                c.outputs.push_back(node++);
                            }
                        c = c.pruned();
                    }
                    c.verify(e.direct);
                    r["full_circuit"] = c.json();
                    r["total_cost"] = c.cost();
                }
                result["records"].arr.push_back(r);
                result.write(o.output);
                std::cout << Json::object(
                                 {{"profile", name},
                                  {"method", method},
                                  {"path", path},
                                  {"status", r["status"]},
                                  {"total", r.has("total_cost") ? r["total_cost"] : Json()}})
                                 .dump(-1)
                          << '\n'
                          << std::flush;
                require(r["status"].string() != "error" && r["status"].string() != "cancelled",
                        "backend error (partial comparison saved)");
            }
    return result;
}
} // namespace ringlab
