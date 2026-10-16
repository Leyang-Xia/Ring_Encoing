#include "ringlab/experiments.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>

namespace ringlab {
namespace {
int integer(const Json &j, const std::string &key) {
    return j[key].int32();
}
std::vector<std::uint64_t> numbers(const Json &j) {
    std::vector<std::uint64_t> v;
    for (auto &x : j.arr) {
        auto n = x.integer();
        require(n >= 0, "negative field value");
        v.push_back(n);
    }
    return v;
}
Matrix cauchy(const std::vector<std::uint64_t> &x, const std::vector<std::uint64_t> &y, int w,
              std::uint64_t p) {
    std::set<std::uint64_t> unique(x.begin(), x.end());
    unique.insert(y.begin(), y.end());
    require(unique.size() == x.size() + y.size(), "Cauchy points must be disjoint and distinct");
    Matrix a;
    for (auto b : x) {
        std::vector<std::uint64_t> row;
        for (auto c : y) {
            require(b < (std::uint64_t(1) << w) && c < (std::uint64_t(1) << w),
                    "field point out of range");
            row.push_back(power(b ^ c, (std::uint64_t(1) << w) - 2, w, p));
        }
        a.push_back(row);
    }
    return a;
}
void dimensions(const Matrix &a, int rows, int cols) {
    require(int(a.size()) == rows, "coefficient row count");
    for (auto &r : a)
        require(int(r.size()) == cols, "coefficient column count");
}
void verify_named(const Json &j, const std::string &key, const Rows &targets) {
    Circuit::from_json(j[key]).verify(targets);
}
void check_measurement(const Json &j, const Experiment &e) {
    auto status = j["status"].string();
    if (status == "verified") {
        require(j.has("circuit"), "verified measurement missing circuit");
        const auto &path = j["path"].string();
        const auto &targets = path == "direct" ? e.direct
                              : path == "main" ? e.middle
                                               : e.projection;
        require(path == "direct" || path == "main" || path == "projection",
                "unknown measurement path");
        auto c = Circuit::from_json(j["circuit"]);
        require(c.cost() == integer(j, "cost"), "measurement cost mismatch");
        c.verify(targets);
    } else
        require((status == "timeout" || status == "cancelled") && j["cost"].null() &&
                    !j.has("circuit"),
                "invalid unfinished measurement");
}
void random_vectors(const Experiment &e, const Circuit &c, int count) {
    SeedRandom random(20261016);
    for (int sample = 0; sample < count; ++sample) {
        std::vector<std::uint64_t> data;
        Mask word;
        for (int j = 0; j < e.k; ++j) {
            auto v = random.below(std::size_t(1) << e.w);
            data.push_back(v);
            for (int bit = 0; bit < e.w; ++bit)
                if ((v >> bit) & 1)
                    word.toggle(j * e.w + bit);
        }
        std::vector<int> expected;
        std::vector<std::uint64_t> ring(e.m);
        for (int i = 0; i < e.m; ++i)
            for (int j = 0; j < e.k; ++j)
                ring[i] ^= multiply(e.lifted[i][j], data[j], e.n, e.modulus);
        for (int i = 0; i < e.m; ++i) {
            std::uint64_t want = 0, projected = 0;
            for (int j = 0; j < e.k; ++j)
                want ^= multiply(e.a[i][j], data[j], e.w, e.p);
            for (int j = 0; j < e.m; ++j)
                projected ^= multiply(e.post[i][j], remainder(ring[j], e.p), e.w, e.p);
            require(want == projected, "random quotient-ring arithmetic");
            for (int bit = 0; bit < e.w; ++bit)
                expected.push_back((want >> bit) & 1);
        }
        require(c.evaluate(word) == expected, "random circuit result");
    }
}
Json check_histogram(const Json &data, const Experiment &e, const std::string &hist_key,
                     const std::string &seed_key, const std::string &circuit_key) {
    std::map<int, int> histogram;
    for (int seed = 0; seed < 100; ++seed) {
        auto generated = synthesize(e.direct, e.w * e.k, seed);
        ++histogram[generated.cost()];
        if (seed == integer(data, seed_key))
            require(generated.same(Circuit::from_json(data[circuit_key])),
                    "historical CSE selected circuit mismatch");
    }
    require(histogram.size() == data[hist_key].obj.size(), "CSE histogram support mismatch");
    auto report = Json::object();
    for (auto [cost, count] : histogram) {
        require(data[hist_key][std::to_string(cost)].integer() == count,
                "historical CSE histogram mismatch");
        report[std::to_string(cost)] = count;
    }
    return report;
}
Json verify_small(const Json &data, const Json &comparisons, bool exhaustive) {
    auto e = build_case(data);
    require(e.w == 4 && e.p == 31 && e.modulus == 33, "small example definition");
    auto points = numbers(data["points"]);
    require(points.size() == 6, "small Cauchy point count");
    require(cauchy({points[0], points[1]}, {points[2], points[3], points[4], points[5]}, 4, 31) ==
                e.a,
            "small Cauchy matrix");
    require(check_mds(e.a, 4, 31) == 14, "small MDS minors");
    for (auto key : {"direct", "direct_default", "ring"})
        verify_named(data, key, e.direct);
    for (auto &[name, c] : data["uber_certificates"].obj)
        Circuit::from_json(c).verify(e.direct);
    auto values = Circuit::from_json(data["ring"]).verify(e.direct);
    auto middle_nodes = data["ring"]["middle_outputs"].arr;
    require(middle_nodes.size() == e.middle.size(), "small middle outputs");
    for (std::size_t i = 0; i < middle_nodes.size(); ++i)
        require(values.at(middle_nodes[i].integer()) == e.middle[i], "small middle masks");
    auto histogram = check_histogram(data, e, "direct_histogram", "direct_seed", "direct");
    GreedyCircuit greedy(16);
    auto z = greedy.synthesize_masks(e.middle);
    int main = greedy.ops.size();
    std::vector<std::vector<int>> expressions;
    for (int i = 0; i < 2; ++i)
        for (int bit = 0; bit < 4; ++bit)
            expressions.push_back({z[5 * i + bit], z[5 * i + 4]});
    auto out = greedy.synthesize(expressions);
    auto rebuilt = greedy.export_circuit(out);
    require(rebuilt.same(Circuit::from_json(data["ring"])) && main == 27 && rebuilt.cost() == 35,
            "small staged CSE reproduction");
    int best_direct = std::numeric_limits<int>::max(), best_ring = best_direct;
    for (auto &record : comparisons["records"].arr)
        if (record["status"].string() == "verified") {
            bool direct = record["path"].string() == "direct";
            auto c = Circuit::from_json(record["circuit"]);
            c.verify(direct ? e.direct : e.middle);
            auto full = Circuit::from_json(record["full_circuit"]);
            full.verify(e.direct);
            require(full.cost() == integer(record, "total_cost"), "small comparison total");
            if (!direct)
                require(full.cost() == c.cost() + 8, "small projection cost");
            (direct ? best_direct : best_ring) =
                std::min(direct ? best_direct : best_ring, full.cost());
        }
    Json minors = Json::array();
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j)
            minors.arr.emplace_back(multiply(e.a[0][i], e.a[1][j], 4, 31) ^
                                    multiply(e.a[0][j], e.a[1][i], 4, 31));
    auto uber = Json::object();
    for (auto &[name, c] : data["uber_certificates"].obj)
        uber[name] = c["cost"];
    auto report = Json::object({{"matrix_identity", true},
                                {"all_saved_circuits_verified", true},
                                {"cauchy_matrix", true},
                                {"nonzero_2x2_minors", minors},
                                {"direct_cse_default", data["direct_default"]["cost"]},
                                {"direct_cse_100_seeds_best", 36},
                                {"direct_cse_histogram", histogram},
                                {"uber_costs", uber},
                                {"ring_main", main},
                                {"ring_projection", 8},
                                {"ring_total", 35},
                                {"xsets_best_direct", best_direct},
                                {"xsets_best_ring_with_projection", best_ring}});
    if (exhaustive) {
        int distance = 7;
        auto ring = Circuit::from_json(data["ring"]), direct = Circuit::from_json(data["direct"]);
        for (unsigned word = 0; word < 65536; ++word) {
            std::vector<int> bits;
            int weight = 0;
            for (int j = 0; j < 4; ++j)
                weight += ((word >> (4 * j)) & 15) != 0;
            for (int i = 0; i < 2; ++i) {
                std::uint64_t expected = 0, actual = 0;
                for (int j = 0; j < 4; ++j) {
                    auto v = (word >> (4 * j)) & 15;
                    expected ^= multiply(e.a[i][j], v, 4, 31);
                    actual ^= multiply(e.lifted[i][j], v, 5, 33);
                }
                require(remainder(actual, 31) == expected, "small exhaustive ring result");
                weight += expected != 0;
                for (int bit = 0; bit < 4; ++bit)
                    bits.push_back((expected >> bit) & 1);
            }
            require(ring.evaluate(Mask(word)) == bits && direct.evaluate(Mask(word)) == bits,
                    "small exhaustive circuit");
            if (word)
                distance = std::min(distance, weight);
        }
        require(distance == 3, "small code minimum distance");
        report["exhaustive_inputs"] = 65536;
        report["systematic_code"] = Json::array({6, 4, distance});
    }
    return report;
}
Json verify_original(const Json &data, int samples) {
    auto e = build_case(data);
    require(e.w == 8 && e.n == 17 && e.k == 8 && e.m == 3, "original example shape");
    for (int i = 0; i < e.m; ++i)
        for (int j = 0; j < e.k; ++j) {
            require(e.lifted[i][j] == (std::uint64_t(1) << (i * j)), "original monomial lifts");
            require((e.a[i][j] ^ polynomial_product(e.p, data["G"][i][j].integer())) ==
                        e.lifted[i][j],
                    "original lift certificate");
        }
    auto hist = check_histogram(data, e, "direct_cse_histogram", "direct_cse_seed", "direct_cse");
    require(hist.has("217"), "original CSE bound");
    for (auto key : {"ring", "direct_cse", "ring_cse"})
        verify_named(data, key, e.direct);
    auto main = synthesize(e.middle, e.w * e.k);
    require(main.cost() == 152, "original main CSE");
    require(rows_json(data["projection_input_masks"]) == e.virtuals, "original virtual inputs");
    auto projection = Circuit::from_json(data["projection"]);
    projection.verify(e.projection);
    auto rebuilt = compose(main, projection, e.virtuals, e.middle, e.direct);
    require(rebuilt.circuit.same(Circuit::from_json(data["ring"])) && rebuilt.circuit.cost() == 203,
            "original full circuit reproduction");
    auto middle = Circuit::from_json(data["ring"]).verify(e.direct);
    for (std::size_t i = 0; i < e.middle.size(); ++i)
        require(middle.at(data["ring"]["middle_outputs"][i].integer()) == e.middle[i],
                "original middle output");
    auto baselines = Json::object();
    for (auto &[method, r] : data["direct_xsets"].obj) {
        verify_named(r, "circuit", e.direct);
        baselines[method] = r["circuit"]["cost"];
    }
    for (auto &[name, r] : data["direct_uber"].obj)
        if (r["status"].string() == "verified")
            verify_named(r, "circuit", e.direct);
    for (auto key : {"ring_mw_pipeline", "ring_uber_pipeline"}) {
        const auto &r = data[key];
        require(rows_json(r["projection_input_masks"]) == e.virtuals, "pipeline virtual inputs");
        auto cm = Circuit::from_json(r["main"]), cp = Circuit::from_json(r["projection"]);
        cm.verify(e.middle);
        cp.verify(e.projection);
        auto full = Circuit::from_json(r["ring"]);
        full.verify(e.direct);
        require(compose(cm, cp, e.virtuals, e.middle, e.direct).circuit.same(full),
                "pipeline composition");
        auto vals = full.verify(e.direct);
        for (std::size_t i = 0; i < e.middle.size(); ++i)
            require(vals.at(r["ring"]["middle_outputs"][i].integer()) == e.middle[i],
                    "pipeline middle outputs");
    }
    auto minors = check_mds(e.a, e.w, e.p);
    require(minors == 164, "original MDS count");
    random_vectors(e, rebuilt.circuit, samples);
    return Json::object(
        {{"matrix_identity", true},
         {"all_saved_circuits_verified", true},
         {"main_xors", 152},
         {"projection_xors", 51},
         {"ring_total", 203},
         {"direct_cse_100_seeds", 217},
         {"direct_xsets", baselines},
         {"complete_ring_pipelines", Json::object({{"MW_SQ", 203}, {"UBER_XSET", 204}})},
         {"checked_nonzero_minors", minors},
         {"systematic_code", Json::array({11, 8, 4})},
         {"random_data_vectors", samples}});
}
void check_factored(const Json &r, const Json &original, const Experiment &e) {
    if (r.has("factored_baselines") && !r["factored_baselines"].arr.empty()) {
        require(e.k == 8 && e.m == 3 && e.p == std::uint64_t(integer(original, "P")),
                "factored original dimensions");
        auto old = matrix_json(original["A"]);
        std::reverse(old.begin(), old.end());
        auto shifts = numbers(r["case"]["row_shifts"]);
        Matrix post(e.m, std::vector<std::uint64_t>(e.m));
        for (int i = 0; i < e.m; ++i) {
            post[i][i] = remainder(std::uint64_t(1) << shifts[i], e.p);
            for (int j = 0; j < e.k; ++j)
                require(e.a[i][j] == multiply(post[i][i], old[i][j], e.w, e.p),
                        "factored field identity");
        }
        auto core = Circuit::from_json(original["direct_xsets"]["UBER_XSET"]["circuit"]);
        auto outs = core.outputs;
        for (int i = 0; i < 3; ++i)
            for (int bit = 0; bit < 8; ++bit)
                core.outputs[i * 8 + bit] = outs[(2 - i) * 8 + bit];
        auto ct = expand(old, e.w, e.w, e.p);
        core.verify(ct);
        auto pt = expand(post, e.w, e.w, e.p);
        for (auto &f : r["factored_baselines"].arr)
            if (f.has("circuit")) {
                auto pc = Circuit::from_json(f["post"]);
                pc.verify(pt);
                auto fc = Circuit::from_json(f["circuit"]);
                fc.verify(e.direct);
                require(compose(core, pc, ct, ct, e.direct).circuit.same(fc),
                        "factored baseline reconstruction");
            }
    }
    if (r.has("field_decomposition")) {
        auto &f = r["field_decomposition"];
        auto pre = matrix_json(f["matrices"][0]), core = matrix_json(f["matrices"][1]),
             post = matrix_json(f["matrices"][2]);
        dimensions(pre, e.k, e.k);
        dimensions(post, e.m, e.m);
        for (int i = 0; i < e.k; ++i)
            for (int j = 0; j < e.k; ++j)
                if (i != j)
                    require(pre[i][j] == 0, "diagonal preconditioner");
        for (int i = 0; i < e.m; ++i)
            for (int j = 0; j < e.m; ++j)
                if (i != j)
                    require(post[i][j] == 0, "diagonal postconditioner");
        require(matrix_product(post, matrix_product(core, pre, e.w, e.p), e.w, e.p) == e.a,
                "field decomposition identity");
        std::vector<Circuit> stages;
        require(f["stages"].arr.size() == 3, "field stage count");
        for (std::size_t i = 0; i < 3; ++i) {
            auto stage = Circuit::from_json(f["stages"][i]);
            stage.verify(expand(matrix_json(f["matrices"][i]), e.w, e.w, e.p));
            stages.push_back(stage);
        }
        auto full = Circuit::from_json(f["circuit"]);
        full.verify(e.direct);
        require(compose_layers(stages, e.direct).same(full), "field decomposition reconstruction");
    }
}
} // namespace
Json verify_record(const Json &r, int samples, const Json *original) {
    require(samples >= 0, "negative sample count");
    auto e = build_case(r);
    if (r.has("factored_baselines") && !r["factored_baselines"].arr.empty())
        require(original != nullptr,
                "factored baseline needs the original certificate in --data-dir");
    check_factored(r, original ? *original : Json(), e);
    require(rows_json(r["projection_input_masks"]) == e.virtuals, "projection virtual input masks");
    auto main = Circuit::from_json(r["main"]), projection = Circuit::from_json(r["projection"]),
         ring = Circuit::from_json(r["ring"]);
    main.verify(e.middle);
    projection.verify(e.projection);
    ring.verify(e.direct);
    auto rebuilt = compose(main, projection, e.virtuals, e.middle, e.direct);
    require(rebuilt.circuit.same(ring), "composed circuit differs from certificate");
    verify_named(r, "direct_cse", e.direct);
    int baseline = integer(r["direct_cse"], "cost");
    for (auto &measurement : r["measurements"].arr) {
        check_measurement(measurement, e);
        if (measurement["status"].string() == "verified" &&
            measurement["path"].string() == "direct")
            baseline = std::min(baseline, integer(measurement, "cost"));
    }
    if (r.has("factored_baselines"))
        for (auto &f : r["factored_baselines"].arr)
            if (f.has("circuit"))
                baseline = std::min(baseline, integer(f["circuit"], "cost"));
    if (r.has("field_decomposition"))
        baseline = std::min(baseline, integer(r["field_decomposition"]["circuit"], "cost"));
    require(baseline == integer(r, "best_completed_direct"),
            "best completed field baseline mismatch");
    int minors = check_mds(e.a, e.w, e.p);
    require(minors == integer(r, "checked_nonzero_minors"), "MDS minor count mismatch");
    random_vectors(e, ring, samples);
    bool practical = r.has("kind");
    if (practical) {
        auto status = r.has("comparison_status") ? r["comparison_status"].string()
                                                 : "completed_heuristic_comparison";
        return Json::object(
            {{"kind", r["kind"]},
             {"comparison_status", status},
             {"systematic_code", Json::array({e.k + e.m, e.k, e.m + 1})},
             {"main_xors", rebuilt.main_cost},
             {"projection_xors", rebuilt.projection_cost},
             {"ring", ring.cost()},
             {"best_completed_direct", baseline},
             {"net_saved",
              status == "completed_heuristic_comparison" ? Json(baseline - ring.cost()) : Json()},
             {"nonzero_minors", minors},
             {"random_vectors", samples}});
    }
    bool vandermonde = is_vandermonde(e.a, e.w, e.p);
    std::string family = vandermonde ? "vandermonde" : "scaled_extension";
    if (r.has("matrix_family"))
        require(r["matrix_family"].string() == family, "matrix family label mismatch");
    double percent = 100. * (baseline - ring.cost()) / baseline;
    require(std::abs(r["percent"].number() - percent) < 1e-10, "gain percentage mismatch");
    return Json::object({{"matrix_family", family},
                         {"ring", ring.cost()},
                         {"best_completed_direct", baseline},
                         {"percent", r["percent"]},
                         {"nonzero_minors", minors},
                         {"random_vectors", samples}});
}
Json verify_file(const Json &data, int samples, const Json *original) {
    if (data.type == Json::Type::Array) {
        auto reports = Json::array();
        for (auto &r : data.arr)
            reports.arr.push_back(verify_record(r, samples, original));
        return reports;
    }
    if (data.has("target_masks") && data.has("status")) {
        require(data["status"].string() == "verified", "matrix analysis was not completed");
        Circuit::from_json(data["circuit"]).verify(rows_json(data["target_masks"]));
        return Json::object({{"status", "verified"}, {"symbolic_verification", true}});
    }
    if (data.has("records") && data.has("definition")) {
        auto e = build_case(data["definition"]);
        for (const auto &record : data["records"].arr) {
            if (record["status"].string() != "verified") {
                require(record["status"].string() == "timeout" ||
                            record["status"].string() == "cancelled",
                        "comparison contains an execution error");
                continue;
            }
            auto c = Circuit::from_json(record["circuit"]);
            c.verify(record["path"].string() == "direct" ? e.direct : e.middle);
            auto full = Circuit::from_json(record["full_circuit"]);
            full.verify(e.direct);
            require(full.cost() == record["total_cost"].integer(), "comparison total cost");
        }
        return Json::object({{"all_saved_circuits_verified", true}});
    }
    if (data.has("path") && data.has("status")) {
        auto e = build_case(data["definition"]);
        check_measurement(data, e);
        return Json::object({{"status", data["status"]}, {"symbolic_verification", true}});
    }
    if (data.has("examples")) {
        auto reports = Json::object();
        for (auto &[name, r] : data["examples"].obj)
            reports[name] = verify_record(r, samples, original);
        return reports;
    }
    return verify_record(data, samples, original);
}
Json verify_all(const std::string &dir, const std::string &group, int samples, bool exhaustive) {
    require(samples >= 0, "negative sample count");
    auto report = Json::object();
    auto read = [&](const std::string &name) { return Json::read(dir + "/" + name + ".json"); };
    if (group == "all" || group == "small")
        report["small"] =
            verify_small(read("ring_lift_example"), read("xsets_comparison"), exhaustive);
    if (group == "all" || group == "original")
        report["original"] = verify_original(read("structured_ring_example"), samples);
    if (group == "all" || group == "family") {
        auto data = read("structured_ring_family"), original = read("structured_ring_example");
        auto e = build_case(data["primary_case"]);
        require(is_vandermonde(e.a, e.w, e.p), "primary matrix must be Vandermonde");
        for (auto &r : data["primary_additional_runs"].arr) {
            check_measurement(r, e);
            if (r.has("ring"))
                verify_named(r, "ring", e.direct);
        }
        auto examples = Json::object();
        for (auto &[name, r] : data["examples"].obj) {
            examples[name] = verify_record(r, samples, &original);
        }
        report["family"] = Json::object(
            {{"primary_additional_runs", int(data["primary_additional_runs"].arr.size())},
             {"examples", examples}});
    }
    if (group == "all" || group == "practical")
        report["practical"] = verify_file(read("practical_mds_examples"), samples);
    require(!report.obj.empty(), "unknown verification group");
    return report;
}
} // namespace ringlab
