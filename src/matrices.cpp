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
Matrix identity(int n) {
    Matrix a(n, std::vector<std::uint64_t>(n));
    for (int i = 0; i < n; ++i)
        a[i][i] = 1;
    return a;
}
void dimensions(const Matrix &a, int rows, int cols) {
    require(int(a.size()) == rows, "coefficient row count");
    for (auto &r : a)
        require(int(r.size()) == cols, "coefficient column count");
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
Matrix source_matrix(const Json &j, int w, std::uint64_t p, int k, int m) {
    auto kind = j["kind"].string();
    auto &definition = j["definition"];
    if (kind == "cauchy")
        return cauchy(numbers(definition["x"]), numbers(definition["y"]), w, p);
    require(kind == "systematic_rs" || kind == "systematic_rs_moments", "unknown matrix kind");
    auto points = numbers(definition["points"]);
    require(int(points.size()) == k + m &&
                std::set<std::uint64_t>(points.begin(), points.end()).size() == points.size(),
            "RS evaluation points");
    Matrix v;
    for (auto point : points) {
        require(point < (std::uint64_t(1) << w), "RS point out of range");
        std::vector<std::uint64_t> row;
        for (int j = 0; j < k; ++j)
            row.push_back(power(point, j, w, p));
        v.push_back(row);
    }
    Matrix top(v.begin(), v.begin() + k);
    auto g = matrix_product(v, inverse_matrix(top, w, p), w, p);
    require(Matrix(g.begin(), g.begin() + k) == identity(k), "RS systematic block");
    Matrix a(g.begin() + k, g.end());
    std::vector<std::uint64_t> data(points.begin(), points.begin() + k),
        parity(points.begin() + k, points.end());
    auto c = cauchy(parity, data, w, p);
    std::vector<std::uint64_t> rw, cw;
    for (auto b : parity) {
        std::uint64_t v = 1;
        for (auto d : data)
            v = multiply(v, d ^ b, w, p);
        rw.push_back(v);
    }
    for (int i = 0; i < k; ++i) {
        std::uint64_t v = 1;
        for (int j = 0; j < k; ++j)
            if (i != j)
                v = multiply(v, data[i] ^ data[j], w, p);
        cw.push_back(power(v, (std::uint64_t(1) << w) - 2, w, p));
    }
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < k; ++j)
            require(a[i][j] == multiply(rw[i], multiply(c[i][j], cw[j], w, p), w, p),
                    "RS Lagrange identity");
    return a;
}
void project_matrices(Experiment &e, const Json *saved_order) {
    e.direct = expand(e.a, e.w, e.w, e.p);
    e.middle = expand(e.lifted, e.w, e.n, e.modulus);
    std::unordered_map<Mask, int, MaskHash> index;
    for (auto &mask : e.middle)
        if (!mask.empty() && !index.count(mask)) {
            index[mask] = e.virtuals.size();
            e.virtuals.push_back(mask);
        }
    // The original certificate puts free unit signals first. Later experiments
    // use first-occurrence order. Both are valid, but gate indices depend on it.
    if (saved_order) {
        auto requested = rows_json(*saved_order);
        std::set<Mask> unique(requested.begin(), requested.end());
        require(unique.size() == requested.size() && requested.size() == e.virtuals.size(),
                "projection input order has duplicates or missing signals");
        for (const auto &mask : requested)
            require(index.count(mask), "projection input is not a middle output");
        e.virtuals = std::move(requested);
        index.clear();
        for (std::size_t i = 0; i < e.virtuals.size(); ++i)
            index[e.virtuals[i]] = i;
    }
    for (int output = 0; output < e.m; ++output)
        for (int bit = 0; bit < e.w; ++bit) {
            Mask row;
            for (int i = 0; i < e.m; ++i)
                for (int j = 0; j < e.n; ++j) {
                    auto mask = e.middle[i * e.n + j];
                    auto coefficient = multiply(e.post[output][i],
                                                remainder(std::uint64_t(1) << j, e.p), e.w, e.p);
                    if (!mask.empty() && ((coefficient >> bit) & 1))
                        row.toggle(index.at(mask));
                }
            e.projection.push_back(std::move(row));
        }
    for (std::size_t i = 0; i < e.projection.size(); ++i) {
        Mask composed;
        for (std::size_t j = 0; j < e.virtuals.size(); ++j)
            if (e.projection[i].bit(j))
                composed ^= e.virtuals[j];
        require(composed == e.direct[i], "field/ring/projection matrix identity");
    }
}
} // namespace
Experiment build_case(const Json &input) {
    const auto &j = input.has("case") ? input["case"] : input;
    Experiment e;
    e.w = integer(j, "w");
    e.p = integer(j, "P");
    require(e.w > 0 && e.w <= 16 && irreducible(e.p, e.w),
            "unsupported field width (1..16) or reducible polynomial");
    bool family = j.has("exponents");
    if (j.has("kind")) {
        auto k = integer(j, "k"), m = integer(j, "m");
        require(k > 0 && m > 0 && k <= (1 << e.w) && m <= (1 << e.w) - k,
                "invalid code dimensions");
    }
    if (family) {
        e.n = integer(j, "n");
        e.k = integer(j, "k");
        e.m = integer(j, "m");
        require(e.n > e.w && e.n < 63 && e.k > 0 && e.m > 0, "invalid family dimensions");
        e.modulus = j.has("ring_modulus") ? std::uint64_t(j["ring_modulus"].integer())
                                          : (std::uint64_t(1) << e.n) | 1;
        auto exponents = numbers(j["exponents"]);
        require(int(exponents.size()) == e.k, "column exponent count");
        std::vector<std::uint64_t> slopes(e.m), shifts(e.m), columns(e.k);
        std::iota(slopes.begin(), slopes.end(), 0);
        if (j.has("slopes"))
            slopes = numbers(j["slopes"]);
        if (j.has("row_shifts"))
            shifts = numbers(j["row_shifts"]);
        if (j.has("column_shifts"))
            columns = numbers(j["column_shifts"]);
        require(int(slopes.size()) == e.m && int(shifts.size()) == e.m &&
                    int(columns.size()) == e.k,
                "shift dimensions");
        for (int i = 0; i < e.m; ++i) {
            std::vector<std::uint64_t> a, lift;
            for (int col = 0; col < e.k; ++col) {
                auto exponent = (slopes[i] * exponents[col] + shifts[i] + columns[col]) % e.n;
                auto v = std::uint64_t(1) << exponent;
                lift.push_back(v);
                a.push_back(remainder(v, e.p));
            }
            e.a.push_back(a);
            e.lifted.push_back(lift);
        }
    } else {
        e.a = j.has("kind") ? source_matrix(j, e.w, e.p, integer(j, "k"), integer(j, "m"))
                            : matrix_json(j["A"]);
        require(!e.a.empty() && !e.a[0].empty(), "empty field matrix");
        e.k = e.a[0].size();
        e.m = e.a.size();
        e.n = j.has("n") ? integer(j, "n") : 5;
        require(e.n > e.w && e.n < 63, "invalid ring dimension");
        e.modulus = j.has("M") ? std::uint64_t(j["M"].integer())
                               : polynomial_product(e.p, j["Q"].integer());
        e.lifted = matrix_json(j.has("lifted") ? j["lifted"] : j["ring"]["lifted"]);
    }
    dimensions(e.a, e.m, e.k);
    dimensions(e.lifted, e.m, e.k);
    if (j.has("k"))
        require(e.k == integer(j, "k"), "declared data dimension mismatch");
    if (j.has("m"))
        require(e.m == integer(j, "m"), "declared parity dimension mismatch");
    require(e.modulus > 1 && Mask(e.modulus).bits() == std::size_t(e.n + 1) &&
                remainder(e.modulus, e.p) == 0,
            "ring modulus must have degree n and be divisible by P");
    e.post = j.has("post_field_matrix") ? matrix_json(j["post_field_matrix"]) : identity(e.m);
    dimensions(e.post, e.m, e.m);
    Matrix field_main = j.has("main_field_matrix") ? matrix_json(j["main_field_matrix"]) : e.a;
    dimensions(field_main, e.m, e.k);
    require(matrix_product(e.post, field_main, e.w, e.p) == e.a, "field postprocessing identity");
    for (int i = 0; i < e.m; ++i)
        for (int col = 0; col < e.k; ++col) {
            require(e.a[i][col] < (std::uint64_t(1) << e.w) &&
                        e.lifted[i][col] < (std::uint64_t(1) << e.n),
                    "coefficient representation out of range");
            require(remainder(e.lifted[i][col], e.p) == field_main[i][col],
                    "lifted coefficient projection");
            if (e.modulus == ((std::uint64_t(1) << e.n) | 1))
                for (int bit = 0; bit < e.w; ++bit)
                    require(cyclic_multiply(e.lifted[i][col], std::uint64_t(1) << bit, e.n) ==
                                multiply(e.lifted[i][col], std::uint64_t(1) << bit, e.n, e.modulus),
                            "independent cyclic multiplication matrix mismatch");
        }
    if (j.has("kind") && j["kind"].string() == "systematic_rs_moments") {
        auto points = numbers(j["definition"]["points"]);
        Matrix left(e.m), right(e.m);
        for (int i = 0; i < e.m; ++i)
            for (int col = 0; col < e.k + e.m; ++col)
                (col < e.k ? left[i] : right[i]).push_back(power(points[col], i + 1, e.w, e.p));
        require(left == field_main && inverse_matrix(right, e.w, e.p) == e.post,
                "systematic RS moment identity");
    }
    if (input.has("A"))
        require(matrix_json(input["A"]) == e.a, "saved source matrix differs from reconstruction");
    if (input.has("lifted"))
        require(matrix_json(input["lifted"]) == e.lifted, "saved lift differs from reconstruction");
    project_matrices(e, input.has("projection_input_masks") ? &input["projection_input_masks"]
                                                            : nullptr);
    return e;
}
bool is_vandermonde(const Matrix &a, int w, std::uint64_t p) {
    if (a.size() < 2 || a[0].empty())
        return false;
    std::set<std::uint64_t> points(a[1].begin(), a[1].end());
    if (points.count(0) || points.size() != a[1].size())
        return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].size() != a[1].size())
            return false;
        for (std::size_t j = 0; j < a[i].size(); ++j)
            if (a[i][j] != power(a[1][j], i, w, p))
                return false;
    }
    return true;
}
} // namespace ringlab
