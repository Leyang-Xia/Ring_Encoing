#include "ringlab/core.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace ringlab {
static int degree(std::uint64_t v) {
    return v ? 63 - __builtin_clzll(v) : -1;
}
std::uint64_t polynomial_product(std::uint64_t a, std::uint64_t b) {
    require(a == 0 || b == 0 || degree(a) + degree(b) < 64, "polynomial product exceeds 64 bits");
    std::uint64_t v = 0;
    while (b) {
        if (b & 1)
            v ^= a;
        b >>= 1;
        a <<= 1;
    }
    return v;
}
std::uint64_t remainder(std::uint64_t a, std::uint64_t p) {
    require(p > 1, "invalid polynomial modulus");
    int d = degree(p);
    while (a && degree(a) >= d)
        a ^= p << (degree(a) - d);
    return a;
}
std::uint64_t cyclic_multiply(std::uint64_t a, std::uint64_t b, int width) {
    require(width > 0 && width < 63, "invalid cyclic ring width");
    const auto mask = (std::uint64_t(1) << width) - 1;
    require(a <= mask && b <= mask, "cyclic operand out of range");
    std::uint64_t result = 0;
    for (int shift = 0; shift < width; ++shift)
        if ((a >> shift) & 1)
            result ^= ((b << shift) | (b >> (width - shift))) & mask;
    return result;
}
std::uint64_t multiply(std::uint64_t a, std::uint64_t b, int width, std::uint64_t p) {
    require(width > 0 && width < 63 && degree(p) == width, "invalid quotient width or modulus");
    a = remainder(a, p);
    b = remainder(b, p);
    std::uint64_t v = 0, top = std::uint64_t(1) << width;
    while (b) {
        if (b & 1)
            v ^= a;
        b >>= 1;
        a <<= 1;
        if (a & top)
            a ^= p;
    }
    return v;
}
std::uint64_t power(std::uint64_t a, std::uint64_t exponent, int w, std::uint64_t p) {
    std::uint64_t v = 1;
    while (exponent) {
        if (exponent & 1)
            v = multiply(v, a, w, p);
        a = multiply(a, a, w, p);
        exponent >>= 1;
    }
    return v;
}
bool irreducible(std::uint64_t p, int w) {
    if (w < 1 || w >= 63 || degree(p) != w)
        return false;
    auto gcd = [](std::uint64_t a, std::uint64_t b) {
        while (b) {
            auto r = b == 1 ? 0 : remainder(a, b);
            a = b;
            b = r;
        }
        return a;
    };
    auto x = remainder(2, p), v = x;
    for (int i = 1; i <= w; ++i) {
        v = multiply(v, v, w, p);
        if (i <= w / 2 && gcd(v ^ x, p) != 1)
            return false;
    }
    return v == x;
}
Rows expand(const Matrix &a, int iw, int ow, std::uint64_t p) {
    require(!a.empty() && !a[0].empty(), "empty coefficient matrix");
    require(iw > 0 && iw < 63 && ow > 0 && ow < 63, "invalid expansion width");
    Rows rows(a.size() * ow);
    for (std::size_t i = 0; i < a.size(); ++i) {
        require(a[i].size() == a[0].size(), "ragged matrix");
        for (std::size_t j = 0; j < a[i].size(); ++j)
            for (int b = 0; b < iw; ++b) {
                auto col = multiply(a[i][j], std::uint64_t(1) << b, ow, p);
                for (int bit = 0; bit < ow; ++bit)
                    if ((col >> bit) & 1)
                        rows[i * ow + bit].toggle(j * iw + b);
            }
    }
    return rows;
}
Matrix matrix_json(const Json &j) {
    require(j.type == Json::Type::Array, "expected matrix array");
    Matrix a;
    for (auto &row : j.arr) {
        require(row.type == Json::Type::Array, "expected matrix row");
        std::vector<std::uint64_t> r;
        for (auto &x : row.arr) {
            auto n = x.integer();
            require(n >= 0, "negative coefficient");
            r.push_back(n);
        }
        a.push_back(std::move(r));
    }
    return a;
}
Json to_json(const Matrix &a) {
    auto j = Json::array();
    for (auto &r : a) {
        auto row = Json::array();
        for (auto x : r)
            row.arr.emplace_back(x);
        j.arr.push_back(std::move(row));
    }
    return j;
}
Rows rows_json(const Json &j) {
    require(j.type == Json::Type::Array, "expected masks array");
    Rows rows;
    for (auto &x : j.arr) {
        require(x.type == Json::Type::Number, "expected numeric mask");
        rows.push_back(Mask::decimal(x.text));
    }
    return rows;
}
Json to_json(const Rows &rows) {
    auto j = Json::array();
    for (auto &mask : rows)
        j.arr.push_back(Json::numeric(mask.decimal()));
    return j;
}
Matrix matrix_product(const Matrix &a, const Matrix &b, int w, std::uint64_t p) {
    require(!a.empty() && !b.empty() && !b[0].empty(), "empty matrix product");
    for (auto &r : a)
        require(r.size() == b.size(), "matrix product dimensions");
    for (auto &r : b)
        require(r.size() == b[0].size(), "ragged matrix product");
    Matrix c(a.size(), std::vector<std::uint64_t>(b[0].size()));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t k = 0; k < b.size(); ++k)
            for (std::size_t j = 0; j < b[0].size(); ++j)
                c[i][j] ^= multiply(a[i][k], b[k][j], w, p);
    return c;
}
Matrix inverse_matrix(const Matrix &a, int w, std::uint64_t p) {
    require(!a.empty(), "empty inverse");
    int n = a.size();
    Matrix b = a;
    for (int i = 0; i < n; ++i) {
        require(int(b[i].size()) == n, "nonsquare inverse");
        b[i].resize(2 * n);
        b[i][n + i] = 1;
    }
    for (int j = 0; j < n; ++j) {
        int pivot = j;
        while (pivot < n && !b[pivot][j])
            ++pivot;
        require(pivot < n, "singular matrix");
        std::swap(b[j], b[pivot]);
        auto inv = power(b[j][j], (std::uint64_t(1) << w) - 2, w, p);
        for (auto &x : b[j])
            x = multiply(x, inv, w, p);
        for (int i = 0; i < n; ++i)
            if (i != j) {
                auto coefficient = b[i][j];
                for (int k = 0; k < 2 * n; ++k)
                    b[i][k] ^= multiply(coefficient, b[j][k], w, p);
            }
    }
    Matrix out;
    for (auto &row : b)
        out.emplace_back(row.begin() + n, row.end());
    return out;
}
bool full_rank(Matrix a, int w, std::uint64_t p) {
    int n = a.size();
    for (auto &r : a)
        require(int(r.size()) == n, "rank requires square matrix");
    for (int j = 0; j < n; ++j) {
        int pivot = j;
        while (pivot < n && !a[pivot][j])
            ++pivot;
        if (pivot == n)
            return false;
        std::swap(a[j], a[pivot]);
        auto inv = power(a[j][j], (std::uint64_t(1) << w) - 2, w, p);
        for (auto &x : a[j])
            x = multiply(x, inv, w, p);
        for (int i = j + 1; i < n; ++i) {
            auto f = a[i][j];
            for (int k = j; k < n; ++k)
                a[i][k] ^= multiply(f, a[j][k], w, p);
        }
    }
    return true;
}
static void combinations(int n, int k, const std::function<void(const std::vector<int> &)> &f) {
    std::vector<int> v;
    std::function<void(int)> visit = [&](int start) {
        if (int(v.size()) == k) {
            f(v);
            return;
        }
        for (int i = start; i <= n - (k - int(v.size())); ++i) {
            v.push_back(i);
            visit(i + 1);
            v.pop_back();
        }
    };
    visit(0);
}
int check_mds(const Matrix &a, int w, std::uint64_t p) {
    require(!a.empty() && !a[0].empty(), "empty MDS matrix");
    for (auto &r : a)
        require(r.size() == a[0].size(), "ragged MDS matrix");
    int count = 0;
    for (int size = 1; size <= int(std::min(a.size(), a[0].size())); ++size)
        combinations(a.size(), size, [&](auto &rows) {
            combinations(a[0].size(), size, [&](auto &cols) {
                Matrix sub;
                for (auto i : rows) {
                    std::vector<std::uint64_t> row;
                    for (auto j : cols)
                        row.push_back(a[i][j]);
                    sub.push_back(row);
                }
                require(full_rank(sub, w, p), "non-MDS matrix: singular minor");
                ++count;
            });
        });
    return count;
}

} // namespace ringlab
