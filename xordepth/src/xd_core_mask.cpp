#include "xd/core.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace xd {
BitMask::BitMask(std::uint64_t v) {
    if (v)
        words.push_back(v);
}
void BitMask::trim() {
    while (!words.empty() && !words.back())
        words.pop_back();
}
BitMask BitMask::unit(std::size_t b) {
    BitMask m;
    m.toggle(b);
    return m;
}
void BitMask::toggle(std::size_t b) {
    require(b < 1000000, "mask too wide");
    if (words.size() <= b / 64)
        words.resize(b / 64 + 1);
    words[b / 64] ^= std::uint64_t(1) << (b % 64);
    trim();
}
bool BitMask::bit(std::size_t b) const {
    return b / 64 < words.size() && ((words[b / 64] >> (b % 64)) & 1);
}
std::size_t BitMask::bits() const {
    return words.empty() ? 0 : 64 * (words.size() - 1) + 64 - __builtin_clzll(words.back());
}
int BitMask::count() const {
    int n = 0;
    for (auto w : words)
        n += __builtin_popcountll(w);
    return n;
}
BitMask &BitMask::operator^=(const BitMask &b) {
    words.resize(std::max(words.size(), b.words.size()));
    for (std::size_t i = 0; i < b.words.size(); ++i)
        words[i] ^= b.words[i];
    trim();
    return *this;
}
bool BitMask::operator<(const BitMask &b) const {
    if (words.size() != b.words.size())
        return words.size() < b.words.size();
    return std::lexicographical_compare(words.rbegin(), words.rend(), b.words.rbegin(),
                                        b.words.rend());
}
std::size_t BitMaskHash::operator()(const BitMask &m) const {
    std::size_t h = m.words.size();
    for (auto w : m.words)
        h ^= std::hash<std::uint64_t>{}(w) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}
BitMask BitMask::decimal(const std::string &s) {
    require(!s.empty() && s.size() < 301100, "invalid mask integer");
    BitMask m;
    // Work in base 2^32 so decimal multiplication is strictly standard C++17.
    std::vector<std::uint32_t> limbs;
    for (char c : s) {
        require(c >= '0' && c <= '9', "mask must be a nonnegative integer");
        std::uint64_t carry = c - '0';
        for (auto &limb : limbs) {
            carry += std::uint64_t(limb) * 10;
            limb = std::uint32_t(carry);
            carry >>= 32;
        }
        if (carry)
            limbs.push_back(std::uint32_t(carry));
    }
    m.words.resize((limbs.size() + 1) / 2);
    for (std::size_t i = 0; i < limbs.size(); ++i)
        m.words[i / 2] |= std::uint64_t(limbs[i]) << ((i % 2) * 32);
    m.trim();
    return m;
}
std::string BitMask::decimal() const {
    std::vector<std::uint32_t> limbs;
    for (auto word : words) {
        limbs.push_back(std::uint32_t(word));
        limbs.push_back(std::uint32_t(word >> 32));
    }
    while (!limbs.empty() && !limbs.back())
        limbs.pop_back();
    std::string result;
    while (!limbs.empty()) {
        std::uint64_t rem = 0;
        for (std::size_t i = limbs.size(); i-- > 0;) {
            std::uint64_t v = (rem << 32) | limbs[i];
            limbs[i] = std::uint32_t(v / 10);
            rem = v % 10;
        }
        result += char('0' + rem);
        while (!limbs.empty() && !limbs.back())
            limbs.pop_back();
    }
    if (result.empty())
        return "0";
    std::reverse(result.begin(), result.end());
    return result;
}

} // namespace xd
