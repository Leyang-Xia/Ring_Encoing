#include "xd/core.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace xd {

// Growth policy: small masks stay inline; larger requests allocate exactly
// what is needed (the growth paths run once per mask, not in hot loops).
void BitMask::reserve_words(std::uint32_t words) {
    if (words <= capacity_) return;
    const std::uint32_t capacity = words > kInlineWords ? words : kInlineWords;
    std::uint64_t *fresh = capacity > kInlineWords ? new std::uint64_t[capacity]()
                                                   : inline_;
    std::copy(words_, words_ + size_, fresh);
    std::fill(fresh + size_, fresh + capacity, 0);
    if (words_ != inline_) delete[] words_;
    words_ = fresh;
    capacity_ = capacity;
}

BitMask::BitMask(std::uint64_t v) {
    if (v) {
        size_ = 1;
        words_[0] = v;
    }
}

BitMask::BitMask(const BitMask &other)
    : size_(other.size_),
      capacity_(other.size_ > kInlineWords ? other.size_ : kInlineWords) {
    words_ = capacity_ > kInlineWords ? new std::uint64_t[capacity_]() : inline_;
    std::copy(other.words_, other.words_ + size_, words_);
}

BitMask::BitMask(BitMask &&other) noexcept : size_(other.size_) {
    if (other.words_ != other.inline_) {
        words_ = other.words_;
        capacity_ = other.capacity_;
        other.words_ = other.inline_;
        other.capacity_ = kInlineWords;
    } else {
        words_ = inline_;
        capacity_ = kInlineWords;
        std::copy(other.inline_, other.inline_ + kInlineWords, inline_);
    }
    other.size_ = 0;
    // Keep the "words beyond size_ are zero" invariant for the moved-from
    // mask: its inline storage may otherwise hold stale bits.
    std::fill(other.inline_, other.inline_ + kInlineWords, 0);
}

BitMask &BitMask::operator=(const BitMask &other) {
    if (this == &other) return *this;
    if (other.size_ <= capacity_) {
        std::copy(other.words_, other.words_ + other.size_, words_);
        std::fill(words_ + other.size_, words_ + capacity_, 0);
        size_ = other.size_;
    } else {
        BitMask fresh(other);
        *this = std::move(fresh);
    }
    return *this;
}

BitMask &BitMask::operator=(BitMask &&other) noexcept {
    if (this == &other) return *this;
    if (other.words_ != other.inline_) {
        if (words_ != inline_) delete[] words_;
        words_ = other.words_;
        capacity_ = other.capacity_;
    } else {
        if (words_ != inline_) {
            delete[] words_;
            words_ = inline_;
            capacity_ = kInlineWords;
        }
        std::copy(other.inline_, other.inline_ + kInlineWords, inline_);
    }
    size_ = other.size_;
    other.size_ = 0;
    other.words_ = other.inline_;
    other.capacity_ = kInlineWords;
    std::fill(other.inline_, other.inline_ + kInlineWords, 0);
    return *this;
}

BitMask::~BitMask() {
    if (words_ != inline_) delete[] words_;
}

void BitMask::trim() {
    while (size_ > 0 && words_[size_ - 1] == 0) --size_;
}

BitMask BitMask::unit(std::size_t b) {
    BitMask m;
    m.toggle(b);
    return m;
}

void BitMask::toggle(std::size_t b) {
    require(b < 1000000, "mask too wide");
    const std::uint32_t index = static_cast<std::uint32_t>(b / 64);
    if (index >= size_) {
        reserve_words(index + 1);
        size_ = index + 1;
    }
    words_[index] ^= std::uint64_t(1) << (b % 64);
    trim();
}

bool BitMask::bit(std::size_t b) const {
    return b / 64 < size_ && ((words_[b / 64] >> (b % 64)) & 1);
}

std::size_t BitMask::bits() const {
    return size_ == 0 ? 0 : 64 * (size_ - 1) + 64 - __builtin_clzll(words_[size_ - 1]);
}

int BitMask::count() const {
    int n = 0;
    for (std::uint32_t i = 0; i < size_; ++i)
        n += __builtin_popcountll(words_[i]);
    return n;
}

BitMask &BitMask::operator^=(const BitMask &b) {
    if (b.size_ > size_) {
        reserve_words(b.size_);
        size_ = b.size_;
    }
    for (std::uint32_t i = 0; i < b.size_; ++i)
        words_[i] ^= b.words_[i];
    trim();
    return *this;
}

bool BitMask::operator==(const BitMask &b) const {
    if (size_ != b.size_) return false;
    for (std::uint32_t i = 0; i < size_; ++i)
        if (words_[i] != b.words_[i]) return false;
    return true;
}

bool BitMask::operator<(const BitMask &b) const {
    if (size_ != b.size_) return size_ < b.size_;
    for (std::size_t i = size_; i-- > 0;) {
        if (words_[i] != b.words_[i]) return words_[i] < b.words_[i];
    }
    return false;
}

std::size_t BitMaskHash::operator()(const BitMask &m) const {
    std::size_t h = m.size_;
    for (std::uint32_t i = 0; i < m.size_; ++i)
        h ^= std::hash<std::uint64_t>{}(m.words_[i]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
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
    const std::uint32_t words = static_cast<std::uint32_t>((limbs.size() + 1) / 2);
    m.reserve_words(words);
    m.size_ = words;
    for (std::size_t i = 0; i < limbs.size(); ++i)
        m.words_[i / 2] |= std::uint64_t(limbs[i]) << ((i % 2) * 32);
    m.trim();
    return m;
}

std::string BitMask::decimal() const {
    std::vector<std::uint32_t> limbs;
    for (std::uint32_t i = 0; i < size_; ++i) {
        limbs.push_back(std::uint32_t(words_[i]));
        limbs.push_back(std::uint32_t(words_[i] >> 32));
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
