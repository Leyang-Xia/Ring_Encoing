#pragma once
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace ringlab {
// Numbers retain their lexical representation: certificate masks can exceed 2^112.
class Json {
  public:
    enum class Type { Null, Boolean, Number, String, Array, Object };
    Type type = Type::Null;
    std::string text;
    std::vector<Json> arr;
    std::map<std::string, Json> obj;
    Json() = default;
    Json(bool value);
    Json(int value);
    Json(std::int64_t value);
    Json(std::uint64_t value);
    Json(double value);
    Json(const char *value);
    Json(const std::string &value);
    static Json array(std::vector<Json> values = {});
    static Json object(std::map<std::string, Json> values = {});
    static Json numeric(const std::string &value);
    static Json parse(const std::string &input);
    static Json read(const std::string &path);
    std::string dump(int indent = 2) const;
    void write(const std::string &path) const;
    bool has(const std::string &key) const;
    bool null() const {
        return type == Type::Null;
    }
    std::int64_t integer() const;
    int int32() const;
    double number() const;
    const std::string &string() const;
    const Json &operator[](const std::string &key) const;
    Json &operator[](const std::string &key);
    const Json &operator[](std::size_t index) const;
    Json &operator[](std::size_t index);
};
void require(bool condition, const std::string &message);
inline void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
} // namespace ringlab
