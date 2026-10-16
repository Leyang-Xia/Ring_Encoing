#include "ringlab/json.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ringlab {
void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}
Json::Json(bool v) : type(Type::Boolean), text(v ? "true" : "false") {}
Json::Json(int v) : Json(static_cast<std::int64_t>(v)) {}
Json::Json(std::int64_t v) : type(Type::Number), text(std::to_string(v)) {}
Json::Json(std::uint64_t v) : type(Type::Number), text(std::to_string(v)) {}
Json::Json(double v) : type(Type::Number) {
    require(std::isfinite(v), "nonfinite JSON number");
    std::ostringstream s;
    s << std::setprecision(17) << v;
    text = s.str();
}
Json::Json(const char *v) : Json(std::string(v)) {}
Json::Json(const std::string &v) : type(Type::String), text(v) {}
Json Json::array(std::vector<Json> v) {
    Json j;
    j.type = Type::Array;
    j.arr = std::move(v);
    return j;
}
Json Json::object(std::map<std::string, Json> v) {
    Json j;
    j.type = Type::Object;
    j.obj = std::move(v);
    return j;
}
Json Json::numeric(const std::string &v) {
    auto j = parse(v);
    require(j.type == Type::Number, "expected JSON number");
    return j;
}
bool Json::has(const std::string &key) const {
    return type == Type::Object && obj.count(key);
}
std::int64_t Json::integer() const {
    require(type == Type::Number, "expected integer");
    std::size_t n;
    auto v = std::stoll(text, &n);
    require(n == text.size(), "expected integral number");
    return v;
}
double Json::number() const {
    require(type == Type::Number, "expected number");
    return std::stod(text);
}
int Json::int32() const {
    auto value = integer();
    require(value >= std::numeric_limits<int>::min() && value <= std::numeric_limits<int>::max(),
            "integer exceeds supported index range");
    return static_cast<int>(value);
}
const std::string &Json::string() const {
    require(type == Type::String, "expected string");
    return text;
}
const Json &Json::operator[](const std::string &k) const {
    require(type == Type::Object, "expected object for " + k);
    return obj.at(k);
}
Json &Json::operator[](const std::string &k) {
    require(type == Type::Object, "expected object for " + k);
    return obj[k];
}
const Json &Json::operator[](std::size_t i) const {
    require(type == Type::Array, "expected array");
    return arr.at(i);
}
Json &Json::operator[](std::size_t i) {
    require(type == Type::Array, "expected array");
    return arr.at(i);
}
namespace {
struct Parser {
    const std::string &s;
    std::size_t pos = 0;
    void ws() {
        while (pos < s.size() &&
               (s[pos] == ' ' || s[pos] == '\n' || s[pos] == '\r' || s[pos] == '\t'))
            ++pos;
    }
    char take() {
        require(pos < s.size(), "unexpected end of JSON");
        return s[pos++];
    }
    bool consume(char c) {
        ws();
        if (pos < s.size() && s[pos] == c) {
            ++pos;
            return true;
        }
        return false;
    }
    unsigned hex4() {
        unsigned v = 0;
        for (int i = 0; i < 4; ++i) {
            char c = take();
            unsigned d = c >= '0' && c <= '9'   ? c - '0'
                         : c >= 'a' && c <= 'f' ? c - 'a' + 10
                         : c >= 'A' && c <= 'F' ? c - 'A' + 10
                                                : 16;
            require(d < 16, "invalid unicode escape");
            v = (v << 4) | d;
        }
        return v;
    }
    static void utf8(std::string &out, unsigned u) {
        if (u < 0x80)
            out += char(u);
        else if (u < 0x800) {
            out += char(0xc0 | (u >> 6));
            out += char(0x80 | (u & 63));
        } else if (u < 0x10000) {
            out += char(0xe0 | (u >> 12));
            out += char(0x80 | ((u >> 6) & 63));
            out += char(0x80 | (u & 63));
        } else {
            out += char(0xf0 | (u >> 18));
            out += char(0x80 | ((u >> 12) & 63));
            out += char(0x80 | ((u >> 6) & 63));
            out += char(0x80 | (u & 63));
        }
    }
    std::string string() {
        require(take() == '"', "expected quoted string");
        std::string out;
        while (true) {
            unsigned char c = take();
            if (c == '"')
                return out;
            require(c >= 32, "control character in JSON");
            if (c != '\\') {
                out += char(c);
                continue;
            }
            switch (take()) {
            case '"':
                out += '"';
                break;
            case '\\':
                out += '\\';
                break;
            case '/':
                out += '/';
                break;
            case 'b':
                out += '\b';
                break;
            case 'f':
                out += '\f';
                break;
            case 'n':
                out += '\n';
                break;
            case 'r':
                out += '\r';
                break;
            case 't':
                out += '\t';
                break;
            case 'u': {
                unsigned u = hex4();
                if (u >= 0xd800 && u <= 0xdbff) {
                    require(take() == '\\' && take() == 'u', "missing low surrogate");
                    unsigned low = hex4();
                    require(low >= 0xdc00 && low <= 0xdfff, "invalid low surrogate");
                    u = 0x10000 + ((u - 0xd800) << 10) + (low - 0xdc00);
                } else
                    require(u < 0xdc00 || u > 0xdfff, "unpaired low surrogate");
                utf8(out, u);
                break;
            }
            default:
                throw std::runtime_error("invalid JSON escape");
            }
        }
    }
    Json value(unsigned depth = 0) {
        require(depth < 256, "JSON nesting too deep");
        ws();
        require(pos < s.size(), "empty JSON value");
        char c = s[pos];
        if (c == '"')
            return Json(string());
        if (c == '[') {
            ++pos;
            auto j = Json::array();
            if (consume(']'))
                return j;
            do {
                j.arr.push_back(value(depth + 1));
            } while (consume(','));
            require(consume(']'), "expected ]");
            return j;
        }
        if (c == '{') {
            ++pos;
            auto j = Json::object();
            if (consume('}'))
                return j;
            do {
                ws();
                auto k = string();
                require(consume(':'), "expected :");
                auto v = value(depth + 1);
                require(j.obj.emplace(k, std::move(v)).second, "duplicate JSON key: " + k);
            } while (consume(','));
            require(consume('}'), "expected }");
            return j;
        }
        if (c == 'n' || c == 't' || c == 'f') {
            std::string token = c == 'n' ? "null" : c == 't' ? "true" : "false";
            require(s.compare(pos, token.size(), token) == 0, "invalid JSON literal");
            pos += token.size();
            return c == 'n' ? Json() : Json(c == 't');
        }
        std::size_t start = pos;
        if (c == '-')
            ++pos;
        require(pos < s.size(), "missing number");
        if (s[pos] == '0')
            ++pos;
        else {
            require(s[pos] >= '1' && s[pos] <= '9', "invalid JSON value");
            while (pos < s.size() && s[pos] >= '0' && s[pos] <= '9')
                ++pos;
        }
        if (pos < s.size() && s[pos] == '.') {
            ++pos;
            auto b = pos;
            while (pos < s.size() && s[pos] >= '0' && s[pos] <= '9')
                ++pos;
            require(pos > b, "missing fractional digits");
        }
        if (pos < s.size() && (s[pos] == 'e' || s[pos] == 'E')) {
            ++pos;
            if (pos < s.size() && (s[pos] == '+' || s[pos] == '-'))
                ++pos;
            auto b = pos;
            while (pos < s.size() && s[pos] >= '0' && s[pos] <= '9')
                ++pos;
            require(pos > b, "missing exponent digits");
        }
        Json j;
        j.type = Json::Type::Number;
        j.text = s.substr(start, pos - start);
        return j;
    }
};
std::string quote(const std::string &s) {
    std::string out = "\"";
    const char *digits = "0123456789abcdef";
    for (unsigned char c : s) {
        switch (c) {
        case '"':
            out += "\\\"";
            break;
        case '\\':
            out += "\\\\";
            break;
        case '\n':
            out += "\\n";
            break;
        case '\r':
            out += "\\r";
            break;
        case '\t':
            out += "\\t";
            break;
        default:
            if (c < 32) {
                out += "\\u00";
                out += digits[c >> 4];
                out += digits[c & 15];
            } else
                out += char(c);
        }
    }
    return out + '"';
}
void render(const Json &j, std::ostream &out, int indent, int depth) {
    using T = Json::Type;
    if (j.type == T::Null) {
        out << "null";
        return;
    }
    if (j.type == T::String) {
        out << quote(j.text);
        return;
    }
    if (j.type == T::Number || j.type == T::Boolean) {
        out << j.text;
        return;
    }
    bool object = j.type == T::Object;
    out << (object ? '{' : '[');
    std::size_t n = object ? j.obj.size() : j.arr.size(), i = 0;
    auto entry = [&](const std::string *key, const Json &v) {
        if (i++)
            out << ',';
        if (indent >= 0)
            out << '\n' << std::string((depth + 1) * indent, ' ');
        if (key)
            out << quote(*key) << (indent >= 0 ? ": " : ":");
        render(v, out, indent, depth + 1);
    };
    if (object)
        for (auto &[key, v] : j.obj)
            entry(&key, v);
    else
        for (auto &v : j.arr)
            entry(nullptr, v);
    if (n && indent >= 0)
        out << '\n' << std::string(depth * indent, ' ');
    out << (object ? '}' : ']');
}
} // namespace
Json Json::parse(const std::string &input) {
    Parser parser{input};
    auto j = parser.value();
    parser.ws();
    require(parser.pos == input.size(), "trailing JSON content");
    return j;
}
Json Json::read(const std::string &path) {
    std::ifstream file(path);
    require(bool(file), "cannot read " + path);
    std::ostringstream s;
    s << file.rdbuf();
    return parse(s.str());
}
std::string Json::dump(int indent) const {
    require(indent >= -1 && indent <= 16, "invalid indent");
    std::ostringstream s;
    render(*this, s, indent, 0);
    return s.str();
}
void Json::write(const std::string &path) const {
    std::ofstream file(path);
    require(bool(file), "cannot write " + path);
    file << dump() << '\n';
    require(bool(file), "cannot finish writing " + path);
}
} // namespace ringlab
