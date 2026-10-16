#include "ringlab/experiments.hpp"
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <regex>
#include <set>
#include <sstream>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>

#ifndef RINGLAB_XSETS_SHA256
#define RINGLAB_XSETS_SHA256 "2a05bf43ff8909b5dabea0990f522003fb445f018261bfcfbf3a2a5ffce9afde"
#endif
namespace ringlab {
std::string backend_source_hash() {
    return RINGLAB_XSETS_SHA256;
}
namespace {
volatile std::sig_atomic_t cancellation_signal = 0;
volatile std::sig_atomic_t backend_pid = 0;
void cancel_backend(int signal) {
    const int saved_errno = errno;
    cancellation_signal = signal;
    const auto pid = static_cast<pid_t>(backend_pid);
    if (pid > 0) {
        kill(-pid, SIGKILL);
        kill(pid, SIGKILL);
    }
    errno = saved_errno;
}
// execute() runs synchronously. Reap cancelled children before restoring handlers.
struct CancellationGuard {
    struct sigaction old_int{}, old_term{};
    CancellationGuard() {
        cancellation_signal = 0;
        backend_pid = 0;
        struct sigaction action{};
        action.sa_handler = cancel_backend;
        sigemptyset(&action.sa_mask);
        require(sigaction(SIGINT, &action, &old_int) == 0, "cannot install SIGINT handler");
        if (sigaction(SIGTERM, &action, &old_term) != 0) {
            sigaction(SIGINT, &old_int, nullptr);
            throw std::runtime_error("cannot install SIGTERM handler");
        }
    }
    ~CancellationGuard() {
        backend_pid = 0;
        sigaction(SIGINT, &old_int, nullptr);
        sigaction(SIGTERM, &old_term, nullptr);
    }
};
struct Temporary {
    std::string path;
    int fd = -1;
    Temporary() {
        char name[] = "/tmp/ringlab-XXXXXX";
        fd = mkstemp(name);
        require(fd >= 0, "cannot create backend temporary file");
        path = name;
    }
    ~Temporary() {
        if (fd >= 0)
            close(fd);
        if (!path.empty())
            unlink(path.c_str());
    }
    Temporary(const Temporary &) = delete;
    Temporary &operator=(const Temporary &) = delete;
    std::string contents() const {
        std::ifstream f(path, std::ios::binary);
        std::ostringstream s;
        s << f.rdbuf();
        return s.str();
    }
};
Circuit parse_schedule(const std::string &text, const Rows &rows, const Rows &tool_rows,
                       int inputs) {
    Circuit c;
    c.inputs = inputs;
    c.zero_node = -1;
    Values values;
    for (int i = 0; i < inputs; ++i)
        values.emplace(i, Mask::unit(i));
    values.emplace(-1, Mask());
    std::map<int, int> printed;
    std::regex gate("s([0-9]+) = s([0-9]+) \\+ s([0-9]+)"), output("b([0-9]+) = s([0-9]+)"),
        count("Total XOR operations: ([0-9]+)");
    std::istringstream stream(text);
    std::string line;
    std::smatch match;
    int reported = -1;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (std::regex_match(line, match, gate)) {
            int node = std::stoi(match[1]), a = std::stoi(match[2]), b = std::stoi(match[3]);
            require(!values.count(node) && values.count(a) && values.count(b),
                    "backend emitted invalid gate");
            values.emplace(node, values.at(a) ^ values.at(b));
            c.ops.push_back({node, a, b});
        } else if (std::regex_match(line, match, output)) {
            require(printed.emplace(std::stoi(match[1]), std::stoi(match[2])).second,
                    "backend duplicate output");
        } else if (std::regex_search(line, match, count))
            reported = std::stoi(match[1]);
    }
    require(reported == c.cost(), "backend reported XOR count does not match schedule");
    require(printed.size() == tool_rows.size(), "backend output count");
    for (std::size_t i = 0; i < tool_rows.size(); ++i)
        require(values.at(printed.at(i)) == tool_rows[i],
                "backend output differs from requested matrix");
    std::unordered_map<Mask, int, MaskHash> bymask;
    for (int i = 0; i < inputs; ++i)
        bymask[Mask::unit(i)] = i;
    bymask[Mask()] = -1;
    for (auto &op : c.ops)
        bymask[values.at(op[0])] = op[0];
    for (auto &row : rows)
        c.outputs.push_back(bymask.at(row));
    c = c.pruned();
    c.verify(rows);
    return c;
}
} // namespace
ProcessResult execute(const std::vector<std::string> &args, const std::string &input,
                      double timeout) {
    require(!args.empty() && timeout > 0 && std::isfinite(timeout),
            "positive finite timeout and command required");
    Temporary in, out, err;
    std::size_t written = 0;
    while (written < input.size()) {
        auto n = write(in.fd, input.data() + written, input.size() - written);
        if (n < 0 && errno == EINTR)
            continue;
        require(n > 0, "cannot write backend input");
        written += n;
    }
    require(lseek(in.fd, 0, SEEK_SET) >= 0, "cannot rewind backend input");
    std::vector<char *> argv;
    for (auto &arg : args)
        argv.push_back(const_cast<char *>(arg.c_str()));
    argv.push_back(nullptr);
    auto start = std::chrono::steady_clock::now();
    CancellationGuard cancellation;
    pid_t child = fork();
    require(child >= 0, "cannot fork backend");
    if (child == 0) {
        signal(SIGINT, SIG_DFL);
        signal(SIGTERM, SIG_DFL);
        setpgid(0, 0);
        if (dup2(in.fd, STDIN_FILENO) < 0 || dup2(out.fd, STDOUT_FILENO) < 0 ||
            dup2(err.fd, STDERR_FILENO) < 0)
            _exit(126);
        close(in.fd);
        close(out.fd);
        close(err.fd);
        execvp(argv[0], argv.data());
        const char message[] = "could not execute backend\n";
        write(STDERR_FILENO, message, sizeof(message) - 1);
        _exit(127);
    }
    backend_pid = child;
    setpgid(child, child);
    ProcessResult result;
    int status = 0;
    while (true) {
        if (cancellation_signal) {
            kill(-child, SIGKILL);
            kill(child, SIGKILL);
        }
        auto done = waitpid(child, &status, WNOHANG);
        if (done == child)
            break;
        if (done < 0 && errno == EINTR)
            continue;
        if (done < 0) {
            kill(-child, SIGKILL);
            kill(child, SIGKILL);
            while (waitpid(child, &status, 0) < 0 && errno == EINTR) {
            }
            throw std::runtime_error("cannot wait for backend");
        }
        if (std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() >=
            timeout) {
            result.timeout = true;
            kill(-child, SIGKILL);
            kill(child, SIGKILL);
            while (waitpid(child, &status, 0) < 0 && errno == EINTR) {
            }
            break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    result.seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    result.exit_code = WIFEXITED(status)     ? WEXITSTATUS(status)
                       : WIFSIGNALED(status) ? 128 + WTERMSIG(status)
                                             : -1;
    result.cancel_signal = cancellation_signal;
    result.output = out.contents();
    result.error = err.contents();
    return result;
}
Json run_backend(const std::string &binary, const Rows &rows, int inputs, const std::string &method,
                 const std::array<int, 3> &params, double timeout, const std::string &path) {
    require(method == "MW_SQ" || method == "UBER_XSET", "unknown X-Sets method");
    require(inputs > 0 && timeout > 0 && std::isfinite(timeout), "invalid input count or timeout");
    for (auto v : params)
        require(v >= 0, "backend parameters must be nonnegative");
    auto record = Json::object(
        {{"path", path},
         {"method", method},
         {"params", Json::array({params[0], params[1], params[2]})},
         {"timeout", timeout},
         {"xsets_source_sha256", backend_source_hash()},
         {"backend_executable", binary},
         {"source_hash_scope", "bundled X-Sets source; a custom --backend is not authenticated"},
         {"cost", Json()}});
    auto start = std::chrono::steady_clock::now();
    try {
        std::set<Mask> normalized;
        for (auto &row : rows) {
            require(row.bits() <= std::size_t(inputs), "backend target exceeds inputs");
            if (row.count() > 1)
                normalized.insert(row);
        }
        Rows tool_rows(normalized.begin(), normalized.end());
        Circuit c;
        if (tool_rows.empty()) {
            c.inputs = inputs;
            c.zero_node = -1;
            for (auto &row : rows)
                c.outputs.push_back(row.empty() ? -1 : int(row.bits() - 1));
            c.verify(rows);
        } else {
            std::string text;
            for (auto &row : tool_rows) {
                for (int bit = 0; bit < inputs; ++bit)
                    text += row.bit(bit) ? '1' : '0';
                text += '\n';
            }
            auto process = execute({binary, std::to_string(params[0]), std::to_string(params[1]),
                                    std::to_string(params[2]), method},
                                   text, timeout);
            record["exit_code"] = process.exit_code;
            if (process.cancel_signal) {
                record["status"] = "cancelled";
                record["cancel_signal"] = process.cancel_signal;
                record["wall_seconds"] = process.seconds;
                return record;
            }
            if (process.timeout) {
                record["status"] = "timeout";
                record["wall_seconds"] = process.seconds;
                return record;
            }
            require(process.exit_code == 0,
                    "backend exit " + std::to_string(process.exit_code) + ": " + process.error);
            c = parse_schedule(process.output, rows, tool_rows, inputs);
        }
        record["status"] = "verified";
        record["cost"] = c.cost();
        record["circuit"] = c.json();
    } catch (const std::exception &e) {
        record["status"] = "error";
        record["error"] = e.what();
    }
    record["wall_seconds"] =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return record;
}
} // namespace ringlab
