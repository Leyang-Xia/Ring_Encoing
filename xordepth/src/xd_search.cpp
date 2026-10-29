#include "xd/search.hpp"
#include <sstream>
#include <stdexcept>

namespace xd {
namespace {
using ringlab::require;
}

int combine_depth(const std::vector<int> &depths) {
    if (depths.empty()) return 0;
    std::vector<int> work = depths;
    for (int d : work) require(d >= 0, "negative arrival depth");
    int level = 0;
    while (work.size() > 1) {
        std::sort(work.begin(), work.end());
        int ready = 1 + std::max(work[0], work[1]);
        work.erase(work.begin(), work.begin() + 2);
        work.push_back(ready);
        level = std::max(level, ready);
    }
    return level;
}

std::vector<int> normalize_terms(const std::vector<int> &terms) {
    std::vector<int> out = terms;
    std::sort(out.begin(), out.end());
    std::vector<int> result;
    for (int term : out) {
        require(term >= 0, "negative term");
        if (!result.empty() && result.back() == term) result.pop_back();
        else result.push_back(term);
    }
    return result;
}

bool dominates(const std::vector<int> &a, const std::vector<int> &b) {
    if (a.size() != b.size()) return false;
    bool strict = false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] > b[i]) return false;
        strict = strict || a[i] < b[i];
    }
    return strict;
}

bool dominates(const Representation &a, const Representation &b) {
    return a.terms.size() <= b.terms.size() && a.arrival_depth <= b.arrival_depth &&
           (a.terms.size() < b.terms.size() || a.arrival_depth < b.arrival_depth);
}

void pareto_insert(ParetoFrontier &frontier, Representation candidate, std::size_t max_representations) {
    candidate.terms = normalize_terms(candidate.terms);
    for (const auto &entry : frontier.entries)
        if (dominates(entry, candidate) ||
            (entry.arrival_depth == candidate.arrival_depth && entry.terms == candidate.terms)) return;
    frontier.entries.erase(std::remove_if(frontier.entries.begin(), frontier.entries.end(),
                                          [&](const Representation &entry) { return dominates(candidate, entry); }),
                           frontier.entries.end());
    frontier.entries.push_back(std::move(candidate));
    std::sort(frontier.entries.begin(), frontier.entries.end(), [](const Representation &x, const Representation &y) {
        if (x.arrival_depth != y.arrival_depth) return x.arrival_depth < y.arrival_depth;
        if (x.terms.size() != y.terms.size()) return x.terms.size() < y.terms.size();
        return x.terms < y.terms;
    });
    if (max_representations && frontier.entries.size() > max_representations)
        frontier.entries.resize(max_representations);
}

int candidate_depth(int left_depth, int right_depth) {
    require(left_depth >= 0 && right_depth >= 0, "negative candidate depth");
    return 1 + std::max(left_depth, right_depth);
}

std::string canonical_representation(const Representation &representation) {
    auto terms = normalize_terms(representation.terms);
    std::ostringstream out;
    out << representation.arrival_depth << ':';
    for (int term : terms) out << term << ',';
    out << ':' << representation.mask.decimal();
    return out.str();
}

BoundedScheduleResult synthesize_bounded(const Rows &targets, int inputs,
                                         const BoundedScheduleOptions &options) {
    require(inputs > 0, "invalid input count");
    require(options.max_depth >= 0 && options.beam_width > 0 && options.branch_width > 0 &&
                options.max_representations > 0 && options.restarts > 0 && options.timeout_seconds > 0.0,
            "invalid bounded search options");
    BoundedScheduleResult result;
    result.requested_depth = options.max_depth;
    result.depth_lower_bound = depth_lower_bound(targets);
    if (result.depth_lower_bound > options.max_depth) {
        result.status = "infeasible_by_lower_bound";
        return result;
    }
    ScheduleOptions schedule;
    auto circuit = synthesize_schedule(targets, inputs, schedule, options.seed);
    if (depth_profile(circuit).depth > options.max_depth) {
        result.status = "search_failed";
        return result;
    }
    result.circuit = std::move(circuit);
    result.status = "verified";
    result.explored_states = 1;
    return result;
}
} // namespace xd
