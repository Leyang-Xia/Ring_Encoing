#pragma once

#include "xd/depth.hpp"
#include <algorithm>
#include <string>
#include <vector>
#include <unordered_map>

namespace xd {

struct Representation {
    std::vector<int> terms;
    BitMask mask;
    int arrival_depth = 0;
    int future_support = 0;
};

struct ParetoFrontier { std::vector<Representation> entries; };

struct Candidate {
    int left = -1;
    int right = -1;
    int arrival_depth = 0;
    int support = 0;
    int target = -1;
    int representation = -1;
    bool toggle = false;
    bool bridge = false;
    int target2 = -1;
    int representation2 = -1;
    int affected_targets = 0;
    int term_reduction = 0;
    int depth_slack = 0;
    int future_support = 0;
    int toggle_penalty = 0;
    // UBER-style completion progress accumulated while scanning pair refs
    // (near-complete count, total remaining terms). Sorting reads these
    // directly instead of re-looking the pair up in a std::map comparator.
    int progress_first = 0;
    int progress_second = 0;
};

struct NodeRecord {
    BitMask mask;
    int depth = 0;
    int left = -1;
    int right = -1;
};

struct PairRef {
    int target = -1;
    int representation = -1;
    int left = -1;
    int right = -1;
};

class SearchState {
  public:
    static SearchState initial(const BitRows &targets, int inputs,
                               const BoundedScheduleOptions &options);
    std::vector<Candidate> candidates() const;
    SearchState apply(const Candidate &candidate) const;
    SearchState apply_toggle(int left, int right, std::size_t target,
                             std::size_t representation, int future_support = 0) const;
    SearchState apply_bridge(int left, int right, std::size_t left_target,
                             std::size_t left_representation, std::size_t right_target,
                             std::size_t right_representation) const;
    std::vector<Candidate> toggle_candidates() const;
    std::vector<Candidate> bridge_candidates() const;
    bool valid() const { return valid_; }
    bool complete() const;
    int best_cost_lower_bound() const;
    // Search priority only: alternatives can overcount mandatory live gates.
    // Never use this score for branch-and-bound or as proof of optimality.
    int cost_priority() const;
    std::string fingerprint() const;
    std::string signature() const;
    Circuit export_circuit() const;
    const ParetoFrontier &frontier(std::size_t index) const;
    int pair_support(const Candidate &candidate) const;

  private:
    int inputs_ = 0;
    BoundedScheduleOptions options_;
    std::vector<NodeRecord> nodes_;
    std::unordered_map<BitMask, int, BitMaskHash> node_by_mask_;
    std::vector<ParetoFrontier> frontiers_;
    // Pair references grouped by edge mask. Every consumer is order
    // insensitive (candidates() ends in total-order sorts, fingerprint()
    // sorts the collected support counts, pair_support/toggle/bridge only
    // count), so the map's iteration order is not observable.
    std::unordered_map<BitMask, std::vector<PairRef>, BitMaskHash> pair_index_;
    std::unordered_map<std::uint64_t, int> residual_pair_support_;
    bool valid_ = true;

    void rebuild_pair_index();
    void rebuild_residual_pair_index();
    void adjust_residual_pair_index(const ParetoFrontier &frontier, int node_count,
                                    int delta);
    // Copy for apply/apply_toggle/apply_bridge: pair_index_ and
    // residual_pair_support_ are derived indexes that every transition fully
    // rebuilds (or the child is invalid and never read), so copying the
    // parent's entries into each child is pure allocation waste.
    SearchState transition_base() const;
    Representation make_representation(const std::vector<int> &terms) const;
    int node_depth(int node) const;
};

int combine_depth(const std::vector<int> &arrival_depths);
std::size_t completion_round_limit(std::size_t outputs, int max_depth);
std::vector<int> normalize_terms(const std::vector<int> &terms);
bool dominates(const std::vector<int> &a, const std::vector<int> &b);
bool dominates(const Representation &a, const Representation &b);
void pareto_insert(ParetoFrontier &frontier, Representation candidate,
                   std::size_t max_representations = 0);
int candidate_depth(int left_depth, int right_depth);
bool candidate_allowed(int left_depth, int right_depth, int max_depth);
int lookup_node(const std::unordered_map<BitMask, int, BitMaskHash> &nodes,
                const BitMask &mask);
void validate_target_width(const BitRows &targets, int inputs);
std::string canonical_representation(const Representation &representation);

} // namespace xd
