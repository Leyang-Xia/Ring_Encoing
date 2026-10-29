#pragma once

#include "xd/depth.hpp"
#include <algorithm>
#include <string>
#include <vector>
#include <unordered_map>

namespace xd {

struct Representation {
    std::vector<int> terms;
    ringlab::Mask mask;
    int arrival_depth = 0;
};

struct ParetoFrontier { std::vector<Representation> entries; };

struct Candidate {
    int left = -1;
    int right = -1;
    int arrival_depth = 0;
    int support = 0;
};

struct NodeRecord {
    ringlab::Mask mask;
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
    static SearchState initial(const ringlab::Rows &targets, int inputs,
                               const BoundedScheduleOptions &options);
    std::vector<Candidate> candidates() const;
    SearchState apply(const Candidate &candidate) const;
    bool valid() const { return valid_; }
    bool complete() const;
    int best_cost_lower_bound() const;
    std::string signature() const;
    ringlab::Circuit export_circuit() const;
    const ParetoFrontier &frontier(std::size_t index) const;
    int pair_support(const Candidate &candidate) const;

  private:
    int inputs_ = 0;
    BoundedScheduleOptions options_;
    std::vector<NodeRecord> nodes_;
    std::unordered_map<ringlab::Mask, int, ringlab::MaskHash> node_by_mask_;
    std::vector<ParetoFrontier> frontiers_;
    std::unordered_map<ringlab::Mask, std::vector<PairRef>, ringlab::MaskHash> pair_index_;
    bool valid_ = true;

    void rebuild_pair_index();
    Representation make_representation(const std::vector<int> &terms) const;
    int node_depth(int node) const;
};

int combine_depth(const std::vector<int> &arrival_depths);
std::vector<int> normalize_terms(const std::vector<int> &terms);
bool dominates(const std::vector<int> &a, const std::vector<int> &b);
bool dominates(const Representation &a, const Representation &b);
void pareto_insert(ParetoFrontier &frontier, Representation candidate,
                   std::size_t max_representations = 0);
int candidate_depth(int left_depth, int right_depth);
bool candidate_allowed(int left_depth, int right_depth, int max_depth);
int lookup_node(const std::unordered_map<ringlab::Mask, int, ringlab::MaskHash> &nodes,
                const ringlab::Mask &mask);
void validate_target_width(const ringlab::Rows &targets, int inputs);
std::string canonical_representation(const Representation &representation);

} // namespace xd
