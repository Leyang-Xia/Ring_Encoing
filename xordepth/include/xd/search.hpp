#pragma once

#include "xd/depth.hpp"
#include <algorithm>
#include <string>
#include <vector>

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

int combine_depth(const std::vector<int> &arrival_depths);
std::vector<int> normalize_terms(const std::vector<int> &terms);
bool dominates(const std::vector<int> &a, const std::vector<int> &b);
bool dominates(const Representation &a, const Representation &b);
void pareto_insert(ParetoFrontier &frontier, Representation candidate,
                   std::size_t max_representations = 0);
int candidate_depth(int left_depth, int right_depth);
std::string canonical_representation(const Representation &representation);

} // namespace xd
