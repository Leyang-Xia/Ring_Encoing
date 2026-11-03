#include "xd/json.hpp"
#include <iostream>
#include <set>
#include <vector>

using namespace xd;

int main() {
    const auto data = Json::read("data/results/mw_sq_vandermonde_8to17_p0_3_0.json");
    std::vector<std::set<int>> representations;
    for (const auto &mask : data["target_masks"].arr) {
        const auto value = std::stoull(mask.string());
        std::set<int> terms;
        for (int bit = 0; bit < 64; ++bit)
            if ((value >> bit) & 1ULL) terms.insert(bit);
        representations.push_back(std::move(terms));
    }

    int first_unreachable = -1;
    int left_representations = 0;
    int right_representations = 0;
    int step = 0;
    for (const auto &operation : data["raw"]["circuit"]["ops"].arr) {
        const int node = operation[0].int32();
        const int left = operation[1].int32();
        const int right = operation[2].int32();
        bool reachable = false;
        for (auto &terms : representations) {
            if (!terms.count(left) || !terms.count(right)) continue;
            terms.erase(left);
            terms.erase(right);
            terms.insert(node);
            reachable = true;
        }
        if (!reachable) {
            first_unreachable = step;
            for (const auto &terms : representations) {
                left_representations += terms.count(left) != 0;
                right_representations += terms.count(right) != 0;
            }
            break;
        }
        ++step;
    }

    require(first_unreachable == 56, "MW_SQ trace first unreachable step changed");
    require(left_representations == 3 && right_representations == 1,
            "MW_SQ trace operands are no longer split across representations");
    std::cout << "MW_SQ trace replay: first unreachable node 120 at step 56, "
              << "operands split across 3 and 1 representations\n";
}
