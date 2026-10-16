#include "ringlab/experiments.hpp"
#include <iostream>
using namespace ringlab;
template <class F> void rejects(F f, const std::string &name) {
    bool rejected = false;
    try {
        f();
    } catch (const std::exception &) {
        rejected = true;
    }
    require(rejected, "invalid evidence accepted: " + name);
}
int main(int argc, char **argv) {
    try {
        require(argc == 2, "certificate directory required");
        auto data = Json::read(std::string(argv[1]) + "/structured_ring_family.json");
        auto original = Json::read(std::string(argv[1]) + "/structured_ring_example.json");
        auto scaled = data["examples"]["scaled-ten"];
        scaled["field_decomposition"]["stages"][0]["ops"][0][1] = 999999;
        rejects([&] { verify_file(scaled, 0); }, "field decomposition stage");
        scaled = data["examples"]["scaled-eight"];
        scaled["factored_baselines"][0]["circuit"]["ops"][0][1] = 999999;
        rejects([&] { verify_file(scaled, 0, &original); }, "factored direct circuit");
        verify_file(data["examples"]["scaled-eight"], 0, &original);
        auto e = build_case(data["examples"]["scaled-ten"]);
        require(e.direct[0].bits() > 64, "wide matrix regression not exercised");
        auto practical = Json::read(std::string(argv[1]) + "/practical_mds_examples.json");
        auto record = practical["examples"]["systematic-rs-7-4-ring12"];
        auto report = verify_record(record, 4);
        require(report["projection_xors"].integer() == 36 &&
                    record["projection"]["cost"].integer() == 38,
                "live projection reuse");
        require(
            verify_record(practical["examples"]["systematic-rs-17-14-preliminary"], 4)["net_saved"]
                .null(),
            "incomplete baseline became positive claim");
        std::cout << "experiments: field baselines, wide masks, live reuse and incomplete status "
                     "passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
