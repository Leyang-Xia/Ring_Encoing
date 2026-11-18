#include "ringlab/core.hpp"
namespace ringlab {
Rows benchmark_inputs(int outputs,int inputs,const std::string&pattern,std::uint64_t seed) {
    require(outputs>0&&inputs>0&&std::uint64_t(outputs)*inputs<=5000000,"invalid benchmark dimensions");
    require(pattern=="sparse"||pattern=="dense"||pattern=="blocks","unknown benchmark pattern");
    SeedRandom rng(seed);Rows rows;
    for(int i=0;i<outputs;++i) {
        Mask row;
        if(pattern=="blocks") {
            int boundary=(inputs/2/4)*4;
            for(int a=0;a<boundary;a+=4) if(rng.below(10000)<5000)
                for(int j=a;j<a+4;++j) row.toggle(j);
            for(int j=boundary;j<inputs;++j) if(rng.below(10000)<2000) row.toggle(j);
        } else {
            auto threshold=pattern=="sparse"?1500U:5000U;
            for(int j=0;j<inputs;++j) if(rng.below(10000)<threshold) row.toggle(j);
        }
        rows.push_back(row);
    }
    return rows;
}
} // namespace ringlab
