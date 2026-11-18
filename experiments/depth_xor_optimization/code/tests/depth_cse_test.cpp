#include "ringlab/core.hpp"
#include <iostream>
#include <random>
using namespace ringlab;
int main(){try{
 // Both consumers must be rewritten together: each rewrite alone saves no gate.
 Circuit shared;shared.inputs=4;shared.zero_node=4;
 shared.ops={{{5,0,1}},{{6,0,2}},{{7,0,3}},{{8,5,2}},{{9,5,3}}};
 shared.outputs={6,7,8,9};Rows multi{Mask(5),Mask(9),Mask(7),Mask(11)};
 shared.verify(multi);auto reduced=optimize_depth_gates(shared,multi,2);
 require(reduced.cost()==4&&circuit_depth(reduced)==2,"transactional multi-consumer gate elimination");
 require(shared.cost()==5,"optimizer preserves caller circuit");
 bool invalid_initial=false;try{optimize_depth_gates(shared,multi,1);}catch(const std::exception&){invalid_initial=true;}
 require(invalid_initial,"reject initial depth violation");
 auto no_budget=optimize_depth_gates(shared,multi,2,0);require(no_budget.cost()==5,"zero rewrite budget");
 Circuit slack;slack.inputs=5;slack.zero_node=5;
 slack.ops={{{6,0,1}},{{7,0,4}},{{8,7,2}},{{9,7,3}},{{10,1,4}},{{11,6,2}},{{12,6,3}}};
 slack.outputs={8,9,10,11,12};Rows slack_rows{Mask(21),Mask(25),Mask(18),Mask(7),Mask(11)};
 auto slack_reduced=optimize_depth_gates(slack,slack_rows,3);
 require(slack_reduced.cost()<slack.cost()&&circuit_depth(slack_reduced)<=3,"consume timing slack to save gates");
 Rows example{Mask(15),Mask(51),Mask(60)};DepthCseOptions o;o.max_depth=2;
 auto c=synthesize_depth_cse(example,6,o);c.verify(example);require(c.cost()==6&&circuit_depth(c)==2,"shared-pair example");
 o.lookahead=true;o.restarts=3;o.seed=7;auto a=synthesize_depth_cse(example,6,o),b=synthesize_depth_cse(example,6,o);require(a.same(b),"fixed options deterministic");
 Rows aliases{Mask(),Mask::unit(0),Mask::unit(0),Mask::unit(199)};o=DepthCseOptions{};auto alias=synthesize_depth_cse(aliases,200,o);require(alias.cost()==0&&circuit_depth(alias)==0,"zero/identity/repeated aliases");
 Mask wide=Mask::unit(0)^Mask::unit(199);auto w=synthesize_depth_cse({wide,wide},200,o);require(w.cost()==1&&circuit_depth(w)==1,"wide masks and duplicate outputs");
 bool rejected=false;o.max_depth=0;try{synthesize_depth_cse({wide},200,o);}catch(const std::exception&){rejected=true;}require(rejected,"reject impossible depth");
 std::mt19937 rng(123);for(int k=0;k<20;++k){Rows rows;for(int i=0;i<9;++i){Mask m;for(int j=0;j<10;++j)if(rng()%2)m.toggle(j);rows.push_back(m);}for(bool look:{false,true})for(bool refine:{false,true}){o=DepthCseOptions{};o.lookahead=look;o.refine=refine;auto r=synthesize_depth_cse(rows,10,o);r.verify(rows);auto v2=optimize_depth_gates(r,rows,linear_depth_lower_bound(rows),40);v2.verify(rows);require(v2.cost()<=r.cost(),"rewrite gate count monotonicity");require(circuit_depth(v2)<=linear_depth_lower_bound(rows),"rewrite hard limit");require(circuit_depth(r)<=linear_depth_lower_bound(rows),"random case hard depth");}}
 auto generated=benchmark_inputs(300,200,"blocks",123);require(generated.size()==300,"benchmark row count");
 for(auto &row:generated)require(row.bits()<=200,"benchmark column count");
 require(generated==benchmark_inputs(300,200,"blocks",123),"benchmark seed repeatability");
 bool bad_pattern=false;try{benchmark_inputs(10,10,"invalid",0);}catch(const std::exception&){bad_pattern=true;}require(bad_pattern,"benchmark invalid pattern");
 std::cout<<"depth-CSE: sharing, reproducibility, wide masks, aliases, hard limits and refinements passed\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
