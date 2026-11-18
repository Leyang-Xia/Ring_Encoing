#include "ringlab/core.hpp"
#include <algorithm>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <tuple>

namespace ringlab {
namespace {
std::unordered_map<int,int> node_depths(const Circuit &c) {
    std::unordered_map<int,int> depths;
    for (int i=0;i<c.inputs;++i) depths[i]=0;
    if(c.has_zero) depths[c.zero_node]=0;
    for(auto &op:c.ops) {
        require(!depths.count(op[0]) && depths.count(op[1]) && depths.count(op[2]), "invalid XOR dependency");
        depths[op[0]]=1+std::max(depths.at(op[1]),depths.at(op[2]));
    }
    return depths;
}
} // namespace
int circuit_depth(const Circuit &c) {
    auto depths=node_depths(c);int result=0;
    for(int output:c.outputs) result=std::max(result,depths.at(output));
    return result;
}
int linear_depth_lower_bound(const Rows &rows) {
    int result=0;
    for(auto &row:rows) {int d=0;std::size_t cap=1;while(cap<std::size_t(row.count())){cap*=2;++d;}result=std::max(result,d);}
    return result;
}
namespace {
using Pair=std::pair<int,int>;
struct Candidate {int a,b,t,dt,gain,count;long long delta;std::vector<int> uses;int score=0;};
Circuit greedy_depth_cse(const Rows&rows,int n,int D,int variant,unsigned seed){
 std::mt19937 rng(seed);Circuit c;c.inputs=n;c.zero_node=n;std::vector<Mask> masks;std::vector<int>d(n+1,0);
 std::unordered_map<Mask,int,MaskHash> cache;
 for(int j=0;j<n;++j){masks.push_back(Mask::unit(j));cache[masks.back()]=j;}masks.emplace_back();cache[masks.back()]=n;
 std::vector<std::set<int>> expr;for(auto&r:rows){std::set<int> e;for(int j=0;j<n;++j)if(r.bit(j))e.insert(j);expr.push_back(e);}
 auto getnode=[&](int a,int b){Mask m=masks[a]^masks[b];int depth=1+std::max(d[a],d[b]);auto it=cache.find(m);if(variant&&it!=cache.end()&&d[it->second]<=depth)return it->second;
  int t=masks.size();masks.push_back(m);d.push_back(depth);c.ops.push_back({t,a,b});auto old=cache.find(m);if(old==cache.end()||depth<d[old->second])cache[m]=t;return t;};
 // Preserve a feasible tree completion for every output after each replacement.
 while(true){
  std::map<Pair,std::vector<int>> counts;std::vector<long long> W(expr.size(),0);
  for(size_t i=0;i<expr.size();++i){for(int a:expr[i])W[i]+=1LL<<d[a];for(auto ia=expr[i].begin();ia!=expr[i].end();++ia)for(auto ib=std::next(ia);ib!=expr[i].end();++ib){int a=*ia,b=*ib,dt=1+std::max(d[a],d[b]);auto it=cache.find(masks[a]^masks[b]);if(variant&&it!=cache.end()&&d[it->second]<=dt)dt=d[it->second];if(W[i]-(1LL<<d[a])-(1LL<<d[b])+(1LL<<dt)<=(1LL<<D))counts[{a,b}].push_back(i);}}
  std::vector<Candidate> cand;int bestgain=0;
  for(auto &kv:counts){auto [a,b]=kv.first;int dt=1+std::max(d[a],d[b]),t=-1;auto it=cache.find(masks[a]^masks[b]);if(variant&&it!=cache.end()&&d[it->second]<=dt){t=it->second;dt=d[t];}int gain=int(kv.second.size())-(t<0);if(gain<=0)continue;bestgain=std::max(bestgain,gain);Candidate q{a,b,t,dt,gain,int(kv.second.size()),(1LL<<dt)-(1LL<<d[a])-(1LL<<d[b]),kv.second};q.score=gain;cand.push_back(q);}
  if(cand.empty())break;
  auto immediate=[](const Candidate&a,const Candidate&b){return std::make_tuple(-a.gain,a.delta,a.dt,a.a,a.b)<std::make_tuple(-b.gain,b.delta,b.dt,b.a,b.b);};
  std::sort(cand.begin(),cand.end(),immediate);
  if(variant==2){
   // Examine at most eight strong pairs, including equal-score candidates in random order.
   std::shuffle(cand.begin(),cand.end(),rng);std::stable_sort(cand.begin(),cand.end(),[&](const Candidate&a,const Candidate&b){return std::make_tuple(-a.gain,a.delta,a.dt)<std::make_tuple(-b.gain,b.delta,b.dt);});
   cand.erase(std::remove_if(cand.begin(),cand.end(),[&](const Candidate&q){return q.gain<bestgain-1;}),cand.end());if(cand.size()>8)cand.resize(8);
   for(auto&q:cand){std::map<int,int> future;
    for(int i:q.uses){long long wp=W[i]+q.delta;for(int z:expr[i])if(z!=q.a&&z!=q.b&&z!=q.t){int dt2=1+std::max(q.dt,d[z]);if(wp-(1LL<<q.dt)-(1LL<<d[z])+(1LL<<dt2)<=(1LL<<D))++future[z];}}
    int h=0;for(auto &kv:future)h=std::max(h,kv.second);q.score=q.gain+std::max(0,h-1);
   }
  }
  auto key=[](const Candidate&q){return std::make_tuple(q.score,q.gain,-q.delta,-q.dt);};
  auto best=key(cand.front());for(auto&q:cand)best=std::max(best,key(q));std::vector<int> tied;for(int i=0;i<int(cand.size());++i)if(key(cand[i])==best)tied.push_back(i);
  auto q=cand[tied[rng()%tied.size()]];int t=getnode(q.a,q.b);
  for(int i:q.uses){expr[i].erase(q.a);expr[i].erase(q.b);if(expr[i].count(t))expr[i].erase(t);else expr[i].insert(t);}
 }
 for(auto&e:expr){std::priority_queue<Pair,std::vector<Pair>,std::greater<Pair>>h;for(int a:e)h.push({d[a],a});while(h.size()>1){int a=h.top().second;h.pop();int b=h.top().second;h.pop();int t=getnode(a,b);h.push({d[t],t});}c.outputs.push_back(h.empty()?n:h.top().second);}
 auto p=c.pruned();p.verify(rows);require(circuit_depth(p)<=D,"depth constraint violated");return p;
}
std::pair<int,int> metric(const Circuit&c){return {c.cost(),circuit_depth(c)};}
// Reconnect a gate only to earlier equivalent signals; accept a strict
// (live gates, output depth) improvement after full-circuit verification.
Circuit improve(Circuit c,const Rows&rows,int D,int&attempts,int&accepted){
 attempts=accepted=0;
 for(int pass=0;pass<12&&attempts<200;++pass){bool changed=false;auto vals=c.verify(rows);auto dep=node_depths(c);
 std::vector<int> available;std::unordered_map<Mask,int,MaskHash> bymask;
 auto add=[&](int id){available.push_back(id);auto it=bymask.find(vals.at(id));if(it==bymask.end()||dep.at(id)<dep.at(it->second))bymask[vals.at(id)]=id;};
 for(int i=0;i<c.inputs;++i) { add(i); }
 add(c.zero_node);
 for(int k=0;k<int(c.ops.size())&&!changed&&attempts<200;++k){auto op=c.ops[k];int t=op[0];std::vector<std::array<int,3>> candidates;
 for(int a:available){auto it=bymask.find(vals.at(a)^vals.at(t));if(it==bymask.end())continue;int b=it->second;if(a>b)continue;if((a==op[1]&&b==op[2])||(a==op[2]&&b==op[1]))continue;int d=1+std::max(dep.at(a),dep.at(b));if(d>D)continue;candidates.push_back({d,a,b});}
 std::sort(candidates.begin(),candidates.end());if(candidates.size()>4)candidates.resize(4);
 for(auto &q:candidates){++attempts;Circuit trial=c;trial.ops[k]={t,q[1],q[2]};trial=trial.pruned();auto mt=metric(trial);if(mt.second<=D&&mt<metric(c)){trial.verify(rows);c=trial;++accepted;changed=true;break;}}
 if(!changed)add(t);
 }
 if(!changed)break;
 }
 c.verify(rows);require(circuit_depth(c)<=D,"resubstitution depth");return c;
}

} // namespace
Circuit synthesize_depth_cse(const Rows &rows,int inputs,const DepthCseOptions &options) {
    require(inputs>0 && inputs<1000000,"invalid depth-CSE input count");
    require(!rows.empty(),"depth-CSE requires output rows");
    for(auto &row:rows) require(row.bits()<=std::size_t(inputs),"target exceeds input width");
    require(options.restarts>0 && options.restarts<=1000000,"depth-CSE restarts must be positive");
    int lower=linear_depth_lower_bound(rows);
    int D=options.max_depth<0?lower:options.max_depth;
    require(D>=lower && D<=60,"depth-CSE limit must cover row-weight bound and be <=60");
    Circuit best;bool found=false;
    for(int k=0;k<options.restarts;++k) {
        auto c=greedy_depth_cse(rows,inputs,D,options.lookahead?2:1,unsigned(options.seed+k));
        if(!found || metric(c)<metric(best)) {best=c;found=true;}
    }
    if(options.refine) {int attempts=0,accepted=0;best=improve(best,rows,D,attempts,accepted);}
    best.verify(rows);
    require(circuit_depth(best)<=D,"depth-CSE final depth limit");
    return best;
}
} // namespace ringlab

namespace ringlab {
// Remove a shared gate only when every consumer can be reconstructed from
// earlier signals. Unlike single-gate hill climbing, all consumer rewrites
// are committed together, allowing individually neutral rewrites to save a gate.
Circuit optimize_depth_gates(Circuit circuit, const Rows &targets, int limit,
                             int max_trials) {
    circuit = circuit.pruned();
    circuit.verify(targets);
    require(circuit_depth(circuit) <= limit, "initial circuit exceeds depth limit");
    require(max_trials >= 0, "negative rewrite budget");
    int trials = 0;
    for (int pass = 0; pass < 4 && trials < max_trials; ++pass) {
        bool changed = false;
        const auto values = circuit.verify(targets);
        // Latest arrival times from output deadlines allow consuming slack
        // away from the critical path while preserving the global hard bound.
        std::unordered_map<int,int> latest;
        for(int output:circuit.outputs) latest[output]=limit;
        for(auto it=circuit.ops.rbegin();it!=circuit.ops.rend();++it) {
            int deadline=latest.at((*it)[0])-1;
            for(int parent:{(*it)[1],(*it)[2]}) {
                auto old=latest.find(parent);
                if(old==latest.end() || deadline<old->second) latest[parent]=deadline;
            }
        }
        std::unordered_map<int,int> uses;
        std::set<int> outputs(circuit.outputs.begin(), circuit.outputs.end());
        for (auto op : circuit.ops) { ++uses[op[1]]; ++uses[op[2]]; }
        // Low fanout first; attempt every consumer as one transaction.
        std::vector<int> victims;
        for (auto op : circuit.ops)
            if (!outputs.count(op[0])) victims.push_back(op[0]);
        std::stable_sort(victims.begin(), victims.end(), [&](int a,int b) {
            return uses[a] < uses[b];
        });
        for (int victim : victims) {
            if (trials++ >= max_trials) break;
            Circuit next; next.inputs=circuit.inputs;
            next.has_zero=circuit.has_zero; next.zero_node=circuit.zero_node;
            std::unordered_map<int,int> alias, ready;
            std::unordered_map<Mask,int,MaskHash> bymask;
            std::vector<int> available;
            auto add = [&](int id, int depth) {
                alias[id]=id; ready[id]=depth; available.push_back(id);
                auto it=bymask.find(values.at(id));
                if (it==bymask.end() || depth<ready.at(it->second))
                    bymask[values.at(id)]=id;
            };
            for (int i=0;i<circuit.inputs;++i) add(i,0);
            if (circuit.has_zero) add(circuit.zero_node,0);
            bool valid=true;
            for (auto op : circuit.ops) {
                int t=op[0]; if(t==victim) continue;
                const Mask &target=values.at(t);
                auto equal=bymask.find(target);
                if(equal!=bymask.end() && ready.at(equal->second)<=latest.at(t)) {
                    alias[t]=equal->second; continue;
                }
                int left=-1,right=-1,best_depth=latest.at(t)+1;
                if(alias.count(op[1]) && alias.count(op[2])) {
                    left=alias.at(op[1]);right=alias.at(op[2]);
                    best_depth=1+std::max(ready.at(left),ready.at(right));
                } else {
                    for(int a:available) {
                        auto it=bymask.find(values.at(a)^target);
                        if(it==bymask.end()) continue;
                        int b=it->second, depth=1+std::max(ready.at(a),ready.at(b));
                        if(depth<best_depth) {left=a;right=b;best_depth=depth;}
                    }
                }
                if(left<0 || best_depth>latest.at(t)) {valid=false;break;}
                next.ops.push_back({t,left,right}); add(t,best_depth);
            }
            if(!valid) continue;
            for(int output:circuit.outputs) {
                if(!alias.count(output)) {valid=false;break;}
                next.outputs.push_back(alias.at(output));
            }
            if(!valid) continue;
            next=next.pruned(); next.verify(targets);
            require(circuit_depth(next)<=limit,"rewrite exceeds depth limit");
            if(next.cost()<circuit.cost()) {
                circuit=std::move(next);changed=true;break;
            }
        }
        if(!changed) break;
    }
    circuit.verify(targets);
    require(circuit_depth(circuit)<=limit,"final rewrite depth limit");
    return circuit;
}
} // namespace ringlab
