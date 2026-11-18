#include "ringlab/core.hpp"
#include "xd/depth.hpp"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
using namespace ringlab;
using Clock=std::chrono::steady_clock;using Pair=std::pair<int,int>;
struct BudgetExpired{};
void check(Clock::time_point end){if(Clock::now()>=end)throw BudgetExpired{};}
struct Builder{
 Circuit c;std::vector<Mask> masks;std::vector<int>d,group;std::unordered_map<Mask,int,MaskHash> cache;
 Builder(int n,const std::vector<int>&g){c.inputs=n;c.zero_node=n;for(int j=0;j<n;++j){masks.push_back(Mask::unit(j));d.push_back(0);group.push_back(g.empty()?-1:g[j]);cache[masks.back()]=j;}masks.emplace_back();d.push_back(0);group.push_back(-1);cache[masks.back()]=n;}
 int lookup(int a,int b){int dt=1+std::max(d[a],d[b]);auto it=cache.find(masks[a]^masks[b]);return it!=cache.end()&&d[it->second]<=dt?it->second:-1;}
 int add(int a,int b){int old=lookup(a,b);if(old>=0)return old;int t=masks.size();auto m=masks[a]^masks[b];int dt=1+std::max(d[a],d[b]);masks.push_back(m);d.push_back(dt);group.push_back(group[a]>=0&&group[a]==group[b]?group[a]:-1);c.ops.push_back({t,a,b});auto it=cache.find(m);if(it==cache.end()||dt<d[it->second])cache[m]=t;return t;}
 int tree(const std::vector<int>&e){std::priority_queue<Pair,std::vector<Pair>,std::greater<Pair>>h;for(int a:e)h.push({d[a],a});while(h.size()>1){int a=h.top().second;h.pop();int b=h.top().second;h.pop();int t=add(a,b);h.push({d[t],t});}return h.empty()?c.zero_node:h.top().second;}
};
std::vector<std::set<int>> expressions(const Rows&r,int n){std::vector<std::set<int>>e;for(auto&m:r){std::set<int>x;for(int j=0;j<n;++j)if(m.bit(j))x.insert(j);e.push_back(x);}return e;}
Circuit balanced(const Rows&r,int n){Builder b(n,{});for(auto&e:expressions(r,n))b.c.outputs.push_back(b.tree(std::vector<int>(e.begin(),e.end())));auto c=b.c.pruned();c.verify(r);return c;}
// Recursive balanced bisection of the normalized co-occurrence Laplacian.
void spectral_split(const std::vector<int>&ids,const std::vector<std::vector<double>>&W,std::vector<int>&g,int&label,Clock::time_point end){
 check(end);int n=ids.size();if(n<=8){for(int id:ids)g[id]=label;++label;return;}
 std::vector<std::vector<double>>a(n,std::vector<double>(n)),v(n,std::vector<double>(n));std::vector<double>deg(n);
 for(int i=0;i<n;++i){v[i][i]=1;for(int j=0;j<n;++j)deg[i]+=W[ids[i]][ids[j]];}
 for(int i=0;i<n;++i)for(int j=0;j<n;++j)a[i][j]=i==j?(deg[i]>0?1:0):(deg[i]>0&&deg[j]>0?-W[ids[i]][ids[j]]/std::sqrt(deg[i]*deg[j]):0);
 for(int sweep=0;sweep<40;++sweep){check(end);double largest=0;for(int p=0;p<n;++p)for(int q=p+1;q<n;++q){double z=a[p][q];largest=std::max(largest,std::abs(z));if(std::abs(z)<1e-11)continue;double tau=(a[q][q]-a[p][p])/(2*z),t=(tau>=0?1.0:-1.0)/(std::abs(tau)+std::sqrt(1+tau*tau));double c=1/std::sqrt(1+t*t),s=t*c;double app=a[p][p],aqq=a[q][q];a[p][p]=app-t*z;a[q][q]=aqq+t*z;a[p][q]=a[q][p]=0;
 for(int k=0;k<n;++k){if(k!=p&&k!=q){double x=a[k][p],y=a[k][q];a[k][p]=a[p][k]=c*x-s*y;a[k][q]=a[q][k]=s*x+c*y;}double x=v[k][p],y=v[k][q];v[k][p]=c*x-s*y;v[k][q]=s*x+c*y;}}
 if(largest<1e-9)break;}
 std::vector<int>eig(n),order(n);std::iota(eig.begin(),eig.end(),0);std::iota(order.begin(),order.end(),0);std::stable_sort(eig.begin(),eig.end(),[&](int i,int j){return a[i][i]<a[j][j];});int f=eig[1];std::stable_sort(order.begin(),order.end(),[&](int i,int j){return v[i][f]<v[j][f];});std::vector<int>left,right;for(int k=0;k<n;++k)(k<n/2?left:right).push_back(ids[order[k]]);spectral_split(left,W,g,label,end);spectral_split(right,W,g,label,end);
}
std::vector<int> spectral_groups(const Rows&r,int n,Clock::time_point end){std::vector<std::vector<double>>W(n,std::vector<double>(n));for(auto&m:r)for(int i=0;i<n;++i)if(m.bit(i))for(int j=i+1;j<n;++j)if(m.bit(j))W[i][j]=++W[j][i];std::vector<int>ids(n),g(n);std::iota(ids.begin(),ids.end(),0);int label=0;spectral_split(ids,W,g,label,end);return g;}
// Diagnostic spectra; ordinary real-valued Laplacians, NOT GF(2) matrices.
std::vector<double> eigvals(std::vector<std::vector<double>>a,Clock::time_point end){int n=a.size();for(int sweep=0;sweep<40;++sweep){check(end);double largest=0;for(int p=0;p<n;++p)for(int q=p+1;q<n;++q){double z=a[p][q];largest=std::max(largest,std::abs(z));if(std::abs(z)<1e-11)continue;double tau=(a[q][q]-a[p][p])/(2*z),t=(tau>=0?1.0:-1.0)/(std::abs(tau)+std::sqrt(1+tau*tau));double c=1/std::sqrt(1+t*t),s=t*c;a[p][p]-=t*z;a[q][q]+=t*z;a[p][q]=a[q][p]=0;for(int k=0;k<n;++k)if(k!=p&&k!=q){double x=a[k][p],y=a[k][q];a[k][p]=a[p][k]=c*x-s*y;a[k][q]=a[q][k]=s*x+c*y;}}if(largest<1e-9)break;}std::vector<double>out;for(int i=0;i<n;++i)out.push_back(a[i][i]);std::sort(out.begin(),out.end());return out;}
Json lap_features(const Rows&r,int n,Clock::time_point end){std::vector<std::vector<double>>W(n,std::vector<double>(n)),L=W,N=W;std::vector<double>deg(n);int nnz=0;for(auto&m:r)for(int i=0;i<n;++i)if(m.bit(i)){++nnz;for(int j=i+1;j<n;++j)if(m.bit(j)){++W[i][j];++W[j][i];}}
 int isolated=0,components=0;std::vector<int>vis(n);for(int i=0;i<n;++i){for(int j=0;j<n;++j)deg[i]+=W[i][j];if(deg[i]==0)++isolated;}
 for(int i=0;i<n;++i)if(!vis[i]){++components;std::vector<int>todo{i};vis[i]=1;while(!todo.empty()){int a=todo.back();todo.pop_back();for(int j=0;j<n;++j)if(W[a][j]>0&&!vis[j]){vis[j]=1;todo.push_back(j);}}}
 for(int i=0;i<n;++i)for(int j=0;j<n;++j){L[i][j]=i==j?deg[i]:-W[i][j];N[i][j]=i==j?(deg[i]>0?1:0):(deg[i]>0&&deg[j]>0?-W[i][j]/std::sqrt(deg[i]*deg[j]):0);}
 auto a=eigvals(L,end),b=eigvals(N,end);int zeros=0;double tol=1e-8*std::max(1.0,a.back()),pos=0,npos=0;for(double x:a){if(std::abs(x)<=tol)++zeros;else if(pos==0&&x>tol)pos=x;}for(double x:b)if(x>1e-8&&npos==0)npos=x;require(zeros==components,"eigensolver zero multiplicity disagrees with components");
 Json av=Json::array(),bv=Json::array();for(double x:a)av.arr.push_back(x);for(double x:b)bv.arr.push_back(x);
 return Json::object({{"lambda2",std::max(0.0,a[1])},{"smallest_positive",pos},{"normalized_lambda2",std::max(0.0,b[1])},{"normalized_smallest_positive",npos},{"normalized_active_lambda2",n-isolated>=2?std::max(0.0,b[isolated+1]):0.0},{"components",components},{"isolated_vertices",isolated},{"active_vertices",n-isolated},{"density",double(nnz)/(r.size()*n)},{"ordinary_eigenvalues",av},{"normalized_eigenvalues",bv}});
}

void block_phase(Builder&b,std::vector<std::set<int>>&e,int D,std::mt19937&rng,Clock::time_point end){
 const int n=b.c.inputs;
 while(true){check(end);std::vector<int>freq(n);for(auto&x:e)for(int a:x)if(a<n)++freq[a];std::set<std::vector<int>>sets;
 for(int i=0;i<int(e.size());++i){check(end);for(int j=i+1;j<int(e.size());++j){std::vector<int>s;std::set_intersection(e[i].begin(),e[i].end(),e[j].begin(),e[j].end(),std::back_inserter(s));s.erase(std::remove_if(s.begin(),s.end(),[&](int a){return a>=n;}),s.end());if(s.size()>8){std::shuffle(s.begin(),s.end(),rng);std::stable_sort(s.begin(),s.end(),[&](int a,int z){return freq[a]>freq[z];});s.resize(8);std::sort(s.begin(),s.end());}if(s.size()>=3)sets.insert(s);}}
 int best=0;std::vector<std::pair<std::vector<int>,std::vector<int>>>tied;
 for(auto&s:sets){check(end);int dt=0;while((1<<dt)<int(s.size()))++dt;std::vector<int>uses;for(int i=0;i<int(e.size());++i){if(!std::includes(e[i].begin(),e[i].end(),s.begin(),s.end()))continue;long long W=0;for(int a:e[i])W+=1LL<<b.d[a];if(W-int(s.size())+(1LL<<dt)<=(1LL<<D))uses.push_back(i);}int gain=(int(uses.size())-1)*(int(s.size())-1);if(gain>best){best=gain;tied.clear();}if(gain==best&&gain>0)tied.push_back({s,uses});}
 if(best<=0)break;auto q=tied[rng()%tied.size()];int t=b.tree(q.first);for(int i:q.second){for(int a:q.first)e[i].erase(a);e[i].insert(t);}
 }
}
Circuit synth(const Rows&r,int n,int D,const std::string&method,const std::vector<int>&fixed,unsigned seed,Clock::time_point end){
 std::mt19937 rng(seed);std::vector<int>groups=fixed;if(method=="RandomGroups"){std::vector<int>order(n);std::iota(order.begin(),order.end(),0);std::shuffle(order.begin(),order.end(),rng);groups.resize(n);for(int i=0;i<n;++i)groups[order[i]]=i/8;}
 Builder b(n,groups);auto e=expressions(r,n);if(method=="CommonBlocks")block_phase(b,e,D,rng,end);
 int phases=(method=="RandomGroups"||method=="SpectralGroups")?2:1;
 for(int phase=0;phase<phases;++phase){bool local=phases==2&&phase==0;
 while(true){check(end);std::map<Pair,std::vector<int>>counts;
 for(int i=0;i<int(e.size());++i){check(end);long long W=0;for(int a:e[i])W+=1LL<<b.d[a];for(auto ia=e[i].begin();ia!=e[i].end();++ia)for(auto ib=std::next(ia);ib!=e[i].end();++ib){int a=*ia,z=*ib;if(local&&(b.group[a]<0||b.group[a]!=b.group[z]))continue;int t=b.lookup(a,z),dt=t>=0?b.d[t]:1+std::max(b.d[a],b.d[z]);if(W-(1LL<<b.d[a])-(1LL<<b.d[z])+(1LL<<dt)<=(1LL<<D))counts[{a,z}].push_back(i);}}
 int best=0,bd=0;long long bw=0;std::vector<Pair>tied;
 for(auto&kv:counts){int a=kv.first.first,z=kv.first.second,t=b.lookup(a,z),dt=t>=0?b.d[t]:1+std::max(b.d[a],b.d[z]);int gain=int(kv.second.size())-(t<0);long long dw=(1LL<<dt)-(1LL<<b.d[a])-(1LL<<b.d[z]);if(gain<=0)continue;if(gain>best||(gain==best&&std::make_pair(dw,dt)<std::make_pair(bw,bd))){best=gain;bw=dw;bd=dt;tied.clear();}if(gain==best&&dw==bw&&dt==bd)tied.push_back(kv.first);}
 if(best<=0)break;Pair p=tied[rng()%tied.size()];int t=b.add(p.first,p.second);for(int i:counts[p]){e[i].erase(p.first);e[i].erase(p.second);if(e[i].count(t))e[i].erase(t);else e[i].insert(t);}
 }}
 for(auto&x:e){check(end);b.c.outputs.push_back(b.tree(std::vector<int>(x.begin(),x.end())));}check(end);auto c=b.c.pruned();c.verify(r);require(xd::depth_profile(c).depth<=D,"graph trial exceeded depth");return c;
}
Rows parse_rows(const Json&item){Rows r;for(auto&v:item["rows_hex"].arr){auto s=v.string();Mask m;int bit=0;for(int k=int(s.size())-1;k>=2;--k){int a=s[k]>='a'?s[k]-'a'+10:s[k]-'0';for(int j=0;j<4;++j)if(a>>j&1)m.toggle(bit+j);bit+=4;}r.push_back(m);}return r;}
int main(int argc,char**argv){try{require(argc==3||argc==4,"usage graph-experiment matrices.json output-dir [fixed3]");std::filesystem::create_directories(argv[2]);Json cases=Json::array();for(auto&i:Json::read(argv[1]).arr)if(i["seed"].int32()<3)cases.arr.push_back(i);
 for(int n:{16,32,64})for(int seed=0;seed<3;++seed){std::mt19937 rng(123456+n*10+seed);std::uniform_real_distribution<double>u(0,1);Rows rows;for(int i=0;i<n;++i){Mask m;int g=i/8;for(int j=0;j<n;++j)if(u(rng)<(j/8==g?.75:.02))m.toggle(j);rows.push_back(m);}cases.arr.push_back(Json::object({{"n",n},{"seed",seed},{"family","modular"},{"target_masks",to_json(rows)}}));}
 Json exact=Json::array();auto start=Clock::now();
 for(auto&item:cases.arr){int n=item["n"].int32(),seed=item["seed"].int32();auto fam=item["family"].string();Rows rows=item.has("target_masks")?rows_json(item["target_masks"]):parse_rows(item);int D=xd::depth_lower_bound(rows);auto features=lap_features(rows,n,Clock::now()+std::chrono::seconds(2));exact.arr.push_back(Json::object({{"n",n},{"seed",seed},{"family",fam},{"target_masks",to_json(rows)}}));
 for(std::string method:{"GlobalCSE","RandomGroups","SpectralGroups","CommonBlocks","LambdaGate"}){
 auto begin=Clock::now();double budget=argc==4?10.0:(n<=32?.15:.8);int max_trials=argc==4?3:128;auto end=begin+std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(budget));Circuit best=balanced(rows,n);int trials=0,aborted=0;std::vector<int>groups;std::string operative=method;
 try{if(method=="LambdaGate"){auto f=lap_features(rows,n,end);operative=f["active_vertices"].int32()>=2&&f["normalized_active_lambda2"].number()<.25?"SpectralGroups":"GlobalCSE";}if(operative=="SpectralGroups")groups=spectral_groups(rows,n,end);for(int k=0;k<max_trials;++k){check(end);auto c=synth(rows,n,D,operative,groups,seed*100+k,end);++trials;if(std::make_pair(c.cost(),xd::depth_profile(c).depth)<std::make_pair(best.cost(),xd::depth_profile(best).depth))best=c;}}
 catch(BudgetExpired&){aborted=1;}
 best.verify(rows);auto dp=xd::depth_profile(best);require(dp.depth<=D,"best graph depth");Json gs=Json::array();for(int g:groups)gs.arr.push_back(g);
 Json r=Json::object({{"n",n},{"seed",seed},{"family",fam},{"method",method},{"status","verified"},{"D0",D},{"cost",best.cost()},{"depth",dp.depth},{"features",features},{"operative_method",operative},{"lambda_gate_threshold",.25},{"completed_trials",trials},{"aborted_trial",aborted},{"budget_seconds",budget},{"elapsed_seconds",std::chrono::duration<double>(Clock::now()-begin).count()},{"group_size",8},{"max_trials",max_trials},{"comparison_mode",argc==4?"fixed3":"equal_time_cap"},{"spectral_groups",gs},{"target_masks",to_json(rows)},{"circuit",best.json()}});r.write(std::string(argv[2])+"/"+std::to_string(n)+"_"+fam+"_"+std::to_string(seed)+"_"+method+".json");}
 std::cout<<n<<" "<<fam<<" "<<seed<<" done "<<std::chrono::duration<double>(Clock::now()-start).count()<<"s"<<std::endl;
 }
 exact.write(std::string(argv[2])+"/inputs.json");
 }catch(std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
