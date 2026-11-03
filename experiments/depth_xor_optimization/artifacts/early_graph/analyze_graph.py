"""Independent exact circuit audit, spectral validation and summary; no synthesis."""
from pathlib import Path
import json,csv,statistics,collections,sys
import numpy as np
P=Path(__file__).resolve().parent;R=P/(sys.argv[1] if len(sys.argv)>1 else 'graph_results');TAG=sys.argv[2] if len(sys.argv)>2 else 'graph'
inputs={(r['n'],r['family'],r['seed']):r['target_masks'] for r in json.load(open(R/'inputs.json'))}
records=[];spec_checked=set();max_eigen_error=0
for p in sorted(R.glob('*.json')):
 if p.name=='inputs.json':continue
 r=json.load(open(p));key=(r['n'],r['family'],r['seed']);assert r['target_masks']==inputs[key]
 c=r['circuit'];v={j:1<<j for j in range(c['inputs'])};d={j:0 for j in v};v[c['zero_node']]=0;d[c['zero_node']]=0;parents={}
 for t,a,b in c['ops']:
  assert t not in v and a in v and b in v
  v[t]=v[a]^v[b];d[t]=1+max(d[a],d[b]);parents[t]=(a,b)
 assert [v[o] for o in c['outputs']]==inputs[key]
 live=set();stack=list(c['outputs'])
 while stack:
  t=stack.pop()
  if t in parents and t not in live:live.add(t);stack.extend(parents[t])
 assert len(live)==len(c['ops'])==r['cost'];assert max(d[o] for o in c['outputs'])==r['depth']<=r['D0']
 f=r['features']
 if key not in spec_checked:
  n=r['n'];A=np.array([[(x>>j)&1 for j in range(n)] for x in inputs[key]],dtype=float)
  W=A.T@A;np.fill_diagonal(W,0);deg=W.sum(axis=1);L=np.diag(deg)-W
  di=np.zeros(n);di[deg>0]=1/np.sqrt(deg[deg>0]);N=np.diag((deg>0).astype(float))-(di[:,None]*W)*di[None,:]
  a=np.linalg.eigvalsh(L);b=np.linalg.eigvalsh(N)
  err=max(float(np.max(np.abs(a-f['ordinary_eigenvalues']))),float(np.max(np.abs(b-f['normalized_eigenvalues']))));max_eigen_error=max(max_eigen_error,err);assert err<1e-6,(key,err)
  spec_checked.add(key)
 records.append(dict(n=r['n'],family=r['family'],seed=r['seed'],method=r['method'],gates=r['cost'],depth=r['depth'],D0=r['D0'],lambda2=f['lambda2'],smallest_positive=f['smallest_positive'],normalized_lambda2=f['normalized_lambda2'],normalized_active_lambda2=f['normalized_active_lambda2'],normalized_smallest_positive=f['normalized_smallest_positive'],components=f['components'],isolated=f['isolated_vertices'],density=f['density'],completed_trials=r['completed_trials'],budget_seconds=r['budget_seconds'],elapsed_seconds=r['elapsed_seconds'],operative_method=r['operative_method'],max_trials=r['max_trials'],comparison_mode=r.get('comparison_mode','equal_time_cap')))
assert len(records)==180 and len(spec_checked)==36
with open(P/(TAG+'_results.csv'),'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
summary=[]
for n in (16,32,64):
 for fam in ('sparse','dense','structured','modular'):
  for method in ('GlobalCSE','RandomGroups','SpectralGroups','CommonBlocks','LambdaGate'):
   rs=[r for r in records if r['n']==n and r['family']==fam and r['method']==method]
   assert len(rs)==3
   summary.append(dict(n=n,family=fam,method=method,gates=statistics.mean(r['gates'] for r in rs),depth=statistics.mean(r['depth'] for r in rs),lambda2=statistics.mean(r['lambda2'] for r in rs),normalized_lambda2=statistics.mean(r['normalized_lambda2'] for r in rs),normalized_active_lambda2=statistics.mean(r['normalized_active_lambda2'] for r in rs),completed_trials=statistics.mean(r['completed_trials'] for r in rs),elapsed_seconds=statistics.mean(r['elapsed_seconds'] for r in rs)))
json.dump(summary,open(P/(TAG+'_summary.json'),'w'),indent=2)
paired=[]
for key in inputs:
 rs={r['method']:r for r in records if (r['n'],r['family'],r['seed'])==key};base=rs['GlobalCSE']['gates']
 paired.append(dict(n=key[0],family=key[1],seed=key[2],lambda2=rs['GlobalCSE']['normalized_active_lambda2'],spectral_gain=(base-rs['SpectralGroups']['gates'])/base,random_gain=(base-rs['RandomGroups']['gates'])/base,block_gain=(base-rs['CommonBlocks']['gates'])/base,gate_gain=(base-rs['LambdaGate']['gates'])/base,gate_choice=rs['LambdaGate']['operative_method']))
json.dump(paired,open(P/(TAG+'_paired.json'),'w'),indent=2)
print('AUDIT',len(records),'circuits;',len(spec_checked),'spectra; maximum eigenvalue error',max_eigen_error)
print('ZERO COMPLETED TRIALS',sum(r['completed_trials']==0 for r in records))
for n in (16,32,64):
 for fam in ('sparse','dense','structured','modular'):
  ss=[r for r in summary if r['n']==n and r['family']==fam]
  print(n,fam,'lambda2norm',round(ss[0]['normalized_active_lambda2'],4),[(r['method'],round(r['gates'],2),round(r['depth'],2),round(r['completed_trials'],1)) for r in ss])
for method in ('RandomGroups','SpectralGroups','CommonBlocks','LambdaGate'):
 wins=ties=losses=0;net=0
 for key in inputs:
  rs={r['method']:r for r in records if (r['n'],r['family'],r['seed'])==key};a=rs['GlobalCSE']['gates'];b=rs[method]['gates'];wins+=b<a;ties+=b==a;losses+=b>a;net+=a-b
 print('vs Global',method,wins,ties,losses,'net saved',net)
x=np.array([r['lambda2'] for r in paired]);y=np.array([r['spectral_gain'] for r in paired]);print('PEARSON lambda2 vs spectral gain',float(np.corrcoef(x,y)[0,1]));print('GATE CHOICES',collections.Counter(r['gate_choice'] for r in paired))
