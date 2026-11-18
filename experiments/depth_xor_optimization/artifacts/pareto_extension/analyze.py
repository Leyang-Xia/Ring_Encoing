"""Independent verification, Pareto selection and plotting; no synthesis."""
from pathlib import Path
import json,csv,hashlib,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
P=Path(__file__).resolve().parent;B=P.parent/'four_matrices';V=P.parent/'v2_results'
if not B.exists():B=P/'reused_evidence'
if not V.exists():V=P/'reused_evidence'
def copy_evidence(a,b):
 if Path(a).resolve()!=Path(b).resolve():shutil.copyfile(a,b)
def audit_frontier(c,targets):
 v={j:1<<j for j in range(c['inputs'])};d={j:0 for j in v};parents={}
 if c['zero_node'] is not None:v[c['zero_node']]=0;d[c['zero_node']]=0
 for t,a,b in c['ops']:
  assert t not in v and a in v and b in v
  v[t]=v[a]^v[b];d[t]=1+max(d[a],d[b]);parents[t]=(a,b)
 assert [v[o] for o in c['outputs']]==targets
 live=set();stack=list(c['outputs'])
 while stack:
  t=stack.pop()
  if t in parents and t not in live:live.add(t);stack.extend(parents[t])
 assert len(live)==len(c['ops'])==c['cost']
 return len(live),max(d[o] for o in c['outputs'])

def pareto(records):
 pts=sorted(set((r['gates'],r['depth']) for r in records))
 return sorted([p for p in pts if not any(q[0]<=p[0] and q[1]<=p[1] and q!=p for q in pts)],key=lambda p:(p[1],p[0]))
allrows=[];summaries=[];verified=0;statuses=[]
cache=P/'reused_evidence';cache.mkdir(exist_ok=True)
for name in ['vandermonde','rs','cauchy','gf8']:
 lines=(B/(name+'.txt')).read_text().splitlines();targets=[sum(int(x)<<j for j,x in enumerate(line.strip())) for line in lines]
 lower=max((x.bit_count()-1).bit_length() if x else 0 for x in targets);cap=14 if name=="rs" else lower+3
 copy_evidence(B/(name+'.txt'),cache/(name+'.txt'))
 native=json.loads((B/(name+'_native_frontier.json')).read_text());copy_evidence(B/(name+'_native_frontier.json'),cache/(name+'_native_frontier.json'))
 records=[]
 for e in native['candidates']:
  g,d=audit_frontier(e['circuit'],targets);assert g==e['cost'] and d==e['depth'];verified+=1
  records.append(dict(matrix=name,method='Repository',label=e['label'],requested_depth=None,gates=g,depth=d,seconds=None,source=str(B/(name+'_native_frontier.json'))))
 # Include the supported depth-aware schedule variant and all previously
 # completed native bounded circuits; failed searches are not circuit points.
 extra_files=[P/f'native_depth_aware_{name}.json']
 extra_files += [q for q in B.glob(f'beam_{name}_d*.json') if '.verify.' not in q.name]
 for extra in extra_files:
  evidence=json.loads(extra.read_text());copy_evidence(extra,cache/extra.name)
  if evidence.get('status') not in [None,'verified']:continue
  entries=evidence['candidates'] if 'candidates' in evidence else [evidence]
  for e in entries:
   g,d=audit_frontier(e['circuit'],targets);assert g==e['cost'] and d==e['depth'];verified+=1
   label=e.get('label',e.get('method','bounded'))
   records.append(dict(matrix=name,method='Repository',label=extra.stem+':'+label,requested_depth=evidence.get('depth_limit'),gates=g,depth=d,seconds=None,source=str(extra)))
 for D in (list(range(lower,lower+4))+([9,11,14] if name=="rs" else [])):
  offset=D-lower
  for method in ['V1','V2']:
   if offset<2:
    stem=V/(f'v1r16_{name}_d{D}' if method=='V1' else f'standalone_{name}_d{D}')
    for ext in ['json','seconds','verify.json']:copy_evidence(stem.with_suffix('.'+ext),cache/stem.with_suffix('.'+ext).name)
   else:stem=P/f'{method.lower()}_{name}_d{D}'
   path=stem.with_suffix('.json');code=int(stem.with_suffix('.exit').read_text()) if offset>=2 else 0
   r=json.loads(path.read_text()) if path.exists() else {}
   status=r.get('status','missing') if code==0 else ('timeout' if code in [124,137] else 'error')
   statuses.append(dict(matrix=name,method=method,depth_limit=D,status=status,exit_code=code))
   if status!='verified':continue
   g,d=audit_frontier(r['circuit'],targets);assert g==r['cost'] and d==r['depth'] and d<=D
   assert r['seed']==1 and r['seeds']==(16 if method=='V1' else 8)
   if method=='V2':assert r['total_synthesis_restarts']==16 and r['initial_source']=='own synthesis'
   assert [int(x) for x in r['target_masks']]==targets
   assert json.loads(stem.with_suffix('.verify.json').read_text())['status']=='verified';verified+=1
   records.append(dict(matrix=name,method=method,label=r['method'],requested_depth=D,gates=g,depth=d,seconds=float(stem.with_suffix('.seconds').read_text()),source=str(path)))
 allrows.extend(records)
 summary=dict(matrix=name,lower=lower,max_tested_depth=cap,frontiers={},full_repository_frontier=pareto([r for r in records if r['method']=='Repository']))
 for method in ['Repository','V1','V2']:
  summary['frontiers'][method]=pareto([r for r in records if r['method']==method and r['depth']<=cap])
 summary['combined_frontier']=pareto([r for r in records if r['depth']<=cap])
 summaries.append(summary)
with (P/'all_candidates.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
(P/'summary.json').write_text(json.dumps(summaries,indent=2))
(P/'task_status.json').write_text(json.dumps(statuses,indent=2))
# Best observed feasible circuit under each ceiling, across all completed runs.
bycap=[]
for s in summaries:
 for D in range(s['lower'],s['max_tested_depth']+1):
  row=dict(matrix=s['matrix'],depth_ceiling=D)
  for method in ['Repository','V1','V2']:
   eligible=[r for r in allrows if r['matrix']==s['matrix'] and r['method']==method and r['depth']<=D]
   best=min(eligible,key=lambda r:(r['gates'],r['depth']));row[method+'_gates']=best['gates'];row[method+'_actual_depth']=best['depth']
  bycap.append(row)
with (P/'best_under_depth.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(bycap[0]));w.writeheader();w.writerows(bycap)
fig,axes=plt.subplots(2,2,figsize=(11,7.2),layout='constrained')
colors={'Repository':'#666666','V1':'#2978b5','V2':'#d44b2a'};markers={'Repository':'s','V1':'o','V2':'^'}
for ax,s in zip(axes.flat,summaries):
 for method in ['Repository','V1','V2']:
  pts=s['frontiers'][method];x=[p[1] for p in pts];y=[p[0] for p in pts]
  ax.step(x+[s['max_tested_depth']],y+[y[-1]],where='post',color=colors[method],alpha=.85)
  ax.scatter(x,y,color=colors[method],marker=markers[method],s=45,label=method,zorder=3)
 ax.set_title({'vandermonde':'Vandermonde (24 x 64)','rs':'Systematic RS (24 x 32)','cauchy':'Cauchy (24 x 32)','gf8':'GF(256) expanded (80 x 80)'}[s['matrix']])
 ax.set_xlabel('Maximum allowed XOR depth');ax.set_ylabel('Best observed XOR gate count');ax.set_xticks(range(s['lower'],s['max_tested_depth']+1));ax.grid(alpha=.22);ax.legend()
fig.suptitle('Observed Pareto fronts: independent V1 / V2 vs repository candidates',fontsize=14)
fig.savefig(P/'pareto_comparison.png',dpi=180);fig.savefig(P/'pareto_comparison.pdf');plt.close(fig)
print('INDEPENDENTLY_VERIFIED',verified,'TASKS',len(statuses))
for s in summaries:print(json.dumps(s))
