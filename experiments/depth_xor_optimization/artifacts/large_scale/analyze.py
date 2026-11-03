"""Independent verification and reporting of native CLI outputs; no synthesis."""
from pathlib import Path
import json,csv,collections
P=Path(__file__).resolve().parent;R=P/'results_controlled'
inputs={};manifest=[]
for p in sorted((P/'inputs').glob('*.json')):
 r=json.load(open(p));lines=(P/'inputs'/(p.stem+'.txt')).read_text().splitlines();assert len(lines)==r['outputs'] and all(len(x)==r['inputs'] for x in lines)
 rows=[sum((int(x)<<j) for j,x in enumerate(line)) for line in lines];assert rows==[int(x) for x in r['target_masks']]
 D=max(((x.bit_count()-1).bit_length() if x else 0) for x in rows);assert D==r['depth_lower_bound'];inputs[p.stem]=rows
 manifest.append(dict(case=p.stem,outputs=r['outputs'],inputs=r['inputs'],pattern=r['pattern'],seed=r['seed'],density=r['ones']/(r['outputs']*r['inputs']),D0=D,independent_tree_gates=sum(max(0,x.bit_count()-1) for x in rows)))
def audit(c,targets):
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
records=[];taskrecords=[];verified=0
for stem in inputs:
 meta=next(x for x in manifest if x['case']==stem)
 for task in ['frontier','frontier_basic','ours_min','ours_plus1','MW_SQ','UBER_XSET','beam','ours_min30','ours_plus130']:
  base=R/(stem+'_'+task);exitfile=Path(str(base)+'.exit')
  if not exitfile.exists():continue
  code=int(exitfile.read_text());out=Path(str(base)+'.json');r=json.load(open(out)) if out.exists() else {}
  time=Path(str(base)+'.time.txt').read_text().splitlines();seconds=float(time[-1]);status=r.get('status','no_result') if code==0 else ('timeout' if code in (124,137) else 'error')
  budget=30 if task in ('ours_min30','ours_plus130') else 10
  tr=dict(case=stem,task=task,status=status,exit_code=code,elapsed_seconds=seconds,budget_seconds=budget)
  taskrecords.append(tr)
  common={**meta,**tr}
  if status!='verified':
   assert status in ('timeout','search_failed','infeasible_by_lower_bound'),tr
   assert not r.get('circuit'),tr
   records.append(dict(common,algorithm=task,gates=None,depth=None));continue
  assert Path(str(base)+'.verify.exit').read_text().strip()=='0',base
  assert [int(x) for x in r['target_masks']]==inputs[stem]
  entries=r['candidates'] if task.startswith('frontier') else [r]
  for e in entries:
   gates,depth=audit(e['circuit'],inputs[stem]);assert gates==e['cost']
   if 'depth' in e:assert depth==e['depth']
   if 'depth_limit' in e:assert depth<=e['depth_limit']
   verified+=1;records.append(dict(common,algorithm=e.get('label',task),gates=gates,depth=depth))
with open(P/'all_results.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
with open(P/'tasks.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(taskrecords[0]));w.writeheader();w.writerows(taskrecords)
json.dump(manifest,open(P/'input_manifest.json','w'),indent=2)
summary=[]
for m in manifest:
 stem=m['case'];q=[r for r in records if r['case']==stem]
 def get(alg,task=None):
  z=[x for x in q if x['algorithm']==alg and (task is None or x['task']==task)]
  good=[x for x in z if x['status']=='verified']
  return min(good,key=lambda x:(x['gates'],x['depth'])) if good else (z[0] if z else None)
 row=dict(m)
 for tag,alg in [('cse_raw','CSE-chained'),('cse_rebalanced','CSE-chained-rebalanced'),('share2','CSE-share2'),('share3','CSE-share3'),('no_share','CSE-noshare'),('ours_min','ours_min'),('ours_plus1','ours_plus1'),('mw','MW_SQ'),('uber','UBER_XSET'),('beam','beam'),('ours_min30','ours_min30'),('ours_plus130','ours_plus130')]:row[tag]=get(alg)
 summary.append(row)
json.dump(summary,open(P/'summary.json','w'),indent=2)
print('INPUTS',len(inputs),'TASKS',len(taskrecords),'INDEPENDENT_VERIFIED_CIRCUITS',verified,'TASK_STATUS',dict(collections.Counter(r['status'] for r in taskrecords)))
for row in summary:
 print('\n',row['case'],'D0',row['D0'],'density',round(row['density'],3))
 for tag in ['no_share','cse_raw','cse_rebalanced','share3','ours_min','ours_plus1','beam','ours_min30','ours_plus130']:
  r=row[tag]
  if r:print(tag,r['status'],r['gates'],r['depth'],r['elapsed_seconds'])
