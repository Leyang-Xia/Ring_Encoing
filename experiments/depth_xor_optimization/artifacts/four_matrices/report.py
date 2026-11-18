"""Read-only independent circuit audit and result reporting; no synthesis."""
from pathlib import Path
import json,csv,hashlib
P=Path(__file__).resolve().parent
records=[]
for name,ds in [('vandermonde',[5,6]),('rs',[5,6]),('cauchy',[5,6]),('gf8',[6,7])]:
 lines=(P/(name+'.txt')).read_text().splitlines();targets=[sum(int(x)<<j for j,x in enumerate(line.strip())) for line in lines]
 for D in ds:
  for method in ['beam','ours']:
   stem=P/f'{method}_{name}_d{D}';exitcode=int(stem.with_suffix('.exit').read_text());r=json.loads(stem.with_suffix('.json').read_text()) if stem.with_suffix('.json').exists() else {}
   status=r['status'] if exitcode==0 else ('external_timeout_190s' if exitcode==124 else ('system_killed_137' if exitcode==137 else 'execution_error'));g=depth=None
   if status=='verified':
    c=r['circuit'];v={j:1<<j for j in range(c['inputs'])};d={j:0 for j in v};parents={}
    if c['zero_node'] is not None:v[c['zero_node']]=0;d[c['zero_node']]=0
    for t,a,b in c['ops']:
     assert t not in v and a in v and b in v
     v[t]=v[a]^v[b];d[t]=1+max(d[a],d[b]);parents[t]=(a,b)
    assert [v[o] for o in c['outputs']]==targets
    live=set();stack=list(c['outputs'])
    while stack:
     t=stack.pop()
     if t in parents and t not in live:live.add(t);stack.extend(parents[t])
    g=len(live);depth=max(d[o] for o in c['outputs']);assert g==len(c['ops'])==c['cost']==r['cost'] and depth==r['depth'] and depth<=D
   else:assert status in ['timeout','search_failed','system_killed_137','execution_error','external_timeout_190s'] and not r.get('circuit')
   records.append(dict(matrix=name,outputs=len(lines),inputs=len(lines[0].strip()),depth_limit=D,method=method,status=status,gates=g,actual_depth=depth,seconds=float(stem.with_suffix('.seconds').read_text()),explored_states=r.get('explored_states')))
with (P/'comparison.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
json.dump(records,open(P/'summary.json','w'),indent=2)
for r in records:print(r)
print('AUDITED',sum(r['status']=='verified' for r in records))
# Independent verification and same-depth best selection from original candidates.
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
combined=[];baseline_count=0
for name,ds in [('vandermonde',[5,6]),('rs',[5,6]),('cauchy',[5,6]),('gf8',[6,7])]:
 lines=(P/(name+'.txt')).read_text().splitlines();targets=[sum(int(x)<<j for j,x in enumerate(line.strip())) for line in lines]
 candidates=json.loads((P/(name+'_native_frontier.json')).read_text())['candidates']
 measured=[]
 for e in candidates:
  g,d=audit_frontier(e['circuit'],targets);assert g==e['cost'] and d==e['depth'];measured.append((g,d,e['label']));baseline_count+=1
 for D in ds:
  b=min(q for q in measured if q[1]<=D)
  beam=next(r for r in records if r['matrix']==name and r['depth_limit']==D and r['method']=='beam')
  ours=next(r for r in records if r['matrix']==name and r['depth_limit']==D and r['method']=='ours')
  combined.append(dict(matrix=name,depth_limit=D,native_best_gates=b[0],native_best_depth=b[1],native_best_label=b[2],beam_status=beam['status'],beam_gates=beam['gates'],beam_depth=beam['actual_depth'],ours_gates=ours['gates'],ours_depth=ours['actual_depth'],saving_vs_native_pct=100*(b[0]-ours['gates'])/b[0] if ours['gates'] is not None else None))
with (P/'combined_comparison.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(combined[0]));w.writeheader();w.writerows(combined)
print('BASELINE_AUDITED',baseline_count)
for r in combined:print(r)
