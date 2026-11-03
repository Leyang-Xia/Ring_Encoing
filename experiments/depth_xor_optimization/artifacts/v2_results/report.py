"""Independent verification and reporting only; not synthesis."""
from pathlib import Path
import json,csv
P=Path(__file__).resolve().parent
B=P.parent/"four_matrices"
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

rows=[];verified=0
old_table=list(csv.DictReader((B/'combined_comparison.csv').open()))
for name,ds in [('vandermonde',[5,6]),('rs',[5,6]),('cauchy',[5,6]),('gf8',[6,7])]:
 lines=(B/(name+'.txt')).read_text().splitlines();targets=[sum(int(x)<<j for j,x in enumerate(line.strip())) for line in lines]
 for D in ds:
  old=next(x for x in old_table if x['matrix']==name and int(x['depth_limit'])==D)
  result=dict(matrix=name,depth_limit=D,native_best=int(old['native_best_gates']),native_best_depth=int(old['native_best_depth']),beam_status=old['beam_status'],beam_gates=old['beam_gates'],beam_depth=old['beam_depth'],v1_gates=int(old['ours_gates']),v1_depth=int(old['ours_depth']))
  control_path=P/f'v1r16_{name}_d{D}.json'
  control=json.loads(control_path.read_text());cg,cd=audit_frontier(control['circuit'],targets)
  assert cg==control['cost'] and cd==control['depth'] and cd<=D
  assert json.loads(control_path.with_suffix('.verify.json').read_text())['status']=='verified'
  result['v1r16_gates']=cg;result['v1r16_depth']=cd;verified+=1
  for mode in ['standalone','warm']:
   path=P/f'{mode}_{name}_d{D}.json';r=json.loads(path.read_text());assert r['status']=='verified'
   gates,depth=audit_frontier(r['circuit'],targets);assert gates==r['cost'] and depth==r['depth'] and depth<=D
   assert [int(x) for x in r['target_masks']]==targets
   assert r['total_synthesis_restarts']==16
   assert r['original_ours_cost']==int(old['ours_gates'])
   native=json.loads(path.with_suffix('.verify.json').read_text());assert native['status']=='verified'
   result[mode+'_gates']=gates;result[mode+'_depth']=depth;result[mode+'_seconds']=float(path.with_suffix('.seconds').read_text());result[mode+'_source']=r['initial_source'];verified+=1
  assert result['standalone_gates']<=result['v1_gates']
  assert result['warm_gates']<=min(result['standalone_gates'],result['native_best'])
  rows.append(result)
with (P/'comparison.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(P/'summary.json').write_text(json.dumps(rows,indent=2))
print('VERIFIED_NEW_CIRCUITS',verified)
for r in rows:print(r)
