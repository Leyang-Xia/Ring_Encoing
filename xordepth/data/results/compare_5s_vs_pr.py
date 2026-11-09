import json, glob, os

pr_dir = 'data/results/greedy_portfolio_full_20261104'
fin_dir = 'data/results/xordepth_full_incumbent_20261106'

pr = {}
for f in glob.glob(pr_dir + '/*.json'):
    base = os.path.basename(f)[:-5]
    if base.endswith('.v1'):
        name, ver = base[:-3], 'V1'
    elif base.endswith('.v2'):
        name, ver = base[:-3], 'V2'
    else:
        continue
    d = json.load(open(f))
    pr.setdefault(name, {})[ver] = (d.get('cost'), d.get('depth'), d.get('status'))

rows = []
for f in sorted(glob.glob(fin_dir + '/*.json')):
    name = os.path.basename(f)[:-5]
    d = json.load(open(f))
    xd = (d.get('cost'), d.get('depth'), d.get('status'))
    p = pr.get(name)
    if not p:
        continue
    v1 = p.get('V1')
    v2 = p.get('V2')
    best = min([x for x in (v1, v2) if x], key=lambda t: (t[0], t[1]))
    rows.append((name, best, xd, v1, v2))

better_count = 0
worse_count = 0
tie_count = 0
total = 0
print(f"{'matrix':34s} {'PRbest':>8s} {'final':>8s} {'dc/dd':>7s}  {'V1':>9s} {'V2':>9s}")
for name, best, xd, v1, v2 in rows:
    if xd[2] != 'verified':
        continue
    total += 1
    dc = xd[0] - best[0]
    dd = xd[1] - best[1]
    tag = ''
    if dc < 0 or (dc == 0 and dd < 0):
        better_count += 1
        tag = 'BETTER'
    elif dc > 0 or (dc == 0 and dd > 0):
        worse_count += 1
        tag = 'WORSE'
    else:
        tie_count += 1
    print(f"{name:34s} {best[0]:>5d}/{best[1]:<2d} {xd[0]:>5d}/{xd[1]:<2d} {dc:+3d}/{dd:+d}  {v1[0]:>5d}/{v1[1]:<2d} {v2[0]:>5d}/{v2[1]:<2d} {tag}")

print(f"\nverified {total}/57: BETTER {better_count}, TIE {tie_count}, WORSE {worse_count}")
print("total gate delta vs PR best:",
      sum(x[1][0] - x[2][0] for x in rows if x[2][2] == 'verified'))
b1 = sum(1 for n, b, x, v1, v2 in rows if x[2] == 'verified' and (x[0], x[1]) < (v1[0], v1[1]))
b2 = sum(1 for n, b, x, v1, v2 in rows if x[2] == 'verified' and (x[0], x[1]) < (v2[0], v2[1]))
w1 = sum(1 for n, b, x, v1, v2 in rows if x[2] == 'verified' and (x[0], x[1]) > (v1[0], v1[1]))
w2 = sum(1 for n, b, x, v1, v2 in rows if x[2] == 'verified' and (x[0], x[1]) > (v2[0], v2[1]))
print(f"final vs V1: better {b1}, worse {w1} | final vs V2: better {b2}, worse {w2}")
