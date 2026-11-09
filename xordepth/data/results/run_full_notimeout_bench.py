import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
XORDEPTH = os.path.join(ROOT, 'build', 'xordepth')

CONFIG = [
    '--beam', '4', '--branch', '4', '--representations', '4',
    '--restarts', '2', '--greedy-restarts', '2', '--greedy-mwsq',
    '--rewrite-trials', '400', '--max-states', '256', '--seed', '1',
]


def main(out_dir, extra_args=()):
    os.makedirs(out_dir, exist_ok=True)
    ref_dir = os.path.join(ROOT, 'data', 'results', 'xordepth_full_mwsq_20261106')
    manifest = os.path.join(ROOT, 'data', 'matrices', 'full.manifest')
    base = os.path.dirname(manifest)
    with open(manifest) as f:
        entries = [ln.strip() for ln in f
                   if ln.strip() and not ln.strip().startswith('#')]
    rows = []
    total = 0.0
    for rel in entries:
        name = os.path.splitext(os.path.basename(rel))[0]
        mat = os.path.join(base, rel)
        ref_file = os.path.join(ref_dir, name + '.json')
        depth = json.load(open(ref_file)).get('depth_limit')
        out = os.path.join(out_dir, name + '.json')
        cmd = [XORDEPTH, 'bounded', mat, '--depth', str(depth)] + CONFIG \
            + list(extra_args) + ['--no-timeout', '--output', out]
        t0 = time.monotonic()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        wall = time.monotonic() - t0
        total += wall
        ref = json.load(open(ref_file))
        if proc.returncode != 0:
            row = (name, None, None, 'ERROR', wall, ref.get('cost'), ref.get('depth'))
        else:
            j = json.load(open(out))
            row = (name, j.get('cost'), j.get('depth'), j.get('status'), wall,
                   ref.get('cost'), ref.get('depth'))
        rows.append(row)
        print(f"{row[0]:34s} cost={row[1]} depth={row[2]} {row[3]:10s} "
              f"wall={row[4]:6.2f}s ref={row[5]}/{row[6]}", flush=True)
    print(f"TOTAL wall = {total:.2f}s over {len(rows)} matrices")
    same = sum(1 for r in rows if r[1] == r[5] and r[2] == r[6])
    diff = [(r[0], r[5], r[6], r[1], r[2]) for r in rows if (r[1], r[2]) != (r[5], r[6])]
    print(f"same (cost,depth) vs 15s reference: {same}/{len(rows)}")
    for d in diff:
        print("DIFF", d)
    with open(os.path.join(out_dir, 'bench.tsv'), 'w') as f:
        f.write('name\tcost\tdepth\tstatus\twall_seconds\tref_cost\tref_depth\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
