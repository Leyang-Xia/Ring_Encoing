import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
XORDEPTH = os.path.join(ROOT, 'build', 'xordepth')

# Sample protocol: explicit D per matrix (row-weight depth lower bound used in
# the PR sample comparison), seed=1, no wall-clock deadline.
D_TABLE = {
    'dense_10x10_seed123.txt': 3,
    'dense_40x64_seed123.txt': 6,
    'sparse_12x20_seed123.txt': 2,
    'sparse_100x100_seed123.txt': 2,
    'blocks_24x40_seed123.txt': 3,
    'blocks_100x100_seed456.txt': 5,
    'modular_10x10_seed456.txt': 3,
    'modular_64x96_seed456.txt': 6,
    'cauchy_7_4_direct.txt': 5,
    'systematic_rs_7_4_direct.txt': 5,
    'vandermonde_8to17_direct.txt': 5,
    'gf8_10x10_expanded.txt': 6,
}

CONFIG = [
    '--beam', '4', '--branch', '4', '--representations', '4',
    '--restarts', '2', '--greedy-restarts', '2', '--greedy-mwsq',
    '--rewrite-trials', '400', '--max-states', '256', '--seed', '1',
]


def run_all(out_dir, extra_args=()):
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    total = 0.0
    manifest = os.path.join(ROOT, 'data', 'matrices', 'sample.manifest')
    base = os.path.dirname(manifest)
    with open(manifest) as f:
        entries = [ln.strip() for ln in f
                   if ln.strip() and not ln.strip().startswith('#')]
    for rel in entries:
        name = os.path.splitext(os.path.basename(rel))[0]
        mat = os.path.join(base, rel)
        depth = D_TABLE[os.path.basename(rel)]
        out = os.path.join(out_dir, name + '.json')
        cmd = [XORDEPTH, 'bounded', mat, '--depth', str(depth)] + CONFIG \
            + list(extra_args) + ['--no-timeout', '--output', out]
        t0 = time.monotonic()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        wall = time.monotonic() - t0
        total += wall
        if proc.returncode != 0:
            row = (name, None, None, 'ERROR', wall, proc.stderr.strip()[:200])
        else:
            j = json.load(open(out))
            row = (name, j.get('cost'), j.get('depth'), j.get('status'), wall, '')
        rows.append(row)
        print(f"{row[0]:34s} cost={row[1]} depth={row[2]} {row[3]:12s} "
              f"wall={row[4]:7.2f}s {row[5]}", flush=True)
    print(f"TOTAL wall = {total:.2f}s over {len(rows)} matrices")
    with open(os.path.join(out_dir, 'bench.tsv'), 'w') as f:
        f.write('name\tcost\tdepth\tstatus\twall_seconds\n')
        for r in rows:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.3f}\n")
    return total, rows


if __name__ == '__main__':
    out_dir = sys.argv[1]
    extra = sys.argv[2:]
    run_all(out_dir, extra)
