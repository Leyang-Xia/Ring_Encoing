import json
import glob
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = os.path.join(ROOT, 'data', 'results')

# Compare full result JSONs across runs: cost/depth/status plus the complete
# circuit ops list, so "bit-exact" here means the exported circuits agree
# gate for gate, not just on the (cost, depth) summary.
IGNORED = {'wall_seconds', 'elapsed', 'timeout_seconds'}


def load(run, name):
    path = os.path.join(RES, run, name + '.json')
    with open(path) as f:
        return json.load(f)


def fingerprint(j):
    ops = j['circuit']['ops'] if 'circuit' in j else None
    return (j.get('cost'), j.get('depth'), j.get('status'), ops)


def compare(run_a, run_b, label):
    mismatch = 0
    checked = 0
    for path in sorted(glob.glob(os.path.join(RES, run_a, '*.json'))):
        name = os.path.splitext(os.path.basename(path))[0]
        if not os.path.exists(os.path.join(RES, run_b, name + '.json')):
            continue
        checked += 1
        if fingerprint(load(run_a, name)) != fingerprint(load(run_b, name)):
            a, b = load(run_a, name), load(run_b, name)
            mismatch += 1
            print(f"  MISMATCH {name}: {a.get('cost')}/{a.get('depth')} vs "
                  f"{b.get('cost')}/{b.get('depth')}")
    print(f"{label}: {checked} compared, {mismatch} mismatched")
    return mismatch


if __name__ == '__main__':
    bad = 0
    bad += compare('notimeout_opt7_20261109', 'notimeout_threads1_20261109',
                   'pre-change opt7 vs new threads=1')
    bad += compare('notimeout_threads1_20261109', 'notimeout_threads4_20261109',
                   'new threads=1 vs new threads=4')
    sys.exit(1 if bad else 0)
