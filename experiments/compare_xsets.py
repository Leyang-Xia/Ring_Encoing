"""Compile unmodified X-Sets and compare two methods on the fixed F16 example.

Example: python3 experiments/compare_xsets.py --output /tmp/xsets-results.json
Requires g++; every emitted circuit is independently checked symbolically.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile
import time

from ring_lift_example import cyclic_product, field_product, matrix_masks, verify_circuit


PROFILES = [('default', [10, 10, 100]), ('enum2', [10, 2, 0]),
            ('enum3', [10, 3, 0]), ('enum3-thresh0', [0, 3, 0])]


def parse_schedule(output, targets):
    signals = {i: 1 << i for i in range(16)}
    ops = []
    outputs = {}
    for line in output.splitlines():
        match = re.fullmatch(r's(\d+) = s(\d+) \+ s(\d+)', line)
        if match:
            node, a, b = map(int, match.groups())
            assert node not in signals
            signals[node] = signals[a] ^ signals[b]
            ops.append([node, a, b])
        match = re.fullmatch(r'b(\d+) = s(\d+)', line)
        if match:
            target, node = map(int, match.groups())
            outputs[target] = node
    reported = int(re.search(r'Total XOR operations: (\d+)', output).group(1))
    assert len(ops) == reported
    assert [signals[outputs[i]] for i in range(len(targets))] == targets
    parents = {node: (a, b) for node, a, b in ops}
    live = set()

    def visit(node):
        if node in live or node not in parents:
            return
        live.add(node)
        for parent in parents[node]:
            visit(parent)

    for node in outputs.values():
        visit(node)
    circuit = dict(inputs=16, zero_node=None,
                   ops=[op for op in ops if op[0] in live],
                   outputs=[outputs[i] for i in range(len(targets))], cost=len(live))
    verify_circuit(circuit, targets)
    return reported, circuit


def append_projection(main):
    ops = [list(op) for op in main['ops']]
    outputs = []
    node = max([15] + [op[0] for op in ops]) + 1
    for symbol in range(2):
        for bit in range(4):
            ops.append([node, main['outputs'][5 * symbol + bit],
                        main['outputs'][5 * symbol + 4]])
            outputs.append(node)
            node += 1
    return dict(inputs=16, zero_node=None, ops=ops, outputs=outputs, cost=len(ops))


def compare(output_path, timeout):
    root = Path(__file__).resolve().parents[1]
    source = root / 'X-Sets.cpp'
    certificate = json.loads(Path(__file__).with_name('ring_lift_example.json').read_text())
    targets = matrix_masks(certificate['A'], 4, field_product)
    middle = matrix_masks(certificate['ring']['lifted'], 5, cyclic_product)
    result = dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  source='X-Sets.cpp', source_modified=False,
                  compiler=subprocess.check_output(['g++', '--version'], text=True).splitlines()[0],
                  compile_flags=['-std=c++17', '-O3'], timeout_seconds=timeout,
                  projection_cost=8, records=[])
    with tempfile.TemporaryDirectory(prefix='ring-xsets-') as directory:
        binary = Path(directory) / 'X-Sets'
        subprocess.run(['g++', '-std=c++17', '-O3', str(source), '-o', str(binary)], check=True)
        for profile, params in PROFILES:
            for technique in ('MW_SQ', 'UBER_XSET'):
                for path, wanted in [('direct', targets), ('ring', middle)]:
                    text = '\n'.join(''.join(str((mask >> j) & 1) for j in range(16))
                                     for mask in wanted) + '\n'
                    record = dict(profile=profile, params=params, technique=technique, path=path)
                    started = time.perf_counter()
                    try:
                        run = subprocess.run([str(binary), *map(str, params), technique],
                                             input=text, capture_output=True, text=True,
                                             check=True, timeout=timeout)
                        reported, circuit = parse_schedule(run.stdout, wanted)
                        full = append_projection(circuit) if path == 'ring' else circuit
                        verify_circuit(full, targets)
                        record.update(status='verified', reported=reported, circuit=circuit,
                                      full_circuit=full, total_cost=full['cost'])
                    except subprocess.TimeoutExpired:
                        record['status'] = 'timeout'
                    record['seconds'] = time.perf_counter() - started
                    result['records'].append(record)
    if output_path:
        output_path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps([dict(profile=r['profile'], method=r['technique'], path=r['path'],
                           status=r['status'], total=r.get('total_cost'))
                      for r in result['records']], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--timeout', type=float, default=20)
    args = parser.parse_args()
    compare(args.output, args.timeout)
