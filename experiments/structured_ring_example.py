"""Verify the F256 -> 17-dimensional cyclic-ring example (standard library only).

Optional --rerun-projection/--rerun-baselines compile unmodified X-Sets.cpp using g++.
Default verification reproduces CSE, proves matrix equivalence and checks random data.
"""

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import random
import re
import subprocess
import tempfile

from ring_lift_example import GreedyCircuit, polynomial_product, remainder, verify_circuit


def multiply(a, b, width, modulus):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & (1 << width):
            a ^= modulus
    return result


def cyclic_multiply(a, b, width):
    result = 0
    for shift in range(width):
        if (a >> shift) & 1:
            result ^= ((b << shift) | (b >> (width - shift))) & ((1 << width) - 1)
    return result


def expand(matrix, input_width, output_width, product):
    result = []
    for row in matrix:
        columns = [product(a, 1 << bit) for a in row for bit in range(input_width)]
        result.extend(sum(((column >> bit) & 1) << j for j, column in enumerate(columns))
                      for bit in range(output_width))
    return result


def power(a, exponent, width, modulus):
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, a, width, modulus)
        a = multiply(a, a, width, modulus)
        exponent >>= 1
    return result


def full_rank(matrix, width, modulus):
    matrix = [list(row) for row in matrix]
    for column in range(len(matrix)):
        pivot = next((r for r in range(column, len(matrix)) if matrix[r][column]), None)
        if pivot is None:
            return False
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        inverse = power(matrix[column][column], (1 << width) - 2, width, modulus)
        matrix[column] = [multiply(v, inverse, width, modulus) for v in matrix[column]]
        for r in range(column + 1, len(matrix)):
            factor = matrix[r][column]
            matrix[r] = [a ^ multiply(factor, b, width, modulus)
                         for a, b in zip(matrix[r], matrix[column])]
    return True


def evaluate(circuit, word):
    values = {i: (word >> i) & 1 for i in range(circuit['inputs'])}
    if circuit['zero_node'] is not None:
        values[circuit['zero_node']] = 0
    for node, a, b in circuit['ops']:
        values[node] = values[a] ^ values[b]
    return [values[node] for node in circuit['outputs']]


def xsets(binary, rows, inputs, method, params, timeout):
    """Remove free outputs, restore original output order, and verify every gate."""
    tool_rows = sorted(set(row for row in rows if row.bit_count() > 1))
    text = '\n'.join(''.join(str((row >> j) & 1) for j in range(inputs))
                     for row in tool_rows) + '\n'
    run = subprocess.run([str(binary), *map(str, params), method], input=text,
                         text=True, capture_output=True, check=True, timeout=timeout)
    values = {i: 1 << i for i in range(inputs)}
    values[-1] = 0
    ops, printed_outputs = [], {}
    for line in run.stdout.splitlines():
        match = re.fullmatch(r's(\d+) = s(\d+) \+ s(\d+)', line)
        if match:
            node, a, b = map(int, match.groups())
            assert node not in values
            values[node] = values[a] ^ values[b]
            ops.append([node, a, b])
        match = re.fullmatch(r'b(\d+) = s(\d+)', line)
        if match:
            output, node = map(int, match.groups())
            printed_outputs[output] = node
    assert [values[printed_outputs[i]] for i in range(len(tool_rows))] == tool_rows
    assert len(ops) == int(re.search(r'Total XOR operations: (\d+)', run.stdout).group(1))
    by_mask = {mask: node for node, mask in values.items()}
    outputs = [by_mask[row] for row in rows]
    parents = {node: (a, b) for node, a, b in ops}
    live = set()

    def visit(node):
        if node in live or node not in parents:
            return
        live.add(node)
        for parent in parents[node]:
            visit(parent)

    for node in outputs:
        visit(node)
    circuit = dict(inputs=inputs, zero_node=-1, outputs=outputs,
                   ops=[op for op in ops if op[0] in live], cost=len(live))
    verify_circuit(circuit, rows)
    return circuit


def verify(args):
    certificate = json.loads(Path(__file__).with_suffix('.json').read_text())
    w, n, P, Q = (certificate[key] for key in ('w', 'n', 'P', 'Q'))
    matrix, lifted = certificate['A'], certificate['lifted']
    k, m = certificate['k'], certificate['m']
    assert (w, n, k, m) == (8, 17, 8, 3)
    assert polynomial_product(P, Q) == (1 << n) | 1
    assert all(remainder(P, candidate)
               for degree in range(1, w // 2 + 1)
               for candidate in range(1 << degree, 1 << (degree + 1)))
    assert all(lifted[i][j] == 1 << (i * j) for i in range(m) for j in range(k))
    assert [[remainder(a, P) for a in row] for row in lifted] == matrix
    assert all(matrix[i][j] ^ polynomial_product(P, certificate['G'][i][j]) == lifted[i][j]
               for i in range(m) for j in range(k))
    field_product = lambda a, b: multiply(a, b, w, P)
    ring_product = lambda a, b: cyclic_multiply(a, b, n)
    targets = expand(matrix, w, w, field_product)
    middle = expand(lifted, w, n, ring_product)
    D = [power(2, bit, w, P) for bit in range(n)]
    assert D == [remainder(1 << bit, P) for bit in range(n)]
    projected = []
    for symbol in range(m):
        for bit in range(w):
            value = 0
            for j in range(n):
                if (D[j] >> bit) & 1:
                    value ^= middle[n * symbol + j]
            projected.append(value)
    assert projected == targets

    main = GreedyCircuit(w * k)
    z = main.synthesize_masks(middle)
    assert [sum(bool(row) for row in middle[n * i:n * (i + 1)]) for i in range(m)] == [8, 15, 17]
    assert [sum(row.bit_count() for row in middle[n * i:n * (i + 1)]) for i in range(m)] == [64, 64, 64]
    assert len(main.ops) == 152
    assert sum(max(row.bit_count() - 1, 0) for row in middle) == 152
    main_values = verify_circuit(certificate['ring'], targets)
    assert [main_values[node] for node in certificate['ring']['middle_outputs']] == middle
    assert certificate['ring']['main_generated'] == 152
    assert certificate['ring']['projection_generated'] == 51
    assert certificate['ring']['cost'] == 203

    virtual_inputs = certificate['projection_input_masks']
    assert len(virtual_inputs) == len(set(virtual_inputs)) == 40
    virtual_index = {mask: i for i, mask in enumerate(virtual_inputs)}
    projection_targets = []
    for symbol in range(m):
        for bit in range(w):
            value = 0
            for j in range(n):
                mask = middle[n * symbol + j]
                if mask and (D[j] >> bit) & 1:
                    value ^= 1 << virtual_index[mask]
            projection_targets.append(value)
    verify_circuit(certificate['projection'], projection_targets)
    assert certificate['projection']['cost'] == 51
    mapping = {i: main.by_mask[mask] for i, mask in enumerate(virtual_inputs)}
    mapping[-1] = main.input_count
    for node, a, b in certificate['projection']['ops']:
        mapping[node] = main.xor(mapping[a], mapping[b])
    outputs = [mapping[node] for node in certificate['projection']['outputs']]
    rebuilt = main.export(outputs)
    assert rebuilt['ops'] == certificate['ring']['ops']
    assert rebuilt['outputs'] == certificate['ring']['outputs']

    histogram = Counter()
    for seed in range(100):
        c = GreedyCircuit(w * k)
        outputs = c.synthesize_masks(targets, seed)
        circuit = c.export(outputs)
        verify_circuit(circuit, targets)
        histogram[circuit['cost']] += 1
        if seed == certificate['direct_cse_seed']:
            assert circuit['ops'] == certificate['direct_cse']['ops']
    assert {str(cost): count for cost, count in histogram.items()} == certificate['direct_cse_histogram']
    assert min(histogram) == 217
    verify_circuit(certificate['direct_cse'], targets)
    verify_circuit(certificate['ring_cse'], targets)
    for record in certificate['direct_xsets'].values():
        verify_circuit(record['circuit'], targets)
    for key, expected_cost in [('ring_mw_pipeline', 203), ('ring_uber_pipeline', 204)]:
        pipeline = certificate[key]
        verify_circuit(pipeline['main'], middle)
        assert pipeline['projection_input_masks'] == virtual_inputs
        verify_circuit(pipeline['projection'], projection_targets)
        values = verify_circuit(pipeline['ring'], targets)
        assert [values[node] for node in pipeline['ring']['middle_outputs']] == middle
        assert pipeline['main']['cost'] == 152
        assert pipeline['ring']['cost'] == expected_cost
    for record in certificate['direct_uber'].values():
        if record['status'] == 'verified':
            verify_circuit(record['circuit'], targets)
    checked_minors = 0
    for size in range(1, m + 1):
        for row_indices in combinations(range(m), size):
            for columns in combinations(range(k), size):
                assert full_rank([[matrix[i][j] for j in columns] for i in row_indices], w, P)
                checked_minors += 1
    assert checked_minors == 164

    rng = random.Random(20261016)
    for _ in range(args.samples):
        data = [rng.randrange(1 << w) for _ in range(k)]
        expected = []
        ring_result = []
        for i in range(m):
            field_value = ring_value = 0
            for j in range(k):
                field_value ^= field_product(matrix[i][j], data[j])
                ring_value ^= ring_product(lifted[i][j], data[j])
            expected.extend((field_value >> bit) & 1 for bit in range(w))
            ring_result.append(remainder(ring_value, P))
            assert ring_result[-1] == field_value
        word = sum(value << (w * j) for j, value in enumerate(data))
        assert evaluate(certificate['ring'], word) == expected
        for record in certificate['direct_xsets'].values():
            assert evaluate(record['circuit'], word) == expected

    report = dict(matrix_identity=True, all_saved_circuits_verified=True,
                  main_xors=152, projection_xors=51, ring_total=203,
                  direct_cse_100_seeds=min(histogram),
                  direct_xsets={name: record['circuit']['cost']
                                for name, record in certificate['direct_xsets'].items()},
                  complete_ring_pipelines={'MW_SQ': 203, 'UBER_XSET': 204},
                  checked_nonzero_minors=checked_minors, systematic_code=[11, 8, 4],
                  random_data_vectors=args.samples)
    if args.rerun_projection or args.rerun_baselines:
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory(prefix='structured-ring-') as directory:
            binary = Path(directory) / 'X-Sets'
            subprocess.run(['g++', '-std=c++17', '-O3', str(root / 'X-Sets.cpp'),
                            '-o', str(binary)], check=True)
            if args.rerun_projection:
                c = xsets(binary, projection_targets, 40, 'MW_SQ', [0, 3, 0], args.timeout)
                assert c['cost'] == 51
                report['projection_rerun'] = c['cost']
            if args.rerun_baselines:
                report['baseline_reruns'] = {}
                for method, record in certificate['direct_xsets'].items():
                    try:
                        c = xsets(binary, targets, w * k, method, record['params'], args.timeout)
                        report['baseline_reruns'][method] = c['cost']
                    except subprocess.TimeoutExpired:
                        report['baseline_reruns'][method] = 'timeout / inconclusive'
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, default=1000)
    parser.add_argument('--rerun-projection', action='store_true')
    parser.add_argument('--rerun-baselines', action='store_true')
    parser.add_argument('--timeout', type=float, default=180)
    verify(parser.parse_args())
