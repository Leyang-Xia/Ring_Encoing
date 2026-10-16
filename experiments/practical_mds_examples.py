"""Reconstruct and verify systematic RS and Cauchy encoding experiments.

Only the Python standard library is needed for verification. --profile additionally
compiles the repository's unmodified X-Sets.cpp and preserves its result or timeout.
"""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random
import subprocess
import tempfile

from ring_lift_example import GreedyCircuit, remainder, verify_circuit
from structured_ring_example import evaluate, expand, full_rank, multiply, power
from structured_ring_family import check_record, run

CERTIFICATE = Path(__file__).with_suffix('.json')


def inverse_matrix(matrix, w, P):
    n = len(matrix)
    aug = [list(row) + [int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if aug[i][j])
        aug[j], aug[pivot] = aug[pivot], aug[j]
        inverse = power(aug[j][j], (1 << w) - 2, w, P)
        aug[j] = [multiply(a, inverse, w, P) for a in aug[j]]
        for i in range(n):
            if i != j:
                coefficient = aug[i][j]
                aug[i] = [a ^ multiply(coefficient, b, w, P) for a, b in zip(aug[i], aug[j])]
    assert [row[:n] for row in aug] == [[int(i == j) for j in range(n)] for i in range(n)]
    return [row[n:] for row in aug]


def matrix_product(A, B, w, P):
    result = []
    for row in A:
        output = []
        for column in zip(*B):
            value = 0
            for a, b in zip(row, column):
                value ^= multiply(a, b, w, P)
            output.append(value)
        result.append(output)
    return result


def cauchy(x, y, w, P):
    assert len(set(x + y)) == len(x + y)
    return [[power(a ^ b, (1 << w) - 2, w, P) for b in y] for a in x]


def source_matrix(record):
    w, P, k, m = (record[key] for key in ('w', 'P', 'k', 'm'))
    if record['kind'] == 'cauchy':
        definition = record['definition']
        A = cauchy(definition['x'], definition['y'], w, P)
        assert len(A) == m and len(A[0]) == k
        return A
    assert record['kind'] in ['systematic_rs', 'systematic_rs_moments']
    points = record['definition']['points']
    assert len(points) == len(set(points)) == k + m
    V = [[power(a, j, w, P) for j in range(k)] for a in points]
    G = matrix_product(V, inverse_matrix(V[:k], w, P), w, P)
    assert G[:k] == [[int(i == j) for j in range(k)] for i in range(k)]
    A = G[k:]
    # Independently derive the same coefficients from Lagrange interpolation.
    row_weights, column_weights = [], []
    for b in points[k:]:
        value = 1
        for a in points[:k]:
            value = multiply(value, a ^ b, w, P)
        row_weights.append(value)
    for j, a in enumerate(points[:k]):
        value = 1
        for other, b in enumerate(points[:k]):
            if other != j:
                value = multiply(value, a ^ b, w, P)
        column_weights.append(power(value, (1 << w) - 2, w, P))
    C = cauchy(points[k:], points[:k], w, P)
    assert A == [[multiply(row_weights[i], multiply(C[i][j], column_weights[j], w, P), w, P)
                  for j in range(k)] for i in range(m)]
    return A


def matrices(record):
    w, P, k, m = (record[key] for key in ('w', 'P', 'k', 'm'))
    A = source_matrix(record)
    assert A == record['A']
    M, lifted = record['M'], record['lifted']
    n = M.bit_length() - 1
    assert n == record['n'] and remainder(M, P) == 0
    field_main = record.get('main_field_matrix', A)
    post = record.get('post_field_matrix', [[int(i == j) for j in range(m)] for i in range(m)])
    assert matrix_product(post, field_main, w, P) == A
    assert [[remainder(a, P) for a in row] for row in lifted] == field_main
    if record['kind'] == 'systematic_rs_moments':
        points = record['definition']['points']
        H = [[power(a, i + 1, w, P) for a in points] for i in range(m)]
        assert field_main == [row[:k] for row in H]
        assert post == inverse_matrix([row[k:] for row in H], w, P)
    T = expand(A, w, w, lambda a, b: multiply(a, b, w, P))
    B = expand(lifted, w, n, lambda a, b: multiply(a, b, n, M))
    virtual = list(dict.fromkeys(mask for mask in B if mask))
    index = {mask: i for i, mask in enumerate(virtual)}
    D = []
    for output in range(m):
        for bit in range(w):
            mask = 0
            for i in range(m):
                for j in range(n):
                    coefficient = multiply(post[output][i], remainder(1 << j, P), w, P)
                    if B[i * n + j] and (coefficient >> bit) & 1:
                        mask ^= 1 << index[B[i * n + j]]
            D.append(mask)
    for target, row in zip(T, D):
        value = 0
        for j, mask in enumerate(virtual):
            if (row >> j) & 1:
                value ^= mask
        assert value == target
    return A, T, B, virtual, D


def compose_counted(main, projection, virtual, middle, targets):
    circuit = GreedyCircuit(main['inputs'])
    mapping = dict(enumerate(range(main['inputs'])))
    if main['zero_node'] is not None:
        mapping[main['zero_node']] = circuit.input_count
    for node, a, b in main['ops']:
        mapping[node] = circuit.xor(mapping[a], mapping[b])
    z = [mapping[node] for node in main['outputs']]
    main_nodes = set(circuit.ops)
    by_mask = {mask: z[i] for i, mask in enumerate(middle) if mask}
    mapping = {i: by_mask[mask] for i, mask in enumerate(virtual)}
    if projection['zero_node'] is not None:
        mapping[projection['zero_node']] = circuit.input_count
    for node, a, b in projection['ops']:
        mapping[node] = circuit.xor(mapping[a], mapping[b])
    result = circuit.export([mapping[node] for node in projection['outputs']])
    verify_circuit(result, targets)
    live = {op[0] for op in result['ops']}
    result['live_main_xors'] = len(live & main_nodes)
    result['live_projection_xors'] = len(live - main_nodes)
    assert result['cost'] == result['live_main_xors'] + result['live_projection_xors']
    return result


def verify(record, samples):
    w, P, k, m = (record[key] for key in ('w', 'P', 'k', 'm'))
    assert all(remainder(P, q) for degree in range(1, w // 2 + 1)
               for q in range(1 << degree, 1 << (degree + 1)))
    A, T, B, virtual, D = matrices(record)
    assert virtual == record['projection_input_masks']
    verify_circuit(record['main'], B)
    verify_circuit(record['projection'], D)
    verify_circuit(record['ring'], T)
    rebuilt = compose_counted(record['main'], record['projection'], virtual, B, T)
    assert rebuilt['ops'] == record['ring']['ops']
    assert rebuilt['outputs'] == record['ring']['outputs']
    assert rebuilt['cost'] == record['ring']['cost']
    verify_circuit(record['direct_cse'], T)
    for measurement in record['measurements']:
        check_record(measurement, {'direct': T, 'main': B, 'projection': D}[measurement['path']])
    baseline = min([record['direct_cse']['cost']] + [r['cost'] for r in record['measurements']
                   if r['path'] == 'direct' and r['status'] == 'verified'])
    assert baseline == record['best_completed_direct']
    minors = 0
    for size in range(1, min(k, m) + 1):
        for rows in combinations(range(m), size):
            for columns in combinations(range(k), size):
                assert full_rank([[A[i][j] for j in columns] for i in rows], w, P)
                minors += 1
    assert minors == record['checked_nonzero_minors']
    rng = random.Random(20261017)
    for _ in range(samples):
        data = [rng.randrange(1 << w) for _ in range(k)]
        expected = []
        for row in A:
            value = 0
            for a, d in zip(row, data):
                value ^= multiply(a, d, w, P)
            expected.extend((value >> bit) & 1 for bit in range(w))
        assert evaluate(record['ring'], sum(d << (w * j) for j, d in enumerate(data))) == expected
    status = record.get('comparison_status', 'completed_heuristic_comparison')
    return dict(kind=record['kind'], comparison_status=status, systematic_code=[k + m, k, m + 1],
                main_xors=rebuilt['live_main_xors'], projection_xors=rebuilt['live_projection_xors'],
                ring=rebuilt['cost'], best_completed_direct=baseline,
                net_saved=baseline-rebuilt['cost'] if status == 'completed_heuristic_comparison' else None,
                nonzero_minors=minors, random_vectors=samples)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', help='example name from the saved certificate')
    parser.add_argument('--path', choices=['direct', 'main', 'projection'], default='direct')
    parser.add_argument('--method', choices=['MW_SQ', 'UBER_XSET'], default='MW_SQ')
    parser.add_argument('--params', nargs=3, type=int, default=[0, 3, 0])
    parser.add_argument('--timeout', type=float, default=180)
    parser.add_argument('--samples', type=int, default=1000)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.timeout <= 0 or args.samples < 0:
        parser.error('timeout must be positive and samples nonnegative')
    data = json.loads(CERTIFICATE.read_text())
    if not args.profile:
        print(json.dumps({name: verify(r, args.samples) for name, r in data['examples'].items()}, indent=2))
    else:
        if args.output is None:
            parser.error('--profile requires --output')
        record = data['examples'][args.profile]
        _, T, B, V, D = matrices(record)
        rows, inputs = {'direct': (T, record['w'] * record['k']),
                        'main': (B, record['w'] * record['k']), 'projection': (D, len(V))}[args.path]
        source = Path(__file__).resolve().parents[1] / 'X-Sets.cpp'
        with tempfile.TemporaryDirectory(prefix='practical-mds-') as directory:
            binary = Path(directory) / 'X-Sets'
            subprocess.run(['g++', '-std=c++17', '-O3', str(source), '-o', str(binary)], check=True)
            result = run(binary, rows, inputs, args.method, args.params, args.timeout, args.path)
        result['xsets_source_sha256'] = sha256(source.read_bytes()).hexdigest()
        args.output.write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps({key: value for key, value in result.items() if key != 'circuit'}, indent=2))
