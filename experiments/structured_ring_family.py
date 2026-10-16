"""Verify and extend the structured ring-lift experiments (Python 3.10+).

Default: verify saved circuits, matrix identities, MDS minors, and random inputs.
--profile NAME: compile unmodified X-Sets.cpp and rerun a selected optimization.
--search: refine additional recorded matrix configurations, saving every outcome.
"""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random
import subprocess
import tempfile
import time

from ring_lift_example import GreedyCircuit, remainder, verify_circuit
from structured_ring_example import evaluate, expand, full_rank, multiply, power, xsets

CERTIFICATE = Path(__file__).with_suffix('.json')


def synthesize(targets, inputs, seed=0):
    circuit = GreedyCircuit(inputs)
    outputs = circuit.synthesize_masks(targets, seed)
    result = circuit.export(outputs)
    verify_circuit(result, targets)
    return result


def coefficients(case):
    w, n, P, k, m = (case[key] for key in ('w', 'n', 'P', 'k', 'm'))
    slopes = case.get('slopes', list(range(m)))
    shifts = case.get('row_shifts', [0] * m)
    column_shifts = case.get('column_shifts', [0] * k)
    exponents = [[(slopes[i] * e + shifts[i] + column_shifts[j]) % n
                  for j, e in enumerate(case['exponents'])] for i in range(m)]
    lifted = [[1 << exponent for exponent in row] for row in exponents]
    A = [[remainder(a, P) for a in row] for row in lifted]
    return A, lifted


def is_vandermonde(A, width, modulus):
    """Recognize actual unscaled columns (1, x, x^2, ...), not construction labels."""
    if len(A) < 2 or not A[0] or any(a != 1 for a in A[0]):
        return False
    points = A[1]
    return (all(points) and len(set(points)) == len(points) and
            all(row == [power(x, i, width, modulus) for x in points]
                for i, row in enumerate(A)))


def matrices(case):
    """Build both paths independently; verify the composed binary matrices."""
    w, n, P, k, m = (case[key] for key in ('w', 'n', 'P', 'k', 'm'))
    modulus = case.get('ring_modulus', (1 << n) | 1)
    assert modulus.bit_length() == n + 1 and remainder(modulus, P) == 0
    A, lifted = coefficients(case)
    targets = expand(A, w, w, lambda a, b: multiply(a, b, w, P))
    middle = expand(lifted, w, n, lambda a, b: multiply(a, b, n, modulus))
    virtual = list(dict.fromkeys(mask for mask in middle if mask))
    index = {mask: i for i, mask in enumerate(virtual)}
    projection = []
    for symbol in range(m):
        for bit in range(w):
            mask = 0
            for j in range(n):
                if middle[n * symbol + j] and (remainder(1 << j, P) >> bit) & 1:
                    mask ^= 1 << index[middle[n * symbol + j]]
            projection.append(mask)
    composed = []
    for row in projection:
        value = 0
        for j, mask in enumerate(virtual):
            if (row >> j) & 1:
                value ^= mask
        composed.append(value)
    assert composed == targets
    return A, lifted, targets, middle, virtual, projection


def compose(main, projection, virtual, middle, targets):
    circuit = GreedyCircuit(main['inputs'])
    mapping = {i: i for i in range(main['inputs'])}
    if main['zero_node'] is not None:
        mapping[main['zero_node']] = circuit.input_count
    for node, a, b in main['ops']:
        mapping[node] = circuit.xor(mapping[a], mapping[b])
    z = [mapping[node] for node in main['outputs']]
    by_mask = {mask: z[i] for i, mask in enumerate(middle) if mask}
    mapping = {i: by_mask[mask] for i, mask in enumerate(virtual)}
    if projection['zero_node'] is not None:
        mapping[projection['zero_node']] = circuit.input_count
    for node, a, b in projection['ops']:
        mapping[node] = circuit.xor(mapping[a], mapping[b])
    result = circuit.export([mapping[node] for node in projection['outputs']])
    verify_circuit(result, targets)
    result.update(main_cost=main['cost'], projection_cost=projection['cost'])
    return result


def check_record(record, targets):
    if record['status'] == 'verified':
        assert record['cost'] == record['circuit']['cost']
        verify_circuit(record['circuit'], targets)
    else:
        assert record['status'] == 'timeout' and record['cost'] is None
        assert 'circuit' not in record


def compose_layers(stages, targets):
    circuit = GreedyCircuit(stages[0]['inputs'])
    values = list(range(circuit.input_count))
    for stage in stages:
        assert stage['inputs'] == len(values)
        mapping = dict(enumerate(values))
        if stage['zero_node'] is not None:
            mapping[stage['zero_node']] = circuit.input_count
        for node, a, b in stage['ops']:
            mapping[node] = circuit.xor(mapping[a], mapping[b])
        values = [mapping[node] for node in stage['outputs']]
    result = circuit.export(values)
    verify_circuit(result, targets)
    return result


def verify_example(record, samples):
    case = record['case']
    w, k, m, P = (case[key] for key in ('w', 'k', 'm', 'P'))
    assert all(remainder(P, q) for degree in range(1, w // 2 + 1)
               for q in range(1 << degree, 1 << (degree + 1)))
    A, lifted, T, B, V, D = matrices(case)
    assert A == record['A'] and lifted == record['lifted']
    vandermonde = is_vandermonde(A, w, P)
    if 'matrix_family' in record:
        assert record['matrix_family'] == ('vandermonde' if vandermonde else 'scaled_extension')
    assert V == record['projection_input_masks']
    verify_circuit(record['main'], B)
    verify_circuit(record['projection'], D)
    verify_circuit(record['ring'], T)
    rebuilt = compose(record['main'], record['projection'], V, B, T)
    assert rebuilt['ops'] == record['ring']['ops']
    assert rebuilt['outputs'] == record['ring']['outputs']
    assert rebuilt['cost'] == record['ring']['cost']
    verify_circuit(record['direct_cse'], T)
    for r in record['measurements']:
        check_record(r, {'projection': D, 'direct': T, 'main': B}[r['path']])
    factored = record.get('factored_baselines', [])
    if factored:
        original = json.loads(Path(__file__).with_name('structured_ring_example.json').read_text())
        assert k == 8 and m == 3 and P == original['P']
        post_coefficients = [remainder(1 << shift, P) for shift in case['row_shifts']]
        assert A == [[multiply(post_coefficients[i], a, w, P)
                      for a in original['A'][2-i]] for i in range(m)]
        core = dict(original['direct_xsets']['UBER_XSET']['circuit'])
        core['outputs'] = core['outputs'][16:24] + core['outputs'][8:16] + core['outputs'][:8]
        core_targets = expand(list(reversed(original['A'])), w, w,
                              lambda a, b: multiply(a, b, w, P))
        verify_circuit(core, core_targets)
        post_matrix = [[post_coefficients[i] if i == j else 0 for j in range(m)]
                       for i in range(m)]
        post_targets = expand(post_matrix, w, w, lambda a, b: multiply(a, b, w, P))
    for r in factored:
        if 'circuit' in r:
            verify_circuit(r['post'], post_targets)
            verify_circuit(r['circuit'], T)
            rebuilt_factored = compose(core, r['post'], core_targets, core_targets, T)
            assert rebuilt_factored['ops'] == r['circuit']['ops']
            assert rebuilt_factored['outputs'] == r['circuit']['outputs']
    field = record.get('field_decomposition')
    if field:
        pre, core, post = field['matrices']
        assert all(pre[i][j] == 0 for i in range(k) for j in range(k) if i != j)
        assert all(post[i][j] == 0 for i in range(m) for j in range(m) if i != j)
        assert A == [[multiply(post[i][i], multiply(core[i][j], pre[j][j], w, P), w, P)
                      for j in range(k)] for i in range(m)]
        for stage, matrix in zip(field['stages'], field['matrices']):
            verify_circuit(stage, expand(matrix, w, w, lambda a, b: multiply(a, b, w, P)))
        verify_circuit(field['circuit'], T)
        rebuilt_field = compose_layers(field['stages'], T)
        assert rebuilt_field['ops'] == field['circuit']['ops']
        assert rebuilt_field['outputs'] == field['circuit']['outputs']
    baseline = min([record['direct_cse']['cost']] +
                   [r['cost'] for r in record['measurements']
                    if r['path'] == 'direct' and r['status'] == 'verified'] +
                   [r['circuit']['cost'] for r in factored if 'circuit' in r] +
                   ([field['circuit']['cost']] if field else []))
    assert baseline == record['best_completed_direct']
    assert abs(record['percent'] - 100 * (baseline - rebuilt['cost']) / baseline) < 1e-10
    minors = 0
    for size in range(1, m + 1):
        for rows in combinations(range(m), size):
            for columns in combinations(range(k), size):
                assert full_rank([[A[i][j] for j in columns] for i in rows], w, P)
                minors += 1
    assert minors == record['checked_nonzero_minors']
    rng = random.Random(20261016)
    for _ in range(samples):
        data = [rng.randrange(1 << w) for _ in range(k)]
        expected = []
        for row in A:
            value = 0
            for a, d in zip(row, data):
                value ^= multiply(a, d, w, P)
            expected.extend((value >> bit) & 1 for bit in range(w))
        word = sum(value << (w * j) for j, value in enumerate(data))
        assert evaluate(record['ring'], word) == expected
    return dict(matrix_family='vandermonde' if vandermonde else 'scaled_extension',
                ring=rebuilt['cost'], best_completed_direct=baseline,
                percent=record['percent'], nonzero_minors=minors, random_vectors=samples)


def verify(data, samples):
    primary = data['primary_case']
    _, _, T, B, _, D = matrices(primary)
    assert is_vandermonde(coefficients(primary)[0], primary['w'], primary['P'])
    for r in data['primary_additional_runs']:
        check_record(r, {'direct': T, 'main': B, 'projection': D}[r['path']])
        if 'ring' in r:
            verify_circuit(r['ring'], T)
    reports = {name: verify_example(record, samples)
               for name, record in data['examples'].items()}
    return dict(primary_additional_runs=len(data['primary_additional_runs']), examples=reports)


def run(binary, rows, inputs, method, params, timeout, path):
    record = dict(path=path, method=method, params=params, timeout=timeout)
    start = time.time()
    try:
        c = xsets(binary, rows, inputs, method, params, timeout)
        record.update(status='verified', cost=c['cost'], circuit=c)
    except subprocess.TimeoutExpired:
        record.update(status='timeout', cost=None)
    record['wall_seconds'] = time.time() - start
    return record


def profile(data, args, binary):
    if args.profile == 'original':
        case = data['primary_case']
    else:
        case = data['examples'][args.profile]['case']
    _, _, T, B, V, D = matrices(case)
    rows, inputs = {'direct': (T, case['w'] * case['k']),
                    'main': (B, case['w'] * case['k']),
                    'projection': (D, len(V))}[args.path]
    record = run(binary, rows, inputs, args.method, args.params, args.timeout, args.path)
    record.update(case=case)
    return record


def search(data, args, binary):
    """Refine a recorded, explicitly defined configuration; no hidden input costs."""
    cases = [case for phase in data['screening'].values() for case in phase]
    if args.family == 'vandermonde':
        cases = [case for case in cases
                 if is_vandermonde(coefficients(case)[0], case['w'], case['P'])]
    cases.sort(key=lambda case: case['percent'], reverse=True)
    results = []
    for case in cases[:args.limit]:
        A, lifted, T, B, V, D = matrices(case)
        minors = 0
        for size in range(1, case['m'] + 1):
            for rows in combinations(range(case['m']), size):
                for columns in combinations(range(case['k']), size):
                    assert full_rank([[A[i][j] for j in columns] for i in rows], case['w'], case['P'])
                    minors += 1
        main = synthesize(B, case['w'] * case['k'])
        projection = min((synthesize(D, len(V), seed) for seed in range(10)),
                         key=lambda c: c['cost'])
        direct = min((synthesize(T, case['w'] * case['k'], seed) for seed in range(20)),
                     key=lambda c: c['cost'])
        measurements = []
        for path, rows, inputs in [('projection', D, len(V)), ('direct', T, case['w'] * case['k'])]:
            for method in ['MW_SQ', 'UBER_XSET']:
                r = run(binary, rows, inputs, method, args.params, args.timeout, path)
                measurements.append(r)
                if path == 'projection' and r['status'] == 'verified' and r['cost'] < projection['cost']:
                    projection = r['circuit']
        ring = compose(main, projection, V, B, T)
        baseline = min([direct['cost']] + [r['cost'] for r in measurements
                       if r['path'] == 'direct' and r['status'] == 'verified'])
        result = dict(case=case, A=A, lifted=lifted, main=main, projection=projection,
                      projection_input_masks=V, ring=ring, direct_cse=direct,
                      measurements=measurements, best_completed_direct=baseline,
                      percent=100 * (baseline-ring['cost']) / baseline,
                      checked_nonzero_minors=minors,
                      xsets_source_sha256=sha256((Path(__file__).resolve().parents[1] / 'X-Sets.cpp').read_bytes()).hexdigest())
        results.append(result)
        args.output.write_text(json.dumps(results, indent=2) + '\n')
        print(json.dumps(dict(case=case, ring=ring['cost'], baseline=baseline,
                              incomplete=[dict(path=r['path'], method=r['method'])
                                          for r in measurements if r['status'] == 'timeout'])), flush=True)
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--profile', help='original or an example name from the certificate')
    action.add_argument('--search', action='store_true')
    parser.add_argument('--path', choices=['direct', 'main', 'projection'], default='direct')
    parser.add_argument('--method', choices=['MW_SQ', 'UBER_XSET'], default='MW_SQ')
    parser.add_argument('--params', type=int, nargs=3, default=[0, 3, 0])
    parser.add_argument('--timeout', type=float, default=180)
    parser.add_argument('--limit', type=int, default=10)
    parser.add_argument('--family', choices=['vandermonde', 'all'], default='vandermonde',
                        help='search unscaled Vandermonde matrices by default')
    parser.add_argument('--samples', type=int, default=1000)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.timeout <= 0 or args.limit < 1 or args.samples < 0:
        parser.error('timeout and limit must be positive; samples must be nonnegative')
    data = json.loads(CERTIFICATE.read_text())
    if not args.profile and not args.search:
        print(json.dumps(verify(data, args.samples), indent=2))
    else:
        if args.output is None:
            parser.error('--profile and --search require --output to preserve results')
        source = Path(__file__).resolve().parents[1] / 'X-Sets.cpp'
        with tempfile.TemporaryDirectory(prefix='ring-family-') as temp:
            binary = Path(temp) / 'X-Sets'
            subprocess.run(['g++', '-std=c++17', '-O3', str(source), '-o', str(binary)], check=True)
            if args.search:
                search(data, args, binary)
            else:
                result = profile(data, args, binary)
                result['xsets_source_sha256'] = sha256(source.read_bytes()).hexdigest()
                args.output.write_text(json.dumps(result, indent=2) + '\n')
                print(json.dumps({key: value for key, value in result.items() if key != 'circuit'}, indent=2))
