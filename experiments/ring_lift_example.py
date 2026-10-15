"""Reproduce and independently verify one ring-lift XOR circuit. Standard library only.

Run with --exhaustive to check all 65,536 data vectors and the MDS distance.
The direct baselines are finite heuristic results, not optimality claims.
"""

import argparse
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path
import random


def polynomial_product(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        b >>= 1
    return result


def remainder(a, modulus):
    while a and a.bit_length() >= modulus.bit_length():
        a ^= modulus << (a.bit_length() - modulus.bit_length())
    return a


def field_product(a, b):
    """Independent shift-and-reduce arithmetic in F16, not product-then-remainder."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 31
    return result


def cyclic_product(a, b):
    """Independent multiplication in R by XORing cyclic rotations of b."""
    result = 0
    for shift in range(5):
        if (a >> shift) & 1:
            result ^= ((b << shift) | (b >> (5 - shift))) & 31
    return result


def project(value):
    return (value & 15) ^ (15 if value & 16 else 0)


def matrix_masks(matrix, out_width, multiply):
    rows = []
    for row in matrix:
        columns = [multiply(a, 1 << bit) for a in row for bit in range(4)]
        rows.extend(sum(((column >> bit) & 1) << j
                        for j, column in enumerate(columns))
                    for bit in range(out_width))
    return rows


class GreedyCircuit:
    """Pair-frequency CSE, with reuse of already available linear forms.

    Seed zero uses lexicographic ties. Other seeds randomize tied pairs and
    the final fold order. Both direct and ring reconstructions use this code.
    """

    def __init__(self):
        self.masks = [1 << i for i in range(16)] + [0]
        self.by_mask = {mask: i for i, mask in enumerate(self.masks)}
        self.ops = {}

    def xor(self, a, b):
        mask = self.masks[a] ^ self.masks[b]
        if mask in self.by_mask:
            return self.by_mask[mask]
        node = len(self.masks)
        self.masks.append(mask)
        self.by_mask[mask] = node
        self.ops[node] = (a, b)
        return node

    def synthesize(self, expressions, seed=0):
        expressions = [set(expression) - {16} for expression in expressions]
        rng = random.Random(seed)
        while True:
            frequencies = Counter(pair for expression in expressions
                                  for pair in combinations(sorted(expression), 2))
            frequency = max(frequencies.values(), default=0)
            if frequency < 2:
                break
            choices = sorted(pair for pair, count in frequencies.items()
                             if count == frequency)
            a, b = choices[0] if seed == 0 else rng.choice(choices)
            node = self.xor(a, b)
            for expression in expressions:
                if a in expression and b in expression:
                    expression.remove(a)
                    expression.remove(b)
                    if node != 16:
                        expression.symmetric_difference_update({node})
        outputs = []
        for expression in expressions:
            mask = 0
            for i in expression:
                mask ^= self.masks[i]
            if mask in self.by_mask:
                outputs.append(self.by_mask[mask])
                continue
            order = sorted(expression)
            if seed:
                rng.shuffle(order)
            node = 16
            for i in order:
                node = self.xor(node, i)
            outputs.append(node)
        return outputs

    def synthesize_masks(self, targets, seed=0):
        return self.synthesize([{j for j in range(16) if (mask >> j) & 1}
                                for mask in targets], seed)

    def export(self, outputs):
        live = set()

        def visit(node):
            if node in live or node not in self.ops:
                return
            live.add(node)
            for parent in self.ops[node]:
                visit(parent)

        for node in outputs:
            visit(node)
        return dict(inputs=16, zero_node=16,
                    ops=[[node, *self.ops[node]] for node in sorted(live)],
                    outputs=outputs, cost=len(live))


def verify_circuit(circuit, targets):
    """Verify every gate, every final output, and that no counted gate is dead."""
    values = {j: 1 << j for j in range(16)}
    if circuit['zero_node'] is not None:
        assert circuit['zero_node'] not in values
        values[circuit['zero_node']] = 0
    parents = {}
    for node, a, b in circuit['ops']:
        assert node not in values and a in values and b in values
        values[node] = values[a] ^ values[b]
        parents[node] = (a, b)
    assert [values[node] for node in circuit['outputs']] == targets
    assert len(circuit['ops']) == circuit['cost']
    live = set()

    def visit(node):
        if node in live or node not in parents:
            return
        live.add(node)
        for parent in parents[node]:
            visit(parent)

    for node in circuit['outputs']:
        visit(node)
    assert len(live) == len(parents)
    return values


def evaluate(circuit, data):
    values = [(data >> j) & 1 for j in range(16)] + [0] * (
        max(node for node, _, _ in circuit['ops']) - 15)
    for node, a, b in circuit['ops']:
        values[node] = values[a] ^ values[b]
    return [sum(values[circuit['outputs'][4 * i + bit]] << bit for bit in range(4))
            for i in range(2)]


def verify(exhaustive):
    certificate = json.loads(Path(__file__).with_suffix('.json').read_text())
    matrix = certificate['A']
    lifted = certificate['ring']['lifted']
    assert certificate['w'] == 4 and certificate['P'] == 31 and certificate['Q'] == 3
    assert polynomial_product(31, 3) == 33
    assert all(remainder(31, factor) for factor in range(2, 8))
    assert all(min(a.bit_count(), (a ^ 31).bit_count()) <= 2 for a in range(16))
    assert [[project(a) for a in row] for row in lifted] == matrix
    assert all(a.bit_count() <= 2 for row in lifted for a in row)
    points = certificate['points']
    assert len(points) == len(set(points)) == 6
    assert all(field_product(matrix[i][j], points[i] ^ points[2 + j]) == 1
               for i in range(2) for j in range(4))
    minors = [field_product(matrix[0][i], matrix[1][j]) ^
              field_product(matrix[0][j], matrix[1][i]) for i, j in combinations(range(4), 2)]
    assert all(minors) and all(a != 0 for row in matrix for a in row)

    targets = matrix_masks(matrix, 4, field_product)
    middle = matrix_masks(lifted, 5, cyclic_product)
    projected = [middle[5 * i + bit] ^ middle[5 * i + 4]
                 for i in range(2) for bit in range(4)]
    assert projected == targets
    for name in ('direct', 'direct_default', 'ring'):
        values = verify_circuit(certificate[name], targets)
        if name == 'ring':
            assert [values[node] for node in certificate[name]['middle_outputs']] == middle
    for circuit in certificate['uber_certificates'].values():
        verify_circuit(circuit, targets)

    histogram = Counter()
    for seed in range(100):
        c = GreedyCircuit()
        outputs = c.synthesize_masks(targets, seed)
        circuit = c.export(outputs)
        verify_circuit(circuit, targets)
        histogram[circuit['cost']] += 1
        if seed == certificate['direct_seed']:
            assert circuit['ops'] == certificate['direct']['ops']
            assert circuit['outputs'] == certificate['direct']['outputs']
    assert {str(cost): count for cost, count in histogram.items()} == certificate['direct_histogram']
    assert min(histogram) == 36

    c = GreedyCircuit()
    z = c.synthesize_masks(middle)
    main_count = len(c.ops)
    outputs = c.synthesize([{z[5 * i + bit], z[5 * i + 4]}
                            for i in range(2) for bit in range(4)])
    ring = c.export(outputs)
    assert ring['ops'] == certificate['ring']['ops']
    assert ring['outputs'] == certificate['ring']['outputs']
    assert main_count == certificate['ring']['main_cost'] == 27
    assert len(c.ops) - main_count == certificate['ring']['projection_cost'] == 8
    assert ring['cost'] == 35

    report = dict(matrix_identity=True, all_saved_circuits_verified=True,
                  cauchy_matrix=True, nonzero_2x2_minors=minors,
                  direct_cse_default=certificate['direct_default']['cost'],
                  direct_cse_100_seeds_best=min(histogram),
                  direct_cse_histogram=dict(sorted(histogram.items())),
                  uber_costs={level: circuit['cost']
                              for level, circuit in certificate['uber_certificates'].items()},
                  ring_main=main_count, ring_projection=8, ring_total=ring['cost'])
    if exhaustive:
        count = 0
        min_distance = 7
        for data in product(range(16), repeat=4):
            wanted = [field_product(row[0], data[0]) ^ field_product(row[1], data[1]) ^
                      field_product(row[2], data[2]) ^ field_product(row[3], data[3])
                      for row in matrix]
            word = sum(value << (4 * j) for j, value in enumerate(data))
            assert evaluate(certificate['ring'], word) == wanted
            assert evaluate(certificate['direct'], word) == wanted
            for i in range(2):
                value = 0
                for j in range(4):
                    value ^= cyclic_product(lifted[i][j], data[j])
                assert project(value) == wanted[i]
            if any(data):
                min_distance = min(min_distance, sum(v != 0 for v in (*data, *wanted)))
            count += 1
        assert count == 65536 and min_distance == 3
        report.update(exhaustive_inputs=count, systematic_code=[6, 4, min_distance])
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exhaustive', action='store_true')
    verify(parser.parse_args().exhaustive)
