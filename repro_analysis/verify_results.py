#!/usr/bin/env python3
"""Reproduce saved fidelities, overlaps, families, distances, and certificates."""
import argparse
from collections import defaultdict
import gc
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'code'))
sys.path.insert(0, str(ROOT))
os.environ.setdefault('MPLBACKEND', 'Agg')
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(name, '1')


def read(name):
    return json.loads((ROOT / 'repro_analysis' / name).read_text())


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def check_metadata():
    catalogue = read('runs.json')
    records = dict(catalogue['analytic_references'])
    for point in catalogue['parameter_points'].values():
        records.update(point['runs'])
        best = max(point['runs'], key=lambda rid: point['runs'][rid]['fidelity'])
        require(best == point['best_in_sweep'], 'Sweep winner mismatch')
    for rid, record in records.items():
        path = ROOT / record['path']
        data = pickle.loads(path.read_bytes())
        require(data[7] == record['fidelity'], f'{rid}: saved fidelity mismatch')
        if 'sha256' in record:
            require(hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256'], f'{rid}: hash mismatch')
        if 'standalone_plot_script' in record:
            require((ROOT / record['standalone_plot_script']).is_file(), f'{rid}: missing plotting script')
    evaluated = 0
    for path in [ROOT / 'repro_analysis/performance.json', ROOT / 'literature_codes/comparison.json']:
        performance = json.loads(path.read_text())
        for key, relative in [('graph_sha256', 'code/graph_n_photons_new.py'),
                              ('optimizer_sha256', 'code/faster_gen_dual_rail.py'),
                              ('wrapper_sha256', 'repro_analysis/evaluate_performance.py')]:
            require(hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == performance[key], relative)
        for item in performance['evaluations']:
            require(hashlib.sha256((ROOT / item['path']).read_bytes()).hexdigest() == item['source_sha256'], item['path'])
            require(item['decoder_diagnostics']['accepted'] and item['encoder_unchanged'], item['path'])
            evaluated += 1
    print(f'PASS: {len(records)} catalogue records; {evaluated} performance inputs', flush=True)


def check_saved():
    import numpy as np
    from faster_gen_dual_rail import TruncatedLoss
    groups = defaultdict(list)
    data_folders = ('repro_sweep', 'analytic_codes', 'additional_codes',
                    'literature_codes', 'repro_analysis')
    paths = [path for folder in data_folders
             for path in (ROOT / folder).rglob('*.pkl')]
    for path in sorted(paths):
        data = pickle.loads(path.read_bytes())
        require(len(data) == 9, str(path))
        groups[(data[0], data[1])].append(path)
    maximum = 0.0
    count = 0
    for (n, eta), paths in sorted(groups.items()):
        loss = TruncatedLoss(n, eta)
        for path in paths:
            data = pickle.loads(path.read_bytes())
            dimension, decoder, encoder = data[2], data[5], data[6]
            actual = float(np.trace(decoder @ (loss.superoperator @ encoder)).real / dimension**2)
            error = abs(actual - data[7])
            require(np.isfinite(actual) and error < 1e-9, f'{path.name}: fidelity error {error}')
            maximum = max(maximum, error)
            count += 1
        print(f'n={n}, eta={eta}: {count} saved results checked', flush=True)
        del loss
        gc.collect()
    print(f'PASS: {count} saved fidelities; maximum difference {maximum:.3g}', flush=True)


def check_overlaps_and_families():
    import graph_n_photons_new as graph
    catalogue = read('runs.json')
    comparisons = read('comparisons.json')['parameter_points']
    families = read('families.json')['parameter_points']
    records = dict(catalogue['analytic_references'])
    for point in catalogue['parameter_points'].values():
        records.update(point['runs'])
    maximum = 0.0
    count = family_count = 0
    for pid, pairs in comparisons.items():
        ids = set(catalogue['parameter_points'][pid]['runs'])
        parent = {rid: rid for rid in ids}
        def find(rid):
            while parent[rid] != rid:
                rid = parent[rid]
            return rid
        for pair in pairs:
            a, b = pair['source'], pair['target']
            actual = graph.projector_overlap_from_files(ROOT / records[a]['path'], ROOT / records[b]['path'], **pair['rotation'])
            error = abs(actual - pair['overlap'])
            require(error < 1e-8, f'{a}, {b}: overlap error {error}')
            maximum = max(maximum, error)
            count += 1
            if a in ids and b in ids and pair['overlap'] >= 0.95:
                parent[find(a)] = find(b)
        components = defaultdict(set)
        for rid in ids:
            components[find(rid)].add(rid)
        saved = families[pid]
        require({frozenset(c) for c in components.values()} == {frozenset(f['members']) for f in saved.values()}, f'{pid}: families differ')
        for family in saved.values():
            family_count += 1
            require(family['best'] == max(family['members'], key=lambda rid: records[rid]['fidelity']), f'{pid}: family winner')
            require((ROOT / family['plot_script']).is_file(), f'{pid}: missing plot script')
            reference = family['analytic_reference']
            if reference:
                require(any({p['source'], p['target']} == {rid, reference} and p['overlap'] >= 0.95
                            for p in pairs for rid in family['members']), f'{pid}: analytic reference')
        print(f'{pid}: overlaps and families checked', flush=True)
    print(f'PASS: {count} overlaps; {family_count} families; maximum difference {maximum:.3g}', flush=True)


def check_distances():
    import numpy as np
    import graph_n_photons_new as graph
    from group_rep_codes_classes import SpinSector
    document = read('constructions.json')
    tolerance = document['distance_tolerance']
    for entry in document['codes']:
        path = ROOT / entry['path']
        require(hashlib.sha256(path.read_bytes()).hexdigest() == entry['sha256'], entry['path'])
        raw, metadata = graph.spin_projector_from_result_file(path)
        basis = graph.projector_basis((raw + raw.conj().T) / 2, rank=metadata['d'])
        spin = SpinSector(metadata['n'])
        residuals = []
        for rank in range(entry['distance'] + 1):
            errors = []
            for tensor in spin.irreducible_tensors(rank).values():
                compressed = basis.conj().T @ tensor @ basis
                errors.append(float(np.linalg.norm(compressed - np.trace(compressed) * np.eye(metadata['d']) / metadata['d'])))
            residuals.append(max(errors))
        residual = max(residuals[:-1])
        require(residual < tolerance and residuals[-1] > tolerance, f'{path.name}: distance mismatch {residuals}')
        print(f'{path.name}: distance {entry["distance"]}, KL residual {residual:.3g}, next-rank failure {residuals[-1]:.3g}', flush=True)
    print(f'PASS: {len(document["codes"])} constructed-code distances at tolerance {tolerance:g}', flush=True)


def check_certificates():
    import sympy as sp
    from Farkas_certificates import generate_farkas_system, verify_farkas_certificate
    paths = sorted((ROOT / 'Distance Certificates').glob('*.json'))
    for path in paths:
        document = json.loads(path.read_text())
        parameters = document['parameters']
        matrix, objective = generate_farkas_system(parameters['n'], parameters['K'], parameters['d'])
        certificate = sp.Matrix([sp.Rational(x) for x in document['certificate']])
        require(verify_farkas_certificate(matrix, objective, certificate), path.name)
        print(f'{path.name}: exact verification passed', flush=True)
    print(f'PASS: {len(paths)} exact certificates', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choices = ['metadata', 'saved', 'overlaps', 'distances', 'certificates']
    parser.add_argument('checks', nargs='+', choices=choices + ['all'])
    args = parser.parse_args()
    checks = choices if 'all' in args.checks else args.checks
    functions = {'metadata': check_metadata, 'saved': check_saved, 'overlaps': check_overlaps_and_families,
                 'distances': check_distances, 'certificates': check_certificates}
    for check in checks:
        functions[check]()


if __name__ == '__main__':
    main()
