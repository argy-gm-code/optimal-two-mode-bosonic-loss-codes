#!/usr/bin/env python3
"""Evaluate saved encoders at a chosen eta by optimizing only their decoders."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import pickle
import sys
import time
import warnings

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "code"))
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(name, "1")


def evaluate(path, eta, eps):
    import numpy as np
    from graph_n_photons_new import test_code

    data = path.read_bytes()
    n, saved_eta, d, rounds, seed, _, encoder, saved_fe, _ = pickle.loads(data)
    if encoder.shape != ((n + 1)**2, d**2):
        raise ValueError("Saved encoder dimensions disagree with n and logical dimension")
    started = time.monotonic()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        optimizer = test_code(path, eta, eps=eps, return_optimizer=True)
    if not np.array_equal(encoder, optimizer.T_E):
        raise RuntimeError("The saved encoder changed during decoder-only evaluation")
    infos = list(optimizer.decoder_solver_info.values())
    diagnostics = {
        "accepted": True,
        "max_primal_residual": max(v["primal_residual"] for v in infos),
        "max_weighted_dual_residual": max(v["weight"] * v["dual_residual"] for v in infos),
        "sum_weighted_fidelity_gaps": sum(v["weight"] * v["gap"] / d**2 for v in infos),
        "max_cp_error": max(v["cp_error"] for v in infos),
        "max_tp_error": max(v["tp_error"] for v in infos),
        "max_iterations": max(v["iterations"] for v in infos),
        "inaccurate_sectors": [q for q, p in optimizer.decoder_sector_problems.items()
                               if p.status != "optimal"],
        "warnings": sorted({str(w.message) for w in caught
                            if "Defaulting to the SCIPY backend" not in str(w.message)}),
    }
    return {
        "path": path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else str(path),
        "source_sha256": hashlib.sha256(data).hexdigest(),
        "n": int(n), "logical_dimension": int(d), "seed": None if seed is None else int(seed),
        "saved_eta": float(saved_eta), "saved_fidelity": float(saved_fe),
        "evaluation_eta": eta, "eps": eps,
        "fidelity": float(optimizer.Fe),
        "coherent_information": float(optimizer.I_c),
        "encoder_unchanged": True, "decoder_diagnostics": diagnostics,
        "elapsed_seconds": time.monotonic() - started,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("codes", nargs="+", help="Run/reference IDs, or pickle paths relative to the release root")
    parser.add_argument("--eta", required=True, type=float, help="Evaluation transmissivity")
    parser.add_argument("--eps", type=float, default=1e-9, help="SCS tolerance (default: 1e-9)")
    parser.add_argument("--output", type=Path, help="Save JSON here, relative to the working directory")
    args = parser.parse_args(argv)
    if not 0 < args.eta < 1 or not 0 < args.eps <= 1e-5:
        parser.error("Require 0 < eta < 1 and 0 < eps <= 1e-5")
    if args.output and args.output.suffix.lower() != ".json":
        parser.error("--output must be a JSON file; no result pickle is written")
    catalogue = json.loads((HERE / "runs.json").read_text())
    records = dict(catalogue["analytic_references"])
    for point in catalogue["parameter_points"].values():
        records.update(point["runs"])
    paths = [(ROOT / Path(records[x]["path"] if x in records else x).expanduser()).resolve()
             for x in args.codes]
    for path in paths:
        if not path.is_file():
            parser.error(f"Input does not exist: {path}")
    document = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "method": "graph_n_photons_new.test_code; fixed encoder, decoder-only optimization",
        "graph_sha256": hashlib.sha256((ROOT / "code/graph_n_photons_new.py").read_bytes()).hexdigest(),
        "optimizer_sha256": hashlib.sha256((ROOT / "code/faster_gen_dual_rail.py").read_bytes()).hexdigest(),
        "wrapper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "versions": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "cvxpy", "scs")},
        "settings": {"eps": args.eps, "target_gap": args.eps,
                     "target_residual": 10 * args.eps, "target_physical": 10 * args.eps,
                     "strict_solver_mode": True, "max_iters": 100000,
                     "max_solver_retries": 5, "decoder_workers": 1},
        "evaluations": [],
    }
    for selector, path in zip(args.codes, paths):
        print(f"Evaluating {selector} at eta={args.eta}, eps={args.eps:g}", file=sys.stderr, flush=True)
        result = evaluate(path, args.eta, args.eps)
        result["code_id"] = selector if selector in records else None
        document["evaluations"].append(result)
        print(f"  fidelity={result['fidelity']}; accepted={result['decoder_diagnostics']['accepted']}",
              file=sys.stderr, flush=True)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(document, indent=2, allow_nan=False) + "\n")
    if not args.output:
        print(json.dumps(document, indent=2, allow_nan=False))
    return 0 if all(x["decoder_diagnostics"]["accepted"] for x in document["evaluations"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
