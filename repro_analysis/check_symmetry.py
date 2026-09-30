#!/usr/bin/env python3
"""Replay stored best-run group fits with the bundled scientific routines."""
import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_id", nargs="?", help="Best-run ID or parameter-point ID")
    parser.add_argument("--all", action="store_true", help="Check all 76 best runs")
    args = parser.parse_args()
    if args.all == bool(args.run_id):
        parser.error("Specify a run/point ID or --all")
    runs = json.loads((HERE / "runs.json").read_text())
    selected = [(pid, p["best_in_sweep"], p["runs"][p["best_in_sweep"]])
                for pid, p in runs["parameter_points"].items()
                if args.all or args.run_id in (pid, p["best_in_sweep"])]
    if not selected:
        parser.error("No stored best-run symmetry fit for that ID")

    sys.dont_write_bytecode = True
    sys.path.insert(0, str(ROOT / "code"))
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "repro_analysis_matplotlib"))
    import numpy as np
    import graph_n_photons_new as graph

    for _, rid, record in selected:
        fit = record["symmetry_fit"]
        source = ROOT / record["path"]
        if hashlib.sha256(source.read_bytes()).hexdigest() != fit["source_sha256"]:
            parser.exit(1, f"{rid}: source hash changed\n")
        raw, metadata = graph.spin_projector_from_result_file(str(source))
        basis = graph.projector_basis((raw + raw.conj().T) / 2, rank=metadata["d"])
        projector = basis @ basis.conj().T
        if fit["reflect"]:
            projector = projector.conj()
        projector = graph.rotate_projector_euler(projector, *fit["alignment_euler"])
        model = graph.GroupQECModel(str(ROOT / "code" / fit["group"]), n=metadata["n"])
        invariant = graph.reconstruct_invariant_projector(projector, model, fit["allocation"])
        loss = graph.group_commutator_loss(projector, model.physical_representations)
        overlap = float(np.trace(projector @ invariant).real / metadata["d"])
        errors = (abs(loss - record["commutator_loss"]),
                  abs(overlap - record["invariant_overlap"]))
        print(f"{rid}: {record['symmetry']}; loss={loss:.16g}; overlap={overlap:.16g}; "
              f"max difference={max(errors):.3g}")
        if not max(errors) < 1e-10:
            parser.exit(1, f"{rid}: stored diagnostics differ by at least 1e-10\n")
    print(f"Verified {len(selected)} stored group fits; no alignment search or optimization.")


if __name__ == "__main__":
    main()
