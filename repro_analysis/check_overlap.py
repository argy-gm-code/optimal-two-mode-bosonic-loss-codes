#!/usr/bin/env python3
"""Recheck a saved comparison using runs.json and comparisons.json; no alignment search."""
import argparse
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("code_a", help="Run/reference ID or recorded pickle path")
    parser.add_argument("code_b", help="Run/reference ID or recorded pickle path")
    args = parser.parse_args()

    runs = json.loads((HERE / "runs.json").read_text())
    comparisons = json.loads((HERE / "comparisons.json").read_text())
    records = dict(runs["analytic_references"])
    for point in runs["parameter_points"].values():
        records.update(point["runs"])

    selected = []
    for selector in (args.code_a, args.code_b):
        if selector in records:
            selected.append(selector)
            continue
        path = (ROOT / Path(selector).expanduser()).resolve()
        matches = [key for key, record in records.items()
                   if (ROOT / record["path"]).resolve() == path]
        if len(matches) != 1:
            parser.error(f"Unknown or ambiguous code: {selector}. Use an ID from runs.json or a recorded pickle path.")
        selected.append(matches[0])

    pairs = [pair for point_pairs in comparisons["parameter_points"].values()
             for pair in point_pairs
             if {pair["source"], pair["target"]} == set(selected)]
    if len(pairs) != 1:
        parser.error(f"Expected one stored comparison for {selected[0]} and {selected[1]}; found {len(pairs)}. No alignment search was run.")
    pair = pairs[0]
    # Keep the recorded order even if the user supplies the two IDs in reverse.
    file_a = ROOT / records[pair["source"]]["path"]
    file_b = ROOT / records[pair["target"]]["path"]
    for path in (file_a, file_b):
        if not path.is_file():
            parser.error(f"Pickle does not exist: {path}")

    sys.dont_write_bytecode = True
    sys.path.insert(0, str(ROOT / "code"))
    os.environ.setdefault("MPLBACKEND", "Agg")
    import graph_n_photons_new

    try:
        actual = graph_n_photons_new.projector_overlap_from_files(
            file_a, file_b, **pair["rotation"],
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    difference = abs(actual - pair["overlap"])
    print("Stored source (A):", pair["source"], "->", file_a)
    print("Stored target (B):", pair["target"], "->", file_b)
    print("Stored rotations:", json.dumps(pair["rotation"]))
    print(f"Saved overlap:      {pair['overlap']:.16g}")
    print(f"Recomputed overlap: {actual:.16g}")
    print(f"Absolute difference: {difference:.3g}")
    if not difference < 1e-8:
        parser.exit(status=1, message="Overlap differs from the saved value by at least 1e-8; check whether an input pickle changed.\n")


if __name__ == "__main__":
    main()
