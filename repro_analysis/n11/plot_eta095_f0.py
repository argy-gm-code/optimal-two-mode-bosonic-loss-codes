#!/usr/bin/env python3
"""eta095_f0: approximate structural family (0.95 connected links); per-run symmetry retained: C4: rho_3+rho_4.

Reflection -> fixed Euler prealignment -> common visual theta/phi rotation.
Projector orientation is shared; the saved logical basis is not gauge-aligned.
"""
import argparse
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_ANALYTIC_FILE = None
BEST_RUN_ID = 'n11_eta095_s3'
CANONICAL_NUMERICAL = 'n11_eta095_s3'
OUTPUT_ROOT = HERE / "plots" / 'eta095_f0'
VIEW = {'rotate': False, 'theta': 0.0, 'phi': 0.0}
# Human-readable titles printed on the figures; K is code rank/dimension.
NUMERICAL_TITLES = {'projector_sphere': 'Numerical code: projector\n$n=11$, $K=2$',
 'projector_plane': 'Numerical code: projector\n$n=11$, $K=2$',
 'logical0_sphere': 'Numerical code: logical state 0\n$n=11$, $K=2$',
 'logical0_plane': 'Numerical code: logical state 0\n$n=11$, $K=2$',
 'logical1_sphere': 'Numerical code: logical state 1\n$n=11$, $K=2$',
 'logical1_plane': 'Numerical code: logical state 1\n$n=11$, $K=2$'}
ANALYTIC_TITLES = None
# Analytic path and canonical orientation are deliberately separate.
# Standard mode: Each selected entry is sent to plot_file_simple directly, exactly once.
ENTRIES = [{'target': 'n11_eta095_s3',
  'path': 'repro_sweep/n11/n11_eta0.95_d2_r21000_s3.pkl',
  'fidelity': 0.999470543029294,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [0.0, 0.4265087465181524, -2.824036727440917],
  'output': 's3'},
 {'target': 'n11_eta095_s2',
  'path': 'repro_sweep/n11/n11_eta0.95_d2_r21000_s2.pkl',
  'fidelity': 0.999470540664234,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [-0.47777802229621696, 2.511978060729082, 0.3728261559066148],
  'output': 's2'},
 {'target': 'n11_eta095_s1',
  'path': 'repro_sweep/n11/n11_eta0.95_d2_r21000_s1.pkl',
  'fidelity': 0.9994704402646677,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [3.0926102908961877, 1.993636722464162, 2.69886370754237],
  'output': 's1'}]
FAMILY_MINIMUM_OVERLAP = 0.998083182114191
OVERLAPS_TO_BEST = {'n11_eta095_s3': 1.0, 'n11_eta095_s2': 0.9999620581289044, 'n11_eta095_s1': 0.998083182114191}
ANALYTIC_COMPARISONS = [{'path': 'analytic_codes/n11_C4_code.pkl',
  'classification': 'C4 reducible: rho_3 + rho_4',
  'closest_run_id': 'n11_eta095_s1',
  'overlap': 0.8991367806570907,
  'qualifies': False,
  'qualifying_run_ids': [],
  'best_run_overlap': 0.8796633374550372,
  'comparison_scripts': ['repro_analysis/n11/analytic_references/plot_n11_C4_code.py']},
 {'path': 'analytic_codes/n11_2D2_code.pkl',
  'classification': '2D2: rho_5 (pure irrep)',
  'closest_run_id': 'n11_eta095_s1',
  'overlap': 0.9085978008827116,
  'qualifies': False,
  'qualifying_run_ids': [],
  'best_run_overlap': 0.8875832048932001,
  'comparison_scripts': ['repro_analysis/n11/analytic_references/plot_n11_2D2_code.py']},
 {'path': 'analytic_codes/n11_C8_code.pkl',
  'classification': 'C8 reducible: rho_5 + rho_6',
  'closest_run_id': 'n11_eta095_s1',
  'overlap': 0.7556528887320277,
  'qualifies': False,
  'qualifying_run_ids': [],
  'best_run_overlap': 0.748516857513688,
  'comparison_scripts': ['repro_analysis/n11/analytic_references/plot_n11_C8_code.py']}]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--target", default="best")
    parser.add_argument("--analytic-file", type=Path)
    parser.add_argument("--projector-only", action="store_true", help="Plot only the projector sphere and plane")
    parser.add_argument("--interactive", action="store_true", help="Also save interactive sphere and plane HTML")
    parser.add_argument("--no-show", action="store_true", help="Save interactive HTML without opening browser windows")
    args = parser.parse_args()
    if args.list:
        print("Best numerical member:", BEST_RUN_ID)
        print("Structural rule: connected links with overlap >= 0.95; endpoint overlaps may be lower")
        print("Minimum pairwise overlap:", FAMILY_MINIMUM_OVERLAP)
        print("Targets: analytic, best (default), numericals, all, or a run ID")
        if DEFAULT_ANALYTIC_FILE is None:
            print("No analytic counterpart meets the 0.95 criterion for this family.")
        for comparison in ANALYTIC_COMPARISONS:
            print("Analytic counterpart:" if comparison["qualifies"] else "Separate analytic comparison:",
                  comparison["path"], "closest_member=", comparison["closest_run_id"], "overlap=", comparison["overlap"],
                  "scripts=", comparison["comparison_scripts"])
        for entry in ENTRIES:
            path = (args.analytic_file if args.analytic_file else Path(DEFAULT_ANALYTIC_FILE)) if entry["target"] == "analytic" else Path(entry["path"])
            print(entry["target"], "Fe=", entry["fidelity"], "saved_eta=", entry["saved_eta"], "path=", path,
                  "overlap_to_best=", OVERLAPS_TO_BEST.get(entry["target"], "analytic comparison"))
        return
    target = BEST_RUN_ID if args.target == "best" else args.target
    if target == "analytic" and DEFAULT_ANALYTIC_FILE is None:
        parser.exit(message="No analytic counterpart meets the 0.95 criterion for this family; use --list for separate reference scripts.\n")
    selected = [entry for entry in ENTRIES if target == "all" or
                (target == "numericals" and entry["target"] != "analytic") or entry["target"] == target]
    if not selected:
        parser.error("unknown target; use --list")
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(ROOT / "code"))
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "repro_analysis_matplotlib"))
    import graph_n_photons_new
    if args.interactive and args.no_show:
        import plotly.io as pio
        pio.renderers.default = ""
    for entry in selected:
        source = (args.analytic_file if args.analytic_file else Path(DEFAULT_ANALYTIC_FILE)) if entry["target"] == "analytic" else Path(entry["path"])
        source = source if source.is_absolute() else ROOT / source
        if not source.is_file():
            parser.error(f"pickle does not exist: {source}")
        kwargs = dict(
            pkl_path=str(source), out_folder=str(OUTPUT_ROOT / entry["output"]),
            interactive=args.interactive, reflect=entry["reflect"],
            projector_only=args.projector_only,
            alignment_euler=entry["alignment_euler"],
            rotate=VIEW["rotate"], theta=VIEW["theta"], phi=VIEW["phi"],
            titles=ANALYTIC_TITLES if entry["target"] == "analytic" else NUMERICAL_TITLES,
        )
        graph_n_photons_new.plot_file_simple(**kwargs)


if __name__ == "__main__":
    main()
