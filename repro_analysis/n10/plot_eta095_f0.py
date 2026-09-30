#!/usr/bin/env python3
"""eta095_f0: approximate structural family (0.95 connected links); per-run symmetry retained: 2D4: rho_5; C4: rho_1+rho_2.

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
DEFAULT_ANALYTIC_FILE = 'analytic_codes/n10_2D4_code.pkl'
BEST_RUN_ID = 'n10_eta095_s1'
CANONICAL_NUMERICAL = 'n10_eta095_s1'
OUTPUT_ROOT = HERE / "plots" / 'eta095_f0'
ANALYTIC_OUTPUT = ROOT / 'repro_analysis/n10/analytic_references/n10_2D4_code'
VIEW = {'rotate': False, 'theta': 0.0, 'phi': 0.0}
# Human-readable titles printed on the figures; K is code rank/dimension.
NUMERICAL_TITLES = {'projector_sphere': 'Numerical code: projector\n$n=10$, $K=2$',
 'projector_plane': 'Numerical code: projector\n$n=10$, $K=2$',
 'logical0_sphere': 'Numerical code: logical state 0\n$n=10$, $K=2$',
 'logical0_plane': 'Numerical code: logical state 0\n$n=10$, $K=2$',
 'logical1_sphere': 'Numerical code: logical state 1\n$n=10$, $K=2$',
 'logical1_plane': 'Numerical code: logical state 1\n$n=10$, $K=2$'}
ANALYTIC_TITLES = {'projector_sphere': 'Constructed $2D_{4}$ code: projector\n$n=10$, $K=2$',
 'projector_plane': 'Constructed $2D_{4}$ code: projector\n$n=10$, $K=2$',
 'logical0_sphere': 'Constructed $2D_{4}$ code: logical state 0\n$n=10$, $K=2$',
 'logical0_plane': 'Constructed $2D_{4}$ code: logical state 0\n$n=10$, $K=2$',
 'logical1_sphere': 'Constructed $2D_{4}$ code: logical state 1\n$n=10$, $K=2$',
 'logical1_plane': 'Constructed $2D_{4}$ code: logical state 1\n$n=10$, $K=2$'}
# Analytic path and canonical orientation are deliberately separate.
# Standard mode: Each selected entry is sent to plot_file_simple directly, exactly once.
ENTRIES = [{'target': 'analytic',
  'path': 'analytic_codes/n10_2D4_code.pkl',
  'fidelity': 0.9943430410899461,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [-2.561594690529674, 1.0333077974535521e-09, 2.561594688916384],
  'output': 'analytic'},
 {'target': 'n10_eta095_s1',
  'path': 'repro_sweep/n10/n10_eta0.95_d2_r21000_s1.pkl',
  'fidelity': 0.9994063379625839,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [2.51517379524525, 1.7869913198779364, -3.260677584807304],
  'output': 's1'},
 {'target': 'n10_eta095_s0',
  'path': 'repro_sweep/n10/n10_eta0.95_d2_r21000_s0.pkl',
  'fidelity': 0.9994063176237471,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [4.349692199827732, 2.312643948259559, -0.4184243014455036],
  'output': 's0'},
 {'target': 'n10_eta095_s2',
  'path': 'repro_sweep/n10/n10_eta0.95_d2_r21000_s2.pkl',
  'fidelity': 0.9994063174601395,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [2.411156297018863, 2.491882438963763, -0.3601084233573699],
  'output': 's2'},
 {'target': 'n10_eta095_s3',
  'path': 'repro_sweep/n10/n10_eta0.95_d2_r21000_s3.pkl',
  'fidelity': 0.9994058923647077,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [-4.02602316771682, 2.3181612652967165, 0.8322455059475438],
  'output': 's3'},
 {'target': 'n10_eta095_s4',
  'path': 'repro_sweep/n10/n10_eta0.95_d2_r21000_s4.pkl',
  'fidelity': 0.9994012501897215,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [-0.3020384869051286, 2.2396747223347515, 0.5981358583948468],
  'output': 's4'}]
FAMILY_MINIMUM_OVERLAP = 0.9478335210640652
OVERLAPS_TO_BEST = {'n10_eta095_s1': 1.0, 'n10_eta095_s0': 0.9998988845376882, 'n10_eta095_s2': 0.9998987367170971, 'n10_eta095_s3': 0.9978723785850211, 'n10_eta095_s4': 0.9518720698332312}
ANALYTIC_COMPARISONS = [{'path': 'analytic_codes/n10_2D4_code.pkl',
  'classification': '2D4 pure irrep',
  'closest_run_id': 'n10_eta095_s2',
  'overlap': 0.9999748661521702,
  'qualifies': True,
  'qualifying_run_ids': ['n10_eta095_s1', 'n10_eta095_s0', 'n10_eta095_s2', 'n10_eta095_s3'],
  'best_run_overlap': 0.9998750852337654,
  'comparison_scripts': ['repro_analysis/n10/analytic_references/plot_n10_2D4_code.py']}]


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
            pkl_path=str(source), out_folder=str(ANALYTIC_OUTPUT if entry["target"] == "analytic" else OUTPUT_ROOT / entry["output"]),
            interactive=args.interactive, reflect=entry["reflect"],
            projector_only=args.projector_only,
            alignment_euler=entry["alignment_euler"],
            rotate=VIEW["rotate"], theta=VIEW["theta"], phi=VIEW["phi"],
            titles=ANALYTIC_TITLES if entry["target"] == "analytic" else NUMERICAL_TITLES,
        )
        graph_n_photons_new.plot_file_simple(**kwargs)


if __name__ == "__main__":
    main()
