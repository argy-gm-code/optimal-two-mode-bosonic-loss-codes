#!/usr/bin/env python3
"""eta090_f0: approximate structural family (0.95 connected links); per-run symmetry retained: C6: rho_1+rho_3.

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
BEST_RUN_ID = 'n6_eta090_s3'
CANONICAL_NUMERICAL = 'n6_eta090_s3'
OUTPUT_ROOT = HERE / "plots" / 'eta090_f0'
VIEW = {'rotate': False, 'theta': 0.0, 'phi': 0.0}
# Human-readable titles printed on the figures; K is code rank/dimension.
NUMERICAL_TITLES = {'projector_sphere': 'Numerical code: projector\n$n=6$, $K=2$',
 'projector_plane': 'Numerical code: projector\n$n=6$, $K=2$',
 'logical0_sphere': 'Numerical code: logical state 0\n$n=6$, $K=2$',
 'logical0_plane': 'Numerical code: logical state 0\n$n=6$, $K=2$',
 'logical1_sphere': 'Numerical code: logical state 1\n$n=6$, $K=2$',
 'logical1_plane': 'Numerical code: logical state 1\n$n=6$, $K=2$'}
ANALYTIC_TITLES = None
# Analytic path and canonical orientation are deliberately separate.
# Standard mode: Each selected entry is sent to plot_file_simple directly, exactly once.
ENTRIES = [{'target': 'n6_eta090_s3',
  'path': 'repro_sweep/n6/n6_eta0.90_d2_r21000_s3.pkl',
  'fidelity': 0.9818485660673357,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [-1.8177523024879223, 1.5833023413650287, -4.3190569220326385],
  'output': 's3'},
 {'target': 'n6_eta090_s2',
  'path': 'repro_sweep/n6/n6_eta0.90_d2_r21000_s2.pkl',
  'fidelity': 0.9818485657217935,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [2.7660759455565884, 2.8925866422377426, 0.6574857659172129],
  'output': 's2'},
 {'target': 'n6_eta090_s4',
  'path': 'repro_sweep/n6/n6_eta0.90_d2_r21000_s4.pkl',
  'fidelity': 0.9818485655865409,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [-0.24590647257162557, 2.63938870770084, 3.778469056844179],
  'output': 's4'},
 {'target': 'n6_eta090_s0',
  'path': 'repro_sweep/n6/n6_eta0.90_d2_r21000_s0.pkl',
  'fidelity': 0.981848565492681,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [-2.7179155119740988, 1.889223635217386, -1.3665367339775787],
  'output': 's0'},
 {'target': 'n6_eta090_s1',
  'path': 'repro_sweep/n6/n6_eta0.90_d2_r21000_s1.pkl',
  'fidelity': 0.9818485646733509,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [0.25326523525958167, 2.081381685878401, -0.0004069353868163539],
  'output': 's1'}]
FAMILY_MINIMUM_OVERLAP = 0.9999999987566712
OVERLAPS_TO_BEST = {'n6_eta090_s3': 1.0, 'n6_eta090_s2': 0.999999999984848, 'n6_eta090_s4': 0.9999999994247182, 'n6_eta090_s0': 0.9999999999255151, 'n6_eta090_s1': 0.9999999992302244}
ANALYTIC_COMPARISONS = [{'path': 'analytic_codes/n6_C6_code.pkl',
  'classification': 'C6 reducible: rho_1 + rho_3',
  'closest_run_id': 'n6_eta090_s3',
  'overlap': 0.9324585280206621,
  'qualifies': False,
  'qualifying_run_ids': [],
  'best_run_overlap': 0.9324585280206621,
  'comparison_scripts': ['repro_analysis/n6/analytic_references/plot_n6_C6_code.py']}]


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
