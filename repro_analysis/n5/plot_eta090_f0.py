#!/usr/bin/env python3
"""eta090_f0: approximate structural family (0.95 connected links); per-run symmetry retained: C6: rho_2+rho_4; C6: rho_2+rho_6.

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
DEFAULT_ANALYTIC_FILE = 'analytic_codes/n5_C6_code_26.pkl'
BEST_RUN_ID = 'n5_eta090_s2'
CANONICAL_NUMERICAL = 'n5_eta090_s2'
OUTPUT_ROOT = HERE / "plots" / 'eta090_f0'
ANALYTIC_OUTPUT = ROOT / 'repro_analysis/n5/analytic_references/n5_C6_code_26'
VIEW = {'rotate': False, 'theta': 0.0, 'phi': 0.0}
# Human-readable titles printed on the figures; K is code rank/dimension.
NUMERICAL_TITLES = {'projector_sphere': 'Numerical code: projector\n$n=5$, $K=2$',
 'projector_plane': 'Numerical code: projector\n$n=5$, $K=2$',
 'logical0_sphere': 'Numerical code: logical state 0\n$n=5$, $K=2$',
 'logical0_plane': 'Numerical code: logical state 0\n$n=5$, $K=2$',
 'logical1_sphere': 'Numerical code: logical state 1\n$n=5$, $K=2$',
 'logical1_plane': 'Numerical code: logical state 1\n$n=5$, $K=2$'}
ANALYTIC_TITLES = {'projector_sphere': 'Constructed $C_{6}$ code: projector\n$n=5$, $K=2$',
 'projector_plane': 'Constructed $C_{6}$ code: projector\n$n=5$, $K=2$',
 'logical0_sphere': 'Constructed $C_{6}$ code: logical state 0\n$n=5$, $K=2$',
 'logical0_plane': 'Constructed $C_{6}$ code: logical state 0\n$n=5$, $K=2$',
 'logical1_sphere': 'Constructed $C_{6}$ code: logical state 1\n$n=5$, $K=2$',
 'logical1_plane': 'Constructed $C_{6}$ code: logical state 1\n$n=5$, $K=2$'}
# Analytic path and canonical orientation are deliberately separate.
# Standard mode: Each selected entry is sent to plot_file_simple directly, exactly once.
ENTRIES = [{'target': 'analytic',
  'path': 'analytic_codes/n5_C6_code_26.pkl',
  'fidelity': 0.9793261156895041,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [2.098662039606354, 4.1762347090148484e-16, -2.098662039606354],
  'output': 'analytic'},
 {'target': 'n5_eta090_s2',
  'path': 'repro_sweep/n5/n5_eta0.90_d2_r21000_s2.pkl',
  'fidelity': 0.9797096084651045,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [-0.9256640398928924, 0.8543492267844051, 2.6902262809408057],
  'output': 's2'},
 {'target': 'n5_eta090_s4',
  'path': 'repro_sweep/n5/n5_eta0.90_d2_r21000_s4.pkl',
  'fidelity': 0.9797096042261311,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [3.2884136952562626, 1.7292280255384864, -0.7836592902381001],
  'output': 's4'},
 {'target': 'n5_eta090_s3',
  'path': 'repro_sweep/n5/n5_eta0.90_d2_r21000_s3.pkl',
  'fidelity': 0.9797096026457848,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [-0.6695258985234536, 0.5331479506607998, 4.98668824749004],
  'output': 's3'},
 {'target': 'n5_eta090_s1',
  'path': 'repro_sweep/n5/n5_eta0.90_d2_r21000_s1.pkl',
  'fidelity': 0.97970960106887,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [-2.3446158019472367, 1.5609697992501959, 1.2044520953099573],
  'output': 's1'},
 {'target': 'n5_eta090_s0',
  'path': 'repro_sweep/n5/n5_eta0.90_d2_r21000_s0.pkl',
  'fidelity': 0.979709600608413,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [1.7964524712080556, 1.087891535227028, -0.6952558522945745],
  'output': 's0'}]
FAMILY_MINIMUM_OVERLAP = 0.999999998520227
OVERLAPS_TO_BEST = {'n5_eta090_s2': 1.0, 'n5_eta090_s4': 0.9999999989796088, 'n5_eta090_s3': 0.9999999999939653, 'n5_eta090_s1': 0.999999998520227, 'n5_eta090_s0': 0.9999999996571225}
ANALYTIC_COMPARISONS = [{'path': 'analytic_codes/n5_C6_code_26.pkl',
  'classification': 'C6 reducible: rho_2 + rho_6',
  'closest_run_id': 'n5_eta090_s4',
  'overlap': 0.9979261173658326,
  'qualifies': True,
  'qualifying_run_ids': ['n5_eta090_s2',
                         'n5_eta090_s4',
                         'n5_eta090_s3',
                         'n5_eta090_s1',
                         'n5_eta090_s0'],
  'best_run_overlap': 0.9979259182806917,
  'comparison_scripts': ['repro_analysis/n5/analytic_references/plot_n5_C6_code_26.py']},
 {'path': 'analytic_codes/n5_2D4_code.pkl',
  'classification': '2D4: rho_7 (pure irrep)',
  'closest_run_id': 'n5_eta090_s4',
  'overlap': 0.7936271446139135,
  'qualifies': False,
  'qualifying_run_ids': [],
  'best_run_overlap': 0.7936241153512786,
  'comparison_scripts': ['repro_analysis/n5/analytic_references/plot_n5_2D4_code.py']}]


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
