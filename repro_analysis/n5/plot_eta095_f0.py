#!/usr/bin/env python3
"""eta095_f0: approximate structural family (0.95 connected links); per-run symmetry retained: C6: rho_2+rho_4; C6: rho_2+rho_6.

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
BEST_RUN_ID = 'n5_eta095_s3'
CANONICAL_NUMERICAL = 'n5_eta095_s3'
OUTPUT_ROOT = HERE / "plots" / 'eta095_f0'
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
 {'target': 'n5_eta095_s3',
  'path': 'repro_sweep/n5/n5_eta0.95_d2_r21000_s3.pkl',
  'fidelity': 0.9947440593269161,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [1.378908400314614, 0.523083688142802, -1.2751370186544113],
  'output': 's3'},
 {'target': 'n5_eta095_s1',
  'path': 'repro_sweep/n5/n5_eta0.95_d2_r21000_s1.pkl',
  'fidelity': 0.994744058868821,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [4.135133389710689, 1.5655086020069664, -1.2032199937419004],
  'output': 's1'},
 {'target': 'n5_eta095_s4',
  'path': 'repro_sweep/n5/n5_eta0.95_d2_r21000_s4.pkl',
  'fidelity': 0.9947440583947265,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [-5.264838056162205, 1.52922934582875, 0.7490033879957614],
  'output': 's4'},
 {'target': 'n5_eta095_s2',
  'path': 'repro_sweep/n5/n5_eta0.95_d2_r21000_s2.pkl',
  'fidelity': 0.9947440570200283,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [-0.9214238705854233, 0.8636769782140054, 2.6744846533363016],
  'output': 's2'},
 {'target': 'n5_eta095_s0',
  'path': 'repro_sweep/n5/n5_eta0.95_d2_r21000_s0.pkl',
  'fidelity': 0.9947440496336252,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [4.766932727867086, 1.7963325713533471, -0.022007201775163132],
  'output': 's0'}]
FAMILY_MINIMUM_OVERLAP = 0.9999999965529522
OVERLAPS_TO_BEST = {'n5_eta095_s3': 1.0, 'n5_eta095_s1': 0.999999997384997, 'n5_eta095_s4': 0.9999999999261016, 'n5_eta095_s2': 0.9999999992301828, 'n5_eta095_s0': 0.9999999990369272}
ANALYTIC_COMPARISONS = [{'path': 'analytic_codes/n5_C6_code_26.pkl',
  'classification': 'C6 reducible: rho_2 + rho_6',
  'closest_run_id': 'n5_eta095_s2',
  'overlap': 0.9994908391022206,
  'qualifies': True,
  'qualifying_run_ids': ['n5_eta095_s3',
                         'n5_eta095_s1',
                         'n5_eta095_s4',
                         'n5_eta095_s2',
                         'n5_eta095_s0'],
  'best_run_overlap': 0.9994907840886054,
  'comparison_scripts': ['repro_analysis/n5/analytic_references/plot_n5_C6_code_26.py']},
 {'path': 'analytic_codes/n5_2D4_code.pkl',
  'classification': '2D4: rho_7 (pure irrep)',
  'closest_run_id': 'n5_eta095_s1',
  'overlap': 0.7957439650957236,
  'qualifies': False,
  'qualifying_run_ids': [],
  'best_run_overlap': 0.7957311278060197,
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
