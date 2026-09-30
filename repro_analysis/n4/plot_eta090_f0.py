#!/usr/bin/env python3
"""eta090_f0: approximate structural family (0.95 connected links); per-run symmetry retained: 2O: rho_3.

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
DEFAULT_ANALYTIC_FILE = 'analytic_codes/n4_2O_code.pkl'
BEST_RUN_ID = 'n4_eta090_s1'
CANONICAL_NUMERICAL = 'n4_eta090_s1'
OUTPUT_ROOT = HERE / "plots" / 'eta090_f0'
ANALYTIC_OUTPUT = ROOT / 'repro_analysis/n4/analytic_references/n4_2O_code'
VIEW = {'rotate': True, 'theta': 0.7853981633974477, 'phi': 0.9553166181245097}
# Human-readable titles printed on the figures; K is code rank/dimension.
NUMERICAL_TITLES = {'projector_sphere': 'Numerical code: projector\n$n=4$, $K=2$',
 'projector_plane': 'Numerical code: projector\n$n=4$, $K=2$',
 'logical0_sphere': 'Numerical code: logical state 0\n$n=4$, $K=2$',
 'logical0_plane': 'Numerical code: logical state 0\n$n=4$, $K=2$',
 'logical1_sphere': 'Numerical code: logical state 1\n$n=4$, $K=2$',
 'logical1_plane': 'Numerical code: logical state 1\n$n=4$, $K=2$'}
ANALYTIC_TITLES = {'projector_sphere': 'Constructed $2O$ code: projector\n$n=4$, $K=2$',
 'projector_plane': 'Constructed $2O$ code: projector\n$n=4$, $K=2$',
 'logical0_sphere': 'Constructed $2O$ code: logical state 0\n$n=4$, $K=2$',
 'logical0_plane': 'Constructed $2O$ code: logical state 0\n$n=4$, $K=2$',
 'logical1_sphere': 'Constructed $2O$ code: logical state 1\n$n=4$, $K=2$',
 'logical1_plane': 'Constructed $2O$ code: logical state 1\n$n=4$, $K=2$'}
# Analytic path and canonical orientation are deliberately separate.
# Standard mode: Each selected entry is sent to plot_file_simple directly, exactly once.
ENTRIES = [{'target': 'analytic',
  'path': 'analytic_codes/n4_2O_code.pkl',
  'fidelity': 0.972925000001253,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [-1.8309105932142482, 6.234081354616416e-09, 1.830910586834444],
  'output': 'analytic'},
 {'target': 'n4_eta090_s1',
  'path': 'repro_sweep/n4/n4_eta0.90_d2_r21000_s1.pkl',
  'fidelity': 0.9732733654475184,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [1.8227151531774575, 0.5548776266536745, -0.6118246646061233],
  'output': 's1'},
 {'target': 'n4_eta090_s4',
  'path': 'repro_sweep/n4/n4_eta0.90_d2_r21000_s4.pkl',
  'fidelity': 0.9732733652307413,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [-3.866066836269873, 1.1430757793774386, 0.553487200054487],
  'output': 's4'},
 {'target': 'n4_eta090_s2',
  'path': 'repro_sweep/n4/n4_eta0.90_d2_r21000_s2.pkl',
  'fidelity': 0.9732733648201679,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [3.4137928139732505, 0.8696523648813651, -2.647489914811591],
  'output': 's2'},
 {'target': 'n4_eta090_s3',
  'path': 'repro_sweep/n4/n4_eta0.90_d2_r21000_s3.pkl',
  'fidelity': 0.9732733646605283,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [1.1538095325928361, 1.577234933953301, -2.3626671663096825],
  'output': 's3'},
 {'target': 'n4_eta090_s0',
  'path': 'repro_sweep/n4/n4_eta0.90_d2_r21000_s0.pkl',
  'fidelity': 0.9732733629384898,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [2.211567153012107, 1.787313017043448, -1.8255877464820611],
  'output': 's0'}]
FAMILY_MINIMUM_OVERLAP = 0.9999999984230661
OVERLAPS_TO_BEST = {'n4_eta090_s1': 1.0, 'n4_eta090_s4': 0.9999999993288732, 'n4_eta090_s2': 0.9999999991744261, 'n4_eta090_s3': 0.999999999044409, 'n4_eta090_s0': 0.9999999992257264}
ANALYTIC_COMPARISONS = [{'path': 'analytic_codes/n4_2O_code.pkl',
  'classification': '2O pure irrep (CLY code)',
  'closest_run_id': 'n4_eta090_s4',
  'overlap': 0.9940732327784625,
  'qualifies': True,
  'qualifying_run_ids': ['n4_eta090_s1',
                         'n4_eta090_s4',
                         'n4_eta090_s2',
                         'n4_eta090_s3',
                         'n4_eta090_s0'],
  'best_run_overlap': 0.9940732229898555,
  'comparison_scripts': ['repro_analysis/n4/analytic_references/plot_n4_2O_code.py']}]


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
