#!/usr/bin/env python3
"""eta090_f0: approximate structural family (0.95 connected links); per-run symmetry retained: 2O: rho_5.

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
DEFAULT_ANALYTIC_FILE = 'analytic_codes/n13_2O_code.pkl'
BEST_RUN_ID = 'n13_eta090_s4'
CANONICAL_NUMERICAL = 'n13_eta090_s4'
OUTPUT_ROOT = HERE / "plots" / 'eta090_f0'
ANALYTIC_OUTPUT = ROOT / 'repro_analysis/n13/analytic_references/n13_2O_code'
VIEW = {'rotate': True, 'theta': 0.7853981633974477, 'phi': 0.9553166181245097}
# Human-readable titles printed on the figures; K is code rank/dimension.
NUMERICAL_TITLES = {'projector_sphere': 'Numerical code: projector\n$n=13$, $K=2$',
 'projector_plane': 'Numerical code: projector\n$n=13$, $K=2$',
 'logical0_sphere': 'Numerical code: logical state 0\n$n=13$, $K=2$',
 'logical0_plane': 'Numerical code: logical state 0\n$n=13$, $K=2$',
 'logical1_sphere': 'Numerical code: logical state 1\n$n=13$, $K=2$',
 'logical1_plane': 'Numerical code: logical state 1\n$n=13$, $K=2$'}
ANALYTIC_TITLES = {'projector_sphere': 'Constructed $2O$ code: projector\n$n=13$, $K=2$',
 'projector_plane': 'Constructed $2O$ code: projector\n$n=13$, $K=2$',
 'logical0_sphere': 'Constructed $2O$ code: logical state 0\n$n=13$, $K=2$',
 'logical0_plane': 'Constructed $2O$ code: logical state 0\n$n=13$, $K=2$',
 'logical1_sphere': 'Constructed $2O$ code: logical state 1\n$n=13$, $K=2$',
 'logical1_plane': 'Constructed $2O$ code: logical state 1\n$n=13$, $K=2$'}
# Analytic path and canonical orientation are deliberately separate.
# Standard mode: Each selected entry is sent to plot_file_simple directly, exactly once.
ENTRIES = [{'target': 'analytic',
  'path': 'analytic_codes/n13_2O_code.pkl',
  'fidelity': 0.9974664033578942,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [1.587107774809809, 4.412497388161908e-09, -1.5871077733070582],
  'output': 'analytic'},
 {'target': 'n13_eta090_s4',
  'path': 'repro_sweep/n13/n13_eta0.90_d2_r21000_s4.pkl',
  'fidelity': 0.9974689177719939,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [-0.19130907788220686, 2.478748508020803, 5.142475732766739],
  'output': 's4'},
 {'target': 'n13_eta090_s3',
  'path': 'repro_sweep/n13/n13_eta0.90_d2_r21000_s3.pkl',
  'fidelity': 0.9974689173422291,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [-4.622270810099245, 0.6308402117032378, 0.7612919437857573],
  'output': 's3'},
 {'target': 'n13_eta090_s2',
  'path': 'repro_sweep/n13/n13_eta0.90_d2_r21000_s2.pkl',
  'fidelity': 0.9974689173123074,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [4.404440887794168, 1.763984751743464, -0.7534620440319713],
  'output': 's2'},
 {'target': 'n13_eta090_s0',
  'path': 'repro_sweep/n13/n13_eta0.90_d2_r21000_s0.pkl',
  'fidelity': 0.9974689172794596,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [3.2667841656863743, 0.5332859292472879, 1.5492146204787813],
  'output': 's0'},
 {'target': 'n13_eta090_s1',
  'path': 'repro_sweep/n13/n13_eta0.90_d2_r21000_s1.pkl',
  'fidelity': 0.9974689171996644,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [0.00945558206388264, 2.9403304258874927, -2.8279410682236277],
  'output': 's1'}]
FAMILY_MINIMUM_OVERLAP = 0.9999998776662755
OVERLAPS_TO_BEST = {'n13_eta090_s4': 1.0, 'n13_eta090_s3': 0.9999999430526387, 'n13_eta090_s2': 0.9999998776662755, 'n13_eta090_s0': 0.9999999201990202, 'n13_eta090_s1': 0.9999999464944306}
ANALYTIC_COMPARISONS = [{'path': 'analytic_codes/n13_2O_code.pkl',
  'classification': '2O pure irrep',
  'closest_run_id': 'n13_eta090_s0',
  'overlap': 0.9998379891203326,
  'qualifies': True,
  'qualifying_run_ids': ['n13_eta090_s4',
                         'n13_eta090_s3',
                         'n13_eta090_s2',
                         'n13_eta090_s0',
                         'n13_eta090_s1'],
  'best_run_overlap': 0.9998378914084853,
  'comparison_scripts': ['repro_analysis/n13/analytic_references/plot_n13_2O_code.py']}]


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
