#!/usr/bin/env python3
"""eta090_f0: approximate structural family (0.95 connected links); per-run symmetry retained: 2D2: rho_2+rho_3; 2D2: rho_3+rho_4; C4: rho_2 x2.

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
DEFAULT_ANALYTIC_FILE = 'analytic_codes/n18_2D2_code.pkl'
BEST_RUN_ID = 'n18_eta090_s4'
CANONICAL_NUMERICAL = 'n18_eta090_s4'
OUTPUT_ROOT = HERE / "plots" / 'eta090_f0'
ANALYTIC_OUTPUT = ROOT / 'repro_analysis/n18/analytic_references/n18_2D2_code'
VIEW = {'rotate': False, 'theta': 0.0, 'phi': 0.0}
# Human-readable titles printed on the figures; K is code rank/dimension.
NUMERICAL_TITLES = {'projector_sphere': 'Numerical code: projector\n$n=18$, $K=2$',
 'projector_plane': 'Numerical code: projector\n$n=18$, $K=2$',
 'logical0_sphere': 'Numerical code: logical state 0\n$n=18$, $K=2$',
 'logical0_plane': 'Numerical code: logical state 0\n$n=18$, $K=2$',
 'logical1_sphere': 'Numerical code: logical state 1\n$n=18$, $K=2$',
 'logical1_plane': 'Numerical code: logical state 1\n$n=18$, $K=2$'}
ANALYTIC_TITLES = {'projector_sphere': 'Constructed $2D_{2}$ code: projector\n$n=18$, $K=2$',
 'projector_plane': 'Constructed $2D_{2}$ code: projector\n$n=18$, $K=2$',
 'logical0_sphere': 'Constructed $2D_{2}$ code: logical state 0\n$n=18$, $K=2$',
 'logical0_plane': 'Constructed $2D_{2}$ code: logical state 0\n$n=18$, $K=2$',
 'logical1_sphere': 'Constructed $2D_{2}$ code: logical state 1\n$n=18$, $K=2$',
 'logical1_plane': 'Constructed $2D_{2}$ code: logical state 1\n$n=18$, $K=2$'}
# Analytic path and canonical orientation are deliberately separate.
# Standard mode: Each selected entry is sent to plot_file_simple directly, exactly once.
ENTRIES = [{'target': 'analytic',
  'path': 'analytic_codes/n18_2D2_code.pkl',
  'fidelity': 0.9986794735106906,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [3.1415926464396375, 1.5707963197241004, 3.141592647123489],
  'output': 'analytic'},
 {'target': 'n18_eta090_s4',
  'path': 'repro_sweep/n18/n18_eta0.90_d2_r21000_s4.pkl',
  'fidelity': 0.9987385258537915,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [1.4181456694956533, 1.816519445384406, -1.3932776881110873],
  'output': 's4'},
 {'target': 'n18_eta090_s0',
  'path': 'repro_sweep/n18/n18_eta0.90_d2_r21000_s0.pkl',
  'fidelity': 0.9987384698023508,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [-0.17051753378989587, 1.3672862011128144, 2.2008796271358415],
  'output': 's0'},
 {'target': 'n18_eta090_s2',
  'path': 'repro_sweep/n18/n18_eta0.90_d2_r21000_s2.pkl',
  'fidelity': 0.9987384445044466,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [4.918298361665541, 3.0314767962380054, -1.2219970508308002],
  'output': 's2'},
 {'target': 'n18_eta090_s3',
  'path': 'repro_sweep/n18/n18_eta0.90_d2_r21000_s3.pkl',
  'fidelity': 0.998738400960062,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [-2.784018244283829, 1.9640200495627793, 3.319116758938637],
  'output': 's3'},
 {'target': 'n18_eta090_s1',
  'path': 'repro_sweep/n18/n18_eta0.90_d2_r21000_s1.pkl',
  'fidelity': 0.9987375460564476,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [5.115898774729618, 0.9541890995449474, -0.26949698657380683],
  'output': 's1'}]
FAMILY_MINIMUM_OVERLAP = 0.9672129796582472
OVERLAPS_TO_BEST = {'n18_eta090_s4': 1.0, 'n18_eta090_s0': 0.9954199145414182, 'n18_eta090_s2': 0.9965613208078448, 'n18_eta090_s3': 0.9948613298855256, 'n18_eta090_s1': 0.9687389818464851}
ANALYTIC_COMPARISONS = [{'path': 'analytic_codes/n18_2D2_code.pkl',
  'classification': '2D2 reducible: rho_3 + rho_4',
  'closest_run_id': 'n18_eta090_s2',
  'overlap': 0.9845741714684988,
  'qualifies': True,
  'qualifying_run_ids': ['n18_eta090_s4', 'n18_eta090_s0', 'n18_eta090_s2', 'n18_eta090_s3'],
  'best_run_overlap': 0.978034238697603,
  'comparison_scripts': ['repro_analysis/n18/analytic_references/plot_n18_2D2_code.py']}]


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
