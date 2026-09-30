#!/usr/bin/env python3
"""eta095_f0: approximate structural family (0.95 connected links); per-run symmetry retained: 2D2: rho_1+rho_2; 2D2: rho_1+rho_4; 2T: rho_1+rho_2; 2T: rho_1+rho_3.

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
DEFAULT_ANALYTIC_FILE = 'analytic_codes/n12_2T_code.pkl'
BEST_RUN_ID = 'n12_eta095_s0'
CANONICAL_NUMERICAL = 'n12_eta095_s0'
OUTPUT_ROOT = HERE / "plots" / 'eta095_f0'
ANALYTIC_OUTPUT = ROOT / 'repro_analysis/n12/analytic_references/n12_2T_code'
VIEW = {'rotate': True, 'theta': 0.07315757, 'phi': 3.34027524}
# Human-readable titles printed on the figures; K is code rank/dimension.
NUMERICAL_TITLES = {'projector_sphere': 'Numerical code: projector\n$n=12$, $K=2$',
 'projector_plane': 'Numerical code: projector\n$n=12$, $K=2$',
 'logical0_sphere': 'Numerical code: logical state 0\n$n=12$, $K=2$',
 'logical0_plane': 'Numerical code: logical state 0\n$n=12$, $K=2$',
 'logical1_sphere': 'Numerical code: logical state 1\n$n=12$, $K=2$',
 'logical1_plane': 'Numerical code: logical state 1\n$n=12$, $K=2$'}
ANALYTIC_TITLES = {'projector_sphere': 'Constructed $2T$ code: projector\n$n=12$, $K=2$',
 'projector_plane': 'Constructed $2T$ code: projector\n$n=12$, $K=2$',
 'logical0_sphere': 'Constructed $2T$ code: logical state 0\n$n=12$, $K=2$',
 'logical0_plane': 'Constructed $2T$ code: logical state 0\n$n=12$, $K=2$',
 'logical1_sphere': 'Constructed $2T$ code: logical state 1\n$n=12$, $K=2$',
 'logical1_plane': 'Constructed $2T$ code: logical state 1\n$n=12$, $K=2$'}
# Analytic path and canonical orientation are deliberately separate.
# Standard mode: Each selected entry is sent to plot_file_simple directly, exactly once.
ENTRIES = [{'target': 'analytic',
  'path': 'analytic_codes/n12_2T_code.pkl',
  'fidelity': 0.9957620688022321,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [3.4795490998418366, 2.1138607910212497, 1.5583685786308503],
  'output': 'analytic'},
 {'target': 'n12_eta095_s0',
  'path': 'repro_sweep/n12/n12_eta0.95_d2_r21000_s0.pkl',
  'fidelity': 0.999618020444189,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [0.6869746232194367, 2.7604459491809896, -3.905582909061961],
  'output': 's0'},
 {'target': 'n12_eta095_s2',
  'path': 'repro_sweep/n12/n12_eta0.95_d2_r21000_s2.pkl',
  'fidelity': 0.9996180000615174,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [-0.22102978042883292, 1.0339392494253652, 2.082432896062777],
  'output': 's2'},
 {'target': 'n12_eta095_s4',
  'path': 'repro_sweep/n12/n12_eta0.95_d2_r21000_s4.pkl',
  'fidelity': 0.9996179393676252,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [0.23438405759808528, 1.2387702135371723, -1.0738821349175764],
  'output': 's4'},
 {'target': 'n12_eta095_s1',
  'path': 'repro_sweep/n12/n12_eta0.95_d2_r21000_s1.pkl',
  'fidelity': 0.9996155907639961,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [2.1745663209428114, 1.6644439295960178, -0.6510168212397919],
  'output': 's1'},
 {'target': 'n12_eta095_s3',
  'path': 'repro_sweep/n12/n12_eta0.95_d2_r21000_s3.pkl',
  'fidelity': 0.9996153975805109,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [-1.1243322444819426, 1.672300358624805, 5.0878141832591215],
  'output': 's3'}]
FAMILY_MINIMUM_OVERLAP = 0.94006639076971
OVERLAPS_TO_BEST = {'n12_eta095_s0': 1.0, 'n12_eta095_s2': 0.9999661463104191, 'n12_eta095_s4': 0.9996669339038972, 'n12_eta095_s1': 0.9497730201445179, 'n12_eta095_s3': 0.94006639076971}
ANALYTIC_COMPARISONS = [{'path': 'analytic_codes/n12_2T_code.pkl',
  'classification': '2T reducible: rho_1 + rho_2',
  'closest_run_id': 'n12_eta095_s0',
  'overlap': 0.9992397619881066,
  'qualifies': True,
  'qualifying_run_ids': ['n12_eta095_s0', 'n12_eta095_s2', 'n12_eta095_s4'],
  'best_run_overlap': 0.9992397619881066,
  'comparison_scripts': ['repro_analysis/n12/analytic_references/plot_n12_2T_code.py']}]


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
