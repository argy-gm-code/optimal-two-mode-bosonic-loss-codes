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
DEFAULT_ANALYTIC_FILE = None
BEST_RUN_ID = 'n15_eta090_s4'
CANONICAL_NUMERICAL = 'n15_eta090_s4'
OUTPUT_ROOT = HERE / "plots" / 'eta090_f0'
VIEW = {'rotate': False, 'theta': 0.0, 'phi': 0.0}
# Human-readable titles printed on the figures; K is code rank/dimension.
NUMERICAL_TITLES = {'projector_sphere': 'Numerical code: projector\n$n=15$, $K=2$',
 'projector_plane': 'Numerical code: projector\n$n=15$, $K=2$',
 'logical0_sphere': 'Numerical code: logical state 0\n$n=15$, $K=2$',
 'logical0_plane': 'Numerical code: logical state 0\n$n=15$, $K=2$',
 'logical1_sphere': 'Numerical code: logical state 1\n$n=15$, $K=2$',
 'logical1_plane': 'Numerical code: logical state 1\n$n=15$, $K=2$'}
ANALYTIC_TITLES = None
# Analytic path and canonical orientation are deliberately separate.
# Standard mode: Each selected entry is sent to plot_file_simple directly, exactly once.
ENTRIES = [{'target': 'n15_eta090_s4',
  'path': 'repro_sweep/n15/n15_eta0.90_d2_r21000_s4.pkl',
  'fidelity': 0.9978174886411162,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [-1.1102230246251565e-16, 2.761647890826747, -1.5792881036595117],
  'output': 's4'},
 {'target': 'n15_eta090_s2',
  'path': 'repro_sweep/n15/n15_eta0.90_d2_r21000_s2.pkl',
  'fidelity': 0.997817485767305,
  'saved_eta': 0.9,
  'reflect': True,
  'alignment_euler': [3.932604034420938, 2.3410336080665757, -2.2235797225930525],
  'output': 's2'},
 {'target': 'n15_eta090_s0',
  'path': 'repro_sweep/n15/n15_eta0.90_d2_r21000_s0.pkl',
  'fidelity': 0.9978174782311295,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [2.666956147798652, 0.4813375403615283, -1.7630458813043224],
  'output': 's0'},
 {'target': 'n15_eta090_s1',
  'path': 'repro_sweep/n15/n15_eta0.90_d2_r21000_s1.pkl',
  'fidelity': 0.9978172690228511,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [-4.702279891271375, 2.1173805614707253, 0.5702193709120085],
  'output': 's1'},
 {'target': 'n15_eta090_s3',
  'path': 'repro_sweep/n15/n15_eta0.90_d2_r21000_s3.pkl',
  'fidelity': 0.9978168820088293,
  'saved_eta': 0.9,
  'reflect': False,
  'alignment_euler': [4.400871642428669, 2.431596213416029, -1.7301316507712206],
  'output': 's3'}]
FAMILY_MINIMUM_OVERLAP = 0.9935186430186762
OVERLAPS_TO_BEST = {'n15_eta090_s4': 1.0, 'n15_eta090_s2': 0.9999778601306587, 'n15_eta090_s0': 0.9999020229762644, 'n15_eta090_s1': 0.9978134736939537, 'n15_eta090_s3': 0.9935186430186762}
ANALYTIC_COMPARISONS = [{'path': 'analytic_codes/n15_2T_code.pkl',
  'classification': '2T pure irrep: rho_5',
  'closest_run_id': 'n15_eta090_s4',
  'overlap': 0.8684440301697087,
  'qualifies': False,
  'qualifying_run_ids': [],
  'best_run_overlap': 0.8684440301697087,
  'comparison_scripts': ['repro_analysis/n15/analytic_references/plot_n15_2T_code.py']}]


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
