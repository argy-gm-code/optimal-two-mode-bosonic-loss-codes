#!/usr/bin/env python3
"""eta095_f0: approximate structural family (0.95 connected links); per-run symmetry retained: 2T: rho_5; C4: rho_3+rho_4; C6: rho_2+rho_4.

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
DEFAULT_ANALYTIC_FILE = 'analytic_codes/n15_2T_code.pkl'
BEST_RUN_ID = 'n15_eta095_s1'
CANONICAL_NUMERICAL = 'n15_eta095_s1'
OUTPUT_ROOT = HERE / "plots" / 'eta095_f0'
ANALYTIC_OUTPUT = ROOT / 'repro_analysis/n15/analytic_references/n15_2T_code'
VIEW = {'rotate': False, 'theta': 0.0, 'phi': 0.0}
# Human-readable titles printed on the figures; K is code rank/dimension.
NUMERICAL_TITLES = {'projector_sphere': 'Numerical code: projector\n$n=15$, $K=2$',
 'projector_plane': 'Numerical code: projector\n$n=15$, $K=2$',
 'logical0_sphere': 'Numerical code: logical state 0\n$n=15$, $K=2$',
 'logical0_plane': 'Numerical code: logical state 0\n$n=15$, $K=2$',
 'logical1_sphere': 'Numerical code: logical state 1\n$n=15$, $K=2$',
 'logical1_plane': 'Numerical code: logical state 1\n$n=15$, $K=2$'}
ANALYTIC_TITLES = {'projector_sphere': 'Constructed $2T$ code: projector\n$n=15$, $K=2$',
 'projector_plane': 'Constructed $2T$ code: projector\n$n=15$, $K=2$',
 'logical0_sphere': 'Constructed $2T$ code: logical state 0\n$n=15$, $K=2$',
 'logical0_plane': 'Constructed $2T$ code: logical state 0\n$n=15$, $K=2$',
 'logical1_sphere': 'Constructed $2T$ code: logical state 1\n$n=15$, $K=2$',
 'logical1_plane': 'Constructed $2T$ code: logical state 1\n$n=15$, $K=2$'}
# Analytic path and canonical orientation are deliberately separate.
# Standard mode: Each selected entry is sent to plot_file_simple directly, exactly once.
ENTRIES = [{'target': 'analytic',
  'path': 'analytic_codes/n15_2T_code.pkl',
  'fidelity': 0.9998690006307915,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [-0.9538177419019098, 4.29825129755929e-09, 0.9538177279657734],
  'output': 'analytic'},
 {'target': 'n15_eta095_s1',
  'path': 'repro_sweep/n15/n15_eta0.95_d2_r21000_s1.pkl',
  'fidelity': 0.9998690035709338,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [-2.502572080771218, 0.5943946312727232, 3.6613399307392624],
  'output': 's1'},
 {'target': 'n15_eta095_s4',
  'path': 'repro_sweep/n15/n15_eta0.95_d2_r21000_s4.pkl',
  'fidelity': 0.9998686782730225,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [0.3438433037568324, 0.4643088470762981, -1.0703180849696532],
  'output': 's4'},
 {'target': 'n15_eta095_s3',
  'path': 'repro_sweep/n15/n15_eta0.95_d2_r21000_s3.pkl',
  'fidelity': 0.9998670170480581,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [3.396146600688468, 1.2372873945614853, -2.187658861356472],
  'output': 's3'},
 {'target': 'n15_eta095_s7',
  'path': 'repro_sweep/n15/n15_eta0.95_d2_r21000_s7.pkl',
  'fidelity': 0.9998662471247488,
  'saved_eta': 0.95,
  'reflect': False,
  'alignment_euler': [3.402009797306456, 1.1821687467883266, -2.370743613039984],
  'output': 's7'},
 {'target': 'n15_eta095_s5',
  'path': 'repro_sweep/n15/n15_eta0.95_d2_r21000_s5.pkl',
  'fidelity': 0.9998644986766251,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [-0.575221598765989, 0.17133721266897037, 3.617867471170454],
  'output': 's5'},
 {'target': 'n15_eta095_s9',
  'path': 'repro_sweep/n15/n15_eta0.95_d2_r21000_s9.pkl',
  'fidelity': 0.9998637135809407,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [-1.1380421548715731, 0.2289771006265721, 3.0549871941877873],
  'output': 's9'},
 {'target': 'n15_eta095_s2',
  'path': 'repro_sweep/n15/n15_eta0.95_d2_r21000_s2.pkl',
  'fidelity': 0.9998635437633986,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [-1.1348031510241392, 1.7581594156151155, 0.06351193271333022],
  'output': 's2'},
 {'target': 'n15_eta095_s0',
  'path': 'repro_sweep/n15/n15_eta0.95_d2_r21000_s0.pkl',
  'fidelity': 0.9998587012229823,
  'saved_eta': 0.95,
  'reflect': True,
  'alignment_euler': [-4.005217955569584, 2.3499144022564353, 2.233248384705217],
  'output': 's0'}]
FAMILY_MINIMUM_OVERLAP = 0.8826678339823102
OVERLAPS_TO_BEST = {'n15_eta095_s1': 1.0, 'n15_eta095_s4': 0.9962152700567232, 'n15_eta095_s3': 0.9579951076457012, 'n15_eta095_s7': 0.9656789990826204, 'n15_eta095_s5': 0.9527601444761586, 'n15_eta095_s9': 0.9161157103579124, 'n15_eta095_s2': 0.9164388883462757, 'n15_eta095_s0': 0.8967719362943857}
ANALYTIC_COMPARISONS = [{'path': 'analytic_codes/n15_2T_code.pkl',
  'classification': '2T pure irrep: rho_5',
  'closest_run_id': 'n15_eta095_s1',
  'overlap': 0.9748832145217277,
  'qualifies': True,
  'qualifying_run_ids': ['n15_eta095_s1', 'n15_eta095_s4'],
  'best_run_overlap': 0.9748832145217277,
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
