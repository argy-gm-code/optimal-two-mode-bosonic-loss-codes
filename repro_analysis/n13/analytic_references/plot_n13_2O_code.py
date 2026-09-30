#!/usr/bin/env python3
"""One fixed projector viewing frame for this analytic construction.

The saved fidelity belongs to saved_eta. A viewing rotation does not reevaluate
the fidelity or create a different code at another noise strength.
"""
import argparse
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DEFAULT_ANALYTIC_FILE = 'analytic_codes/n13_2O_code.pkl'
OUTPUT_ROOT = HERE / 'n13_2O_code'
REFERENCE = {'fidelity': 0.9974664033578942, 'saved_eta': 0.9}
ANALYTIC_TITLES = {'projector_sphere': 'Constructed $2O$ code: projector\n$n=13$, $K=2$',
 'projector_plane': 'Constructed $2O$ code: projector\n$n=13$, $K=2$',
 'logical0_sphere': 'Constructed $2O$ code: logical state 0\n$n=13$, $K=2$',
 'logical0_plane': 'Constructed $2O$ code: logical state 0\n$n=13$, $K=2$',
 'logical1_sphere': 'Constructed $2O$ code: logical state 1\n$n=13$, $K=2$',
 'logical1_plane': 'Constructed $2O$ code: logical state 1\n$n=13$, $K=2$'}
ANALYTIC_VIEW = {'reflect': False,
 'alignment_euler': [1.587107774809809, 4.412497388161908e-09, -1.5871077733070582],
 'view': {'rotate': True, 'theta': 0.7853981633974477, 'phi': 0.9553166181245097}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--target", choices=("analytic", "all"), default="analytic")
    parser.add_argument("--analytic-file", type=Path)
    parser.add_argument("--projector-only", action="store_true")
    parser.add_argument("--interactive", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--eta", help="Compatibility option; this construction always uses its single fixed view")
    parser.add_argument("--canonical", action="store_true", help="Compatibility alias for the single fixed view")
    args = parser.parse_args()
    source = args.analytic_file if args.analytic_file else Path(DEFAULT_ANALYTIC_FILE)
    if args.list:
        print("Targets: analytic (default), all; both plot the construction once")
        print("Analytic reference:", source, "saved_eta=", REFERENCE["saved_eta"], "Fe=", REFERENCE["fidelity"])
        print("Fixed viewing orientation:", ANALYTIC_VIEW)
        print("Output:", OUTPUT_ROOT)
        print("--eta and --canonical are compatibility options; no additional plot sets are produced.")
        return
    source = source if source.is_absolute() else ROOT / source
    if not source.is_file():
        parser.error(f"pickle does not exist: {source}")
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(ROOT / "code"))
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "repro_analysis_matplotlib"))
    import graph_n_photons_new
    if args.interactive and args.no_show:
        import plotly.io as pio
        pio.renderers.default = ""
    view = ANALYTIC_VIEW["view"]
    graph_n_photons_new.plot_file_simple(
        pkl_path=str(source), out_folder=str(OUTPUT_ROOT), interactive=args.interactive,
        projector_only=args.projector_only, reflect=ANALYTIC_VIEW["reflect"],
        alignment_euler=ANALYTIC_VIEW["alignment_euler"], rotate=view["rotate"],
        theta=view["theta"], phi=view["phi"], titles=ANALYTIC_TITLES,
    )


if __name__ == "__main__":
    main()
