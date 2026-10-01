# Optimal two-mode bosonic loss codes from finite group symmetry

Code and data accompanying the paper by Argyris Giannisis Manes,
Mahadevan Subramanian, and Liang Jiang (2026).

## Install

Verified with Python 3.9.6 and 3.12.8. Run commands from this directory:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r code/requirements.txt
source .venv/bin/activate
```

After activation, `python` uses this environment.

GAP is optional; regenerating the bundled group JSON uses RepnDecomp and repsn.

## Files

- `code/`: optimization, representation theory, plotting, Farkas certificates, and group data/GAP sources.
- `repro_sweep/`: 289 numerical runs at 76 parameter points.
- `analytic_codes/`, `additional_codes/`: 37 constructed codes.
- `literature_codes/`: the two n=25 comparison inputs and their reevaluation record.
- `repro_analysis/`: classifications, overlaps, performance, plotting recipes, and 319 plane plots.
- `Distance Certificates/`: 52 exact certificates.
- `n28_symbolic_2O_code/`: exact symbolic construction and verification.

Pickles contain `(n, eta, logical_dimension, rounds, seed, decoder, encoder,
fidelity, coherent_information)`. Only load trusted pickle files.
`runs.json` records saved results; `performance.json` separately reevaluates fixed
encoders with `eps=1e-9`. The extra n=18 refinement is under `repro_analysis/n18/`.
`constructions.json` identifies the larger family and additional distance constructions.
Numerical KL checks use the stated tolerance; only the symbolic and certificate
checks use exact arithmetic. Family labels use 0.95-overlap connected components;
completed runs and family membership do not certify global optimality.

## Verify saved results

Short IDs such as `n9_eta090_s3` map to pickle paths in
`repro_analysis/runs.json`. The overlap checker also accepts recorded file paths;
it replays a stored comparison rather than searching for a new alignment.

```sh
# Saved results, 898 overlaps/124 families, constructed distances, exact certificates:
python -B repro_analysis/verify_results.py all
# Replay all 76 stored symmetry fits:
python -B repro_analysis/check_symmetry.py --all
# Exact n=28 verification, without regenerating its JSON:
python -B -c "from n28_symbolic_2O_code.n28_symbolic_2O_code import verify_code; verify_code()"
# Example stored overlap:
python -B repro_analysis/check_overlap.py n9_eta090_s3 n9_C6_code
```

## Regenerate plots

```sh
# Plot one family; --list shows inputs. Add --interactive --no-show for HTML:
python -B repro_analysis/n9/plot_eta090_f0.py --target all --projector-only
```

Plot recipes regenerate plane and sphere images in their output folders,
replacing existing plots. Omit `--projector-only` to include logical states.
Use `--target numericals` for all numerical members and `--target analytic` for
the constructed reference. Analytic-only recipes are under
`repro_analysis/n*/analytic_references/`; n=37 has its own recipe.
For any other pickle, call `plot_file_simple(path, interactive=False)` after
adding `code/` to Python's import path. Use the recipes rather than running
`graph_n_photons_new.py` or `group_rep_codes_classes.py` directly.

## Run new calculations

Decoder reevaluations keep the encoder fixed and print JSON. The certificate
and optimization examples save new files under `/tmp`; choose another output
path to keep them permanently.

```sh
# Strict n=11 fixed-encoder reevaluation:
python -B repro_analysis/evaluate_performance.py n11_eta090_s0 --eta 0.90 --eps 1e-9
# Reproduce the n=25 comparison:
python -B repro_analysis/evaluate_performance.py literature_codes/literature_n25_og_2o_eta0.900000.pkl literature_codes/constructed_n25_2o_eta0.900000.pkl --eta 0.90 --eps 1e-9
# Generate a new certificate in a separate output folder:
python -B code/Farkas_certificates.py --n-min 4 --n-max 4 --K 2 --output-directory /tmp/new_certificates
# Start a small new optimization, writing a new result:
python -B code/faster_gen_dual_rail.py --n 4 --eta 0.90 --d 2 --rounds 2 --seed 0 --out /tmp/new_code.pkl --workers 1
```

Hashes identify exact input/source files.
