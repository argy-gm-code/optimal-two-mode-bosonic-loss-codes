"""Generate exact MacWilliams--Farkas distance certificates.

For fixed photon number ``n``, code dimension ``K``, and proposed loss
distance ``d``, this module constructs the matrices in the Farkas system

    U x <= 0,        c.T x < 0.

All entries are SymPy rationals.  A returned certificate is normalized so
that ``c.T x == -1`` and is verified exactly before it can be saved.
"""

from __future__ import annotations

import argparse
import json
from functools import lru_cache
from pathlib import Path
from typing import Iterable

from sympy import (
    ImmutableMatrix,
    Integer,
    Matrix,
    Rational,
    eye,
    ones,
    simplify,
    symbols,
    zeros,
)
from sympy.matrices.matrixbase import MatrixBase
from sympy.physics.wigner import wigner_6j
from sympy.solvers.simplex import InfeasibleLPError, UnboundedLPError, lpmin


DEFAULT_CERTIFICATE_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "Distance Certificates"
)


def _require_integer(name: str, value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer; got {value!r}.")


def _validate_parameters(n: int, K: int, d: int | None = None) -> None:
    _require_integer("n", n)
    _require_integer("K", K)
    if n < 0:
        raise ValueError(f"n must be nonnegative; got {n}.")
    if K <= 0:
        raise ValueError(f"K must be positive; got {K}.")
    if K > n + 1:
        raise ValueError(
            f"A rank-{K} code cannot fit in H_{n}, which has dimension {n + 1}."
        )

    if d is not None:
        _require_integer("d", d)
        if not 1 <= d <= n + 1:
            raise ValueError(f"d must satisfy 1 <= d <= n + 1; got d={d}, n={n}.")


@lru_cache(maxsize=None)
def generate_macwilliams_matrix(n: int) -> ImmutableMatrix:
    """Return the exact SU(2) MacWilliams matrix ``M_n``.

    The row and column indices are ``l, t = 0, ..., n`` and

        (M_n)[l,t] = (-1)^(n+l+t) (2l+1)
                       {n/2 n/2 l; n/2 n/2 t}.

    For these arguments the Racah formula is rational, so every entry is
    required to simplify to a SymPy rational number.
    """

    _require_integer("n", n)
    if n < 0:
        raise ValueError(f"n must be nonnegative; got {n}.")

    dimension = n + 1
    spin = Rational(n, 2)

    def entry(l: int, t: int):
        value = simplify(
            Integer(-1) ** (n + l + t)
            * Integer(2 * l + 1)
            * wigner_6j(spin, spin, l, spin, spin, t)
        )
        if value.is_Rational is not True:
            raise ArithmeticError(
                f"MacWilliams entry ({l}, {t}) did not simplify to a rational: "
                f"{value!r}."
            )
        return value

    return ImmutableMatrix(dimension, dimension, entry)


def generate_farkas_system(
    n: int, K: int, d: int
) -> tuple[ImmutableMatrix, ImmutableMatrix]:
    """Construct the exact matrices ``U`` and ``c`` for ``(n, K, d)``.

    The certificate vector has block order

        (slack relation, sum A, sum B, MacWilliams relation, detection),

    with block sizes ``(n+1, 1, 1, n+1, n+1)``.  Thus ``U`` has shape
    ``(3(n+1), 3n+5)`` and ``c`` has shape ``(3n+5, 1)``.

    The diagonal matrix ``I_d`` selects tensor ranks ``0, ..., d-1``.
    """

    _validate_parameters(n, K, d)

    dimension = n + 1
    identity = eye(dimension)
    one_vector = ones(dimension, 1)
    zero_vector = zeros(dimension, 1)
    zero_matrix = zeros(dimension, dimension)
    macwilliams = Matrix(generate_macwilliams_matrix(n))
    detection = Matrix.diag(
        *([Integer(1)] * d + [Integer(0)] * (dimension - d))
    )

    U = Matrix.vstack(
        Matrix.hstack(
            -identity,
            one_vector,
            zero_vector,
            macwilliams.T,
            -detection,
        ),
        Matrix.hstack(
            K * identity,
            zero_vector,
            one_vector,
            -identity,
            K * detection,
        ),
        Matrix.hstack(
            -identity,
            zero_vector,
            zero_vector,
            zero_matrix,
            zero_matrix,
        ),
    )

    c = -Matrix.vstack(
        zero_vector,
        Matrix([K]),
        Matrix([K**2]),
        zero_vector,
        zero_vector,
    )

    expected_columns = 3 * n + 5
    if U.shape != (3 * dimension, expected_columns):
        raise AssertionError(f"Unexpected U shape {U.shape}.")
    if c.shape != (expected_columns, 1):
        raise AssertionError(f"Unexpected c shape {c.shape}.")
    if not all(value.is_Rational is True for value in U):
        raise ArithmeticError("U contains a non-rational entry.")

    return ImmutableMatrix(U), ImmutableMatrix(c)


def verify_farkas_certificate(
    U: MatrixBase,
    c: MatrixBase,
    certificate: MatrixBase | Iterable,
    *,
    require_normalized: bool = True,
) -> bool:
    """Check a Farkas certificate using exact arithmetic."""

    U = Matrix(U)
    c = Matrix(c)
    certificate = Matrix(certificate)

    if certificate.cols != 1:
        certificate = certificate.reshape(len(certificate), 1)
    if U.cols != certificate.rows or c.shape != (U.cols, 1):
        return False

    inequalities = U * certificate
    if any(value.is_Rational is not True for value in inequalities):
        return False
    if any(value > 0 for value in inequalities):
        return False

    objective = (c.T * certificate)[0]
    if objective.is_Rational is not True:
        return False
    if require_normalized:
        return objective == -1
    return objective < 0


def solve_farkas_system(
    U: MatrixBase, c: MatrixBase
) -> ImmutableMatrix | None:
    """Find an exact certificate for ``U x <= 0, c.T x < 0``.

    The certificate conditions are homogeneous.  Therefore a certificate
    exists if and only if the exact linear program

        minimize    c.T x
        subject to  U x <= 0,
                    -1 <= x_i <= 1

    has a negative optimum.  The box only fixes an arbitrary scale and makes
    the optimization bounded.  A negative optimizer is rescaled exactly to
    satisfy ``c.T x == -1``.

    ``None`` means that no Farkas certificate exists, equivalently that the
    corresponding nonnegative MacWilliams system is feasible.
    """

    U = Matrix(U)
    c = Matrix(c)
    if c.shape != (U.cols, 1):
        raise ValueError(
            f"Expected c to have shape ({U.cols}, 1); got {c.shape}."
        )
    if not all(value.is_Rational is True for value in U) or not all(
        value.is_Rational is True for value in c
    ):
        raise ValueError("Exact certificate solving requires rational U and c.")

    variables = symbols(f"x0:{U.cols}", real=True)
    vector = Matrix(variables)
    constraints = [value <= 0 for value in U * vector]
    constraints.extend(variable <= 1 for variable in variables)
    constraints.extend(variable >= -1 for variable in variables)
    objective = (c.T * vector)[0]

    try:
        optimum, solution = lpmin(objective, constraints)
    except (InfeasibleLPError, UnboundedLPError) as exc:
        raise RuntimeError(
            "The bounded certificate LP should always contain x=0 and cannot "
            "be unbounded."
        ) from exc

    if optimum == 0:
        return None
    if optimum > 0:
        raise ArithmeticError(
            f"Certificate LP returned positive optimum {optimum}, although x=0 "
            "is feasible."
        )

    optimizer = Matrix([solution.get(variable, Integer(0)) for variable in variables])
    optimizer_value = simplify((c.T * optimizer)[0])
    if optimizer_value >= 0 or any(value > 0 for value in U * optimizer):
        raise ArithmeticError("SymPy returned an invalid certificate optimizer.")

    certificate = simplify(optimizer / (-optimizer_value))
    if not verify_farkas_certificate(U, c, certificate):
        raise ArithmeticError("Exact certificate verification failed.")
    return ImmutableMatrix(certificate)


def find_farkas_certificate(n: int, K: int, d: int) -> ImmutableMatrix | None:
    """Generate and solve the exact Farkas system for ``(n, K, d)``."""

    U, c = generate_farkas_system(n, K, d)
    return solve_farkas_system(U, c)


def find_first_infeasible_distance(
    n: int, K: int = 2, start_d: int = 1, *, verbose: bool = False
) -> tuple[int, ImmutableMatrix] | None:
    """Return the first distance whose MacWilliams LP is infeasible.

    Starting at ``start_d``, this increases ``d`` until an exact Farkas
    certificate is found.  A certificate at ``d`` proves that no rank-``K``
    code on ``H_n`` can correct ``d-1`` losses.
    """

    _validate_parameters(n, K)
    _require_integer("start_d", start_d)
    if start_d < 1:
        raise ValueError(f"start_d must be positive; got {start_d}.")

    for d in range(start_d, n + 2):
        if verbose:
            print(f"n={n}, K={K}: testing d={d}")
        certificate = find_farkas_certificate(n, K, d)
        if certificate is not None:
            if verbose:
                print(f"n={n}, K={K}: first infeasible distance is d={d}")
            return d, certificate
    return None


def _rational_string(value) -> str:
    value = simplify(value)
    if value.is_Rational is not True:
        raise ValueError(f"Expected a rational value; got {value!r}.")
    return str(value)


def certificate_filename(n: int, K: int, d: int) -> str:
    """Return the standard JSON filename for a certificate."""

    _validate_parameters(n, K, d)
    return f"farkas_certificate_n{n}_K{K}_d{d}.json"


def save_farkas_certificate(
    certificate: MatrixBase | Iterable,
    n: int,
    K: int,
    d: int,
    output_directory: str | Path | None = None,
    *,
    overwrite: bool = False,
) -> Path:
    """Verify and save an exact certificate as JSON."""

    U, c = generate_farkas_system(n, K, d)
    certificate = Matrix(certificate)
    if certificate.cols != 1:
        certificate = certificate.reshape(len(certificate), 1)
    if not verify_farkas_certificate(U, c, certificate):
        raise ValueError("Refusing to save a certificate that fails exact verification.")

    directory = (
        DEFAULT_CERTIFICATE_DIRECTORY
        if output_directory is None
        else Path(output_directory)
    )
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / certificate_filename(n, K, d)
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Certificate already exists at {path}. Pass overwrite=True to replace it."
        )

    dimension = n + 1
    residual = U * certificate
    payload = {
        "schema_version": 1,
        "parameters": {"n": n, "K": K, "d": d},
        "meaning": (
            "Exact Farkas certificate proving infeasibility of the intrinsic "
            "MacWilliams LP."
        ),
        "normalization": "c^T x = -1",
        "vector_order": [
            {"name": "slack_relation", "start": 0, "length": dimension},
            {"name": "sum_A", "start": dimension, "length": 1},
            {"name": "sum_B", "start": dimension + 1, "length": 1},
            {
                "name": "macwilliams_relation",
                "start": dimension + 2,
                "length": dimension,
            },
            {
                "name": "detection_relation",
                "start": 2 * dimension + 2,
                "length": dimension,
            },
        ],
        "certificate": [_rational_string(value) for value in certificate],
        "exact_checks": {
            "c_transpose_x": _rational_string((c.T * certificate)[0]),
            "maximum_Ux": _rational_string(max(residual)),
            "all_Ux_nonpositive": True,
        },
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


def scan_photon_numbers(
    n_min: int,
    n_max: int,
    K: int = 2,
    start_d: int = 1,
    output_directory: str | Path | None = None,
    *,
    overwrite: bool = False,
    verbose: bool = False,
) -> list[Path]:
    """Find and save the first infeasibility certificate for each ``n``.

    After finding the first infeasible distance at one photon number, begin the
    next photon-number scan at that distance. This uses the monotonicity that
    the first infeasible distance is nondecreasing with ``n``.
    """

    _require_integer("n_min", n_min)
    _require_integer("n_max", n_max)
    _require_integer("K", K)
    if n_min < 0 or n_max < n_min:
        raise ValueError(f"Require 0 <= n_min <= n_max; got {n_min}, {n_max}.")
    if K <= 0:
        raise ValueError(f"K must be positive; got {K}.")
    if n_min + 1 < K:
        raise ValueError(
            f"n_min={n_min} is too small for a rank-{K} code; require n_min >= {K - 1}."
        )

    saved_paths: list[Path] = []
    current_start_d = start_d
    for n in range(n_min, n_max + 1):
        result = find_first_infeasible_distance(
            n, K=K, start_d=current_start_d, verbose=verbose
        )
        if result is None:
            if verbose:
                print(f"n={n}, K={K}: no certificate found for d <= {n + 1}")
            continue

        d, certificate = result
        current_start_d = d
        try:
            path = save_farkas_certificate(
                certificate,
                n,
                K,
                d,
                output_directory=output_directory,
                overwrite=overwrite,
            )
        except FileExistsError:
            directory = (
                DEFAULT_CERTIFICATE_DIRECTORY
                if output_directory is None
                else Path(output_directory)
            )
            path = directory / certificate_filename(n, K, d)
            if verbose:
                print(f"Keeping existing certificate {path}")
        else:
            if verbose:
                print(f"Saved {path}")
        saved_paths.append(path)

    return saved_paths


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate exact Farkas certificates for intrinsic SU(2) "
            "MacWilliams distance bounds."
        )
    )
    parser.add_argument("--n-min", type=int, default=1)
    parser.add_argument("--n-max", type=int, default=10)
    parser.add_argument("--K", type=int, default=2, help="Code dimension.")
    parser.add_argument("--start-d", type=int, default=1)
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=DEFAULT_CERTIFICATE_DIRECTORY,
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing certificate JSON files.",
    )
    return parser


def main() -> None:
    args = _build_argument_parser().parse_args()
    paths = scan_photon_numbers(
        n_min=args.n_min,
        n_max=args.n_max,
        K=args.K,
        start_d=args.start_d,
        output_directory=args.output_directory,
        overwrite=args.overwrite,
        verbose=True,
    )
    print(f"Completed scan; {len(paths)} certificate path(s) available.")


if __name__ == "__main__":
    main()
