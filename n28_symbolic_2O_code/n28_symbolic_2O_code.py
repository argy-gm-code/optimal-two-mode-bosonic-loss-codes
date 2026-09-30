"""Exact n=28, 2O:rho_3 codewords from multiplicity-space KL equations.

All vectors use the Fock basis e_m = |m, 28-m>.  Gross's generators act on
the target rho_3 basis as S = diag(1,-1) and
H = [[-1/2,-sqrt(3)/2],[-sqrt(3)/2,1/2]].
"""

import json
from functools import lru_cache
from math import comb, factorial
from pathlib import Path

import sympy as sp


N = 28
LABELS = (2, 6, 10)
CODE_FILE = Path(__file__).with_suffix(".json")
#The physical R_{28}(S) matrix of the S gate
S28 = sp.diag(*[-(-sp.I) ** m for m in range(N + 1)])
S2 = sp.diag((1 - sp.I) / sp.sqrt(2), (1 + sp.I) / sp.sqrt(2))
H2 = -sp.I * sp.Matrix([[1, 1], [1, -1]]) / sp.sqrt(2)


def e(m):
    state = sp.zeros(N + 1, 1)
    state[m] = 1
    return state


def f(m):
    """A Q8-fixed vector; m is even and 0 <= m <= 14."""
    return e(14) if m == 14 else (e(m) + e(N - m)) / sp.sqrt(2)


@lru_cache(None)
def h_column(m):
    # Computes the action of the physical Hadamard matrix
    """R_28(H)e_m from the binomial action on creation operators."""
    amplitudes = []
    for r in range(N + 1):
        coefficient = sum(
            (-1) ** (N - m - r + t) * comb(m, t) * comb(N - m, r - t)
            for t in range(max(0, r - N + m), min(m, r) + 1)
        )
        scale = sp.sqrt(sp.Rational(factorial(r) * factorial(N - r),
                                    factorial(m) * factorial(N - m)))
        amplitudes.append(sp.Rational(coefficient, 2 ** 14) * scale)
    return sp.Matrix(amplitudes)


def h(state):
    """Physical action of Gross's H on a 29-entry Fock vector."""
    return sum((state[m] * h_column(m) for m in range(N + 1)
                if state[m] != 0), sp.zeros(N + 1, 1)).applyfunc(sp.simplify)


def elementary_tensors(rank, n=N):
    """Return T_mu^(rank) on H_n, with T_rank^(rank) = J_+^rank."""
    if not 0 <= rank <= n:
        raise ValueError("rank must lie between 0 and n")
    j_plus = sp.zeros(n + 1)
    for m in range(n):
        j_plus[m + 1, m] = sp.sqrt((m + 1) * (n - m))
    j_minus = j_plus.T
    tensors = {rank: j_plus ** rank}
    for mu in range(rank, -rank, -1):
        commutator = j_minus * tensors[mu] - tensors[mu] * j_minus
        tensors[mu - 1] = (commutator / sp.sqrt((rank + mu) * (rank - mu + 1))).applyfunc(sp.simplify)
    return tensors


def two_o_orbit(rotation):
    """Map each projective 2O orbit point to one exact rotation g*rotation."""
    a, b = rotation[0, 0], rotation[1, 0]
    if sp.simplify(sp.conjugate(a) * a + sp.conjugate(b) * b - 1) != 0:
        raise ValueError("the first column of the rotation must have unit norm")
    orbit = {}
    pending = [sp.Matrix([a, b])]
    for column in pending:
        first = sp.simplify(column[0])
        key = None if first == 0 else sp.simplify(sp.cancel(column[1] / first))
        if key in orbit:
            continue
        second = sp.simplify(column[1])
        orbit[key] = sp.Matrix([[first, -sp.conjugate(second)],
                                [second, sp.conjugate(first)]])
        if len(orbit) > 24:
            raise ValueError("could not identify the projective 2O orbit exactly")
        pending.extend((S2 * column, H2 * column))
    return orbit


def seeds_at_rank(rank, n=N):
    """Return the levelwise-minimal (rotation, seed) pairs at tensor rank rank."""
    tensors = elementary_tensors(rank, n)
    count = (2 * rank + 24) // 24  # ceil((2 rank + 1)/|2O/{+/-I}|)
    seeds, used_points = [], set()
    candidate = 2
    while len(seeds) < count:
        a = 1 / sp.sqrt(1 + candidate ** 2)
        b = candidate * a
        rotation = sp.Matrix([[a, -sp.conjugate(b)], [b, sp.conjugate(a)]])
        orbit = two_o_orbit(rotation)
        candidate += 1
        if len(orbit) != 24 or used_points.intersection(orbit):
            continue
        used_points.update(orbit)
        seed = sum((sp.sqrt(sp.binomial(2 * rank, rank + mu))
                    * a ** (rank + mu) * b ** (rank - mu) * tensors[mu]
                    for mu in range(-rank, rank + 1)), sp.zeros(n + 1))
        seeds.append((rotation, seed))
    return seeds


def f_matrices():
    """Hermitian End(rho_3) basis in the logical basis of build_basis()."""
    x = sp.Matrix([[0, 1], [1, 0]])
    y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    z = sp.diag(1, -1)
    return {"trivial": (sp.eye(2),),
            "rho_2": (y / sp.sqrt(2),),
            "rho_3": (z / sp.sqrt(2), -x / sp.sqrt(2))}


def build_basis():
    """Return the trivial vector, coefficients c, and paired 29x3 bases V,W."""
    f14 = f(14)
    hf14 = h(f14)
    trivial = f14 + hf14 + S28 * hf14  # S3 average, up to a factor of 1/3; this is the vector of the trivial representation
    trivial = (trivial / trivial[14]).applyfunc(sp.simplify)
    c = sp.Matrix([-sp.simplify((f(m).T * trivial)[0]) for m in LABELS])
    V = sp.Matrix.hstack(*(f(m) + c[i] * f14 for i, m in enumerate(LABELS)))
    W = sp.Matrix.hstack(*(
        (-(V[:, i] + 2 * h(V[:, i])) / sp.sqrt(3)).applyfunc(sp.simplify)
        for i in range(3)
    ))
    return trivial, c, V, W


def isometry(c, V, W):
    """Return S_C with columns ordered by (multiplicity index, logical index)."""
    norm_c_squared = (c.T * c)[0]
    gram_inverse_sqrt = (sp.eye(3) +
                         ((1 / sp.sqrt(1 + norm_c_squared) - 1) / norm_c_squared)
                         * c * c.T)
    U, Z = V * gram_inverse_sqrt, W * gram_inverse_sqrt
    return sp.Matrix.hstack(*(column for i in range(3)
                              for column in (U[:, i], Z[:, i])))


def b_matrices(operator, V, W):
    """Return the B matrices for paired columns V,W and the F basis.

    Pass the orthonormal columns from isometry() for the canonical B matrices.
    Passing the raw columns from build_basis() gives their Gram-congruent form.
    """
    blocks = ((V.H * operator * V, V.H * operator * W),
              (W.H * operator * V, W.H * operator * W))
    return {
        name: tuple((sum((sp.conjugate(F[i, j]) * blocks[i][j]
                          for i in range(2) for j in range(2)), sp.zeros(3))
                     / sp.trace(F.H * F)).applyfunc(sp.simplify)
                    for F in family)
        for name, family in f_matrices().items()
    }


def solve_coefficients(c, V, W):
    """Solve the four rank-2,...,5 equations, with a_2 fixed to 1."""
    u, v, w, z = sp.symbols("u v w z", real=True)
    x = sp.Matrix([1, u + sp.I * v, w + sp.I * z])
    rescale = sp.diag(*(1 / entry for entry in c))
    equations = []
    choices = ((2, "rho_3", sp.sqrt(2)),
               (3, "rho_2", sp.I * sp.sqrt(6)),
               (4, "rho_3", sp.sqrt(2)),
               (5, "rho_3", sp.I * sp.sqrt(2)))
    for rank, irrep, factor in choices:
        operator = seeds_at_rank(rank)[0][1]
        B = b_matrices(operator, V, W)
        chosen = B[irrep][0]
        if irrep == "rho_3":
            other = B[irrep][1]
            pivot = next((i, j) for i in range(3) for j in range(3)
                         if chosen[i, j] != 0)
            ratio = sp.simplify(other[pivot] / chosen[pivot])
            assert all(sp.simplify(entry) == 0 for entry in other - ratio * chosen)
            assert all(entry == 0 for entry in B["rho_2"][0])
        else:
            assert all(entry == 0 for matrix in B["rho_3"] for entry in matrix)
        form = (rescale * chosen * rescale / factor).applyfunc(sp.simplify)
        equation = sp.Poly(sp.expand((x.H * form * x)[0]), u, v, w, z)
        _, integer_equation = equation.clear_denoms(convert=True)
        _, integer_equation = integer_equation.primitive()
        equations.append(integer_equation.as_expr())

    basis = sp.groebner(equations, u, v, w, z, order="lex", method="f5b")
    quartic = basis.polys[-1].as_expr()
    positive_roots = [root for root in sp.solve(quartic, z)
                      if root.is_positive is True]
    if len(positive_roots) != 1:
        raise ValueError("expected one positive real root of the eliminant")
    z_value = positive_roots[0]
    solution = sp.solve([p.as_expr().subs(z, z_value) for p in basis.polys[:-1]],
                        (u, v, w), dict=True)[0]
    raw = rescale * sp.Matrix([1, solution[u] + sp.I * solution[v],
                               solution[w] + sp.I * z_value])
    return raw / raw[0], quartic


def candidate_codewords(coefficients):
    """Return normalized |0_L>, |1_L>, and their common raw norm squared.

    The coefficients are (a_2,a_6,a_10).  A generic choice has rho_3
    symmetry but need not satisfy the loss-correction conditions.
    """
    if len(coefficients) != 3:
        raise ValueError("expected (a_2, a_6, a_10)")
    a = sp.Matrix([sp.sympify(value) for value in coefficients])
    if all(value == 0 for value in a):
        raise ValueError("the coefficient vector must be nonzero")
    _, c, V, W = build_basis()
    center = (c.T * a)[0]
    norm_squared = (a.H * a)[0] + sp.conjugate(center) * center
    return V * a / sp.sqrt(norm_squared), W * a / sp.sqrt(norm_squared), norm_squared


def verify(trivial, c, V, W):
    """Check the exact coefficients, generator action, and pair Gram matrix."""
    zero = lambda matrix: all(sp.simplify(x) == 0 for x in matrix)
    expected = [sp.Rational(437, 238), sp.Rational(741, 595),
                sp.Rational(1287, 1190)]
    assert all(sp.simplify(c[i] ** 2 - expected[i]) == 0 for i in range(3))
    assert zero(S28 * trivial - trivial) and zero(h(trivial) - trivial)
    assert zero(S28 * V - V) and zero(S28 * W + W)
    assert zero(V.H * trivial)
    gram = sp.eye(3) + c * c.T
    assert zero(V.H * V - gram) and zero(W.H * W - gram)
    assert zero(V.H * W)
    for i in range(3):
        assert zero(h(W[:, i]) + sp.sqrt(3) * V[:, i] / 2 - W[:, i] / 2)


def derive_code(path=CODE_FILE):
    """Derive and save exact coefficients in the paired (v_m, w_m) basis."""
    trivial, c, V, W = build_basis()
    verify(trivial, c, V, W)
    S_C = isometry(c, V, W)
    assert all(sp.simplify(entry) == 0 for entry in S_C.H * S_C - sp.eye(6))
    coefficients, quartic = solve_coefficients(c, V, W)
    for operator in elementary_tensors(1).values():
        B = b_matrices(operator, V, W)
        assert all(entry == 0 for name in ("rho_2", "rho_3")
                   for matrix in B[name] for entry in matrix)
    seed_rotations = []
    for rank in range(2, 6):
        rotation, seed = seeds_at_rank(rank)[0]
        seed_rotations.append((rank, rotation))
        B = b_matrices(seed, V, W)
        assert all(sp.simplify((coefficients.H * matrix * coefficients)[0]) == 0
                   for name in ("rho_2", "rho_3") for matrix in B[name])

    path = Path(path)
    path.write_text(json.dumps({
        "n": N,
        "group": "2O",
        "irrep": "rho_3",
        "basis": "|0_L> = V a / norm, |1_L> = W a / norm; V,W from build_basis()",
        "coefficients": {f"a_{m}": sp.sstr(coefficients[i])
                         for i, m in enumerate(LABELS)},
    }, indent=2) + "\n", encoding="utf-8")

    x6 = sp.simplify(c[1] * coefficients[1] / c[0])
    x10 = sp.simplify(c[2] * coefficients[2] / c[0])
    print("n=28, 2O:rho_3 code in the paired basis (v_m, w_m), m=2,6,10.")
    print("The basis coefficients a_m give |0_L> proportional to sum a_m v_m")
    print("and |1_L> proportional to sum a_m w_m, with a_2 = 1.")
    print("Exact basis checks passed; c_m^2 =", [sp.simplify(x**2) for x in c])
    print("KL seeds used to solve the code:")
    for rank, rotation in seed_rotations:
        print(f"  rank {rank}: R_{N}(U) J_+^{rank} R_{N}(U)^dagger, U = {rotation}")
    print("Define u+i*v = c_6*a_6/c_2 and w+i*z = c_10*a_10/c_2.")
    print("Eliminant equation for z = Im(c_10 a_10 / (c_2 a_2)):", quartic, "= 0")
    print("u =", sp.simplify(sp.re(x6)))
    print("v/z =", sp.simplify(sp.im(x6) / sp.im(x10)))
    print("w =", sp.simplify(sp.re(x10)))
    print("z =", sp.simplify(sp.im(x10)))
    print("Recover v by multiplying the printed v/z by z.")
    print("(a_2, a_6, a_10) = (1, (c_2/c_6)(u+i*v), (c_2/c_10)(w+i*z))")
    print("Exact isometry and reduced seed KL checks passed through rank 5.")
    print("Saved exact basis coefficients to", path)
    return candidate_codewords(coefficients)[:2]


def verify_code(path=CODE_FILE):
    """Check the saved code directly against every tensor through rank 5."""
    trivial, c, V, W = build_basis()
    verify(trivial, c, V, W)
    path = Path(path)
    if path.exists():
        data = json.loads(path.read_text(encoding="utf-8"))
        if (data["n"], data["group"], data["irrep"]) != (N, "2O", "rho_3"):
            raise ValueError("saved code has the wrong n, group, or irrep")
        coefficients = sp.Matrix([sp.sympify(data["coefficients"][f"a_{m}"])
                                  for m in LABELS])
        print("Verifying saved basis coefficients from", path)
    else:
        coefficients, _ = solve_coefficients(c, V, W)
        print("No saved file; verifying freshly derived basis coefficients.")

    # Both columns have the same norm, so raw-vector KL checks need no square root.
    code = sp.Matrix.hstack(V * coefficients, W * coefficients)
    gram = code.H * code
    assert sp.simplify(gram[0, 0]) != 0
    assert sp.simplify(gram[0, 0] - gram[1, 1]) == 0
    assert sp.simplify(gram[0, 1]) == 0
    for rank in range(6):
        for mu, tensor in elementary_tensors(rank).items():
            compressed = code.H * tensor * code
            if any(sp.simplify(entry) != 0 for entry in
                   (compressed[0, 1], compressed[1, 0],
                    compressed[0, 0] - compressed[1, 1])):
                raise AssertionError(f"rank {rank}, component {mu} fails KL")

    rank_six = elementary_tensors(6)[6]
    assert sp.simplify((code[:, 0].H * rank_six * code[:, 1])[0]).is_zero is False
    print("Direct KL checks passed for all 36 tensors of ranks 0-5; rank 6 fails.")


if __name__ == "__main__":
    derive_code()
    verify_code()
