import cvxpy as cp
import numpy as np
import math 
from math import comb
import scipy

from scipy.stats import unitary_group
from scipy.special import factorial, gammaln
from scipy.linalg import expm

from scipy.sparse import csr_matrix, kron as spkron, csc_matrix
from scipy.special import gammaln
import matplotlib.pyplot as plt

from faster_gen_dual_rail import *
import os
import re

from qutip import basis, Qobj, ket2dm
from qutip.wigner import spin_wigner

from sympy.physics import wigner as spw 

import pickle 
import argparse 

import itertools
import json, math, cmath, numpy as np

import sympy 
import time 



paulix = np.array([[0,1],[1,0]])
pauliy = np.array([[0,-1j],[1j,0]])
pauliz = np.array([[1,0],[0,-1]])
I2 = np.array([[1,0],[0,1]])



# Takes partial trace of operator over system A (trace_over = 0)
# or B (trace_over = 1)
def partial_trace(operator, dimA, dimB, trace_over=1):
    if operator.shape != (dimA*dimB, dimA*dimB):
        return "Operator is incompatible with the given dimensions"
    operator_tensor = operator.reshape(dimA, dimB, dimA, dimB)
    return np.trace(operator_tensor, axis1 = trace_over, axis2=trace_over + 2)



def row_reduce_multiplicity_matrices(*matrix_families, tol=1e-10):
    """Return an orthonormal basis for the span of the B-matrices."""
    if not matrix_families:
        raise ValueError("At least one matrix family is required")
    if tol < 0:
        raise ValueError("tol must be nonnegative")

    arrays = [np.asarray(family, dtype=complex) for family in matrix_families]
    matrix_shape = arrays[0].shape[-2:]
    if len(matrix_shape) != 2 or matrix_shape[0] != matrix_shape[1]:
        raise ValueError("Multiplicity matrices must be square")
    if any(array.ndim < 2 or array.shape[-2:] != matrix_shape for array in arrays):
        raise ValueError("All multiplicity matrices must have the same shape")

    matrix_size = matrix_shape[0]
    rows = np.concatenate(
        [array.reshape(-1, matrix_size**2) for array in arrays], axis=0
    )

    row_norms = np.linalg.norm(rows, axis=1)
    nonzero_rows = row_norms > tol
    if not np.any(nonzero_rows):
        return np.empty((0, matrix_size, matrix_size), dtype=complex)

    normalized_rows = rows[nonzero_rows] / row_norms[nonzero_rows, None]
    _, singular_values, vh = np.linalg.svd(
        normalized_rows,
        full_matrices=False,
    )
    numerical_rank = np.count_nonzero(
        singular_values > tol * singular_values[0]
    )

    return vh[:numerical_rank].reshape(-1, matrix_size, matrix_size)


def find_multiplicity_vector(
    exact_B,
    next_B=None,
    restarts=20,
):
    exact_B = np.asarray(exact_B, dtype=complex)
    dimension = exact_B.shape[-1]
    rng = np.random.default_rng()

    exact_B = exact_B.reshape(-1, dimension, dimension)
    if len(exact_B) == 0:
        exact_B = np.zeros((1, dimension, dimension), dtype=complex)
    if next_B is not None:
        next_B = np.asarray(next_B, dtype=complex)
        next_B = next_B.reshape(-1, dimension, dimension)

    hermitian_B = np.concatenate(
        [
            (exact_B + exact_B.conj().transpose(0, 2, 1)) / 2,
            (exact_B - exact_B.conj().transpose(0, 2, 1)) / (2j),
        ]
    )
    coordinates = np.concatenate(
        [
            hermitian_B.real.reshape(len(hermitian_B), -1),
            hermitian_B.imag.reshape(len(hermitian_B), -1),
        ],
        axis=1,
    )
    _, singular_values, vh = np.linalg.svd(coordinates, full_matrices=False)
    rank = np.sum(singular_values > 1e-10 * singular_values[0])
    coordinates = vh[:rank]
    hermitian_B = (
        coordinates[:, : dimension**2]
        + 1j * coordinates[:, dimension**2 :]
    ).reshape(-1, dimension, dimension)

    def unpack(parameters):
        return parameters[:dimension] + 1j * parameters[dimension:]

    def expectations(vector, matrices):
        return np.einsum("i,kij,j->k", vector.conj(), matrices, vector)

    def objective(real_vector):
        if next_B is None:
            return 0.0
        values = expectations(unpack(real_vector), next_B)
        return np.sum(np.abs(values) ** 2)

    def constraints(real_vector):
        vector = unpack(real_vector)
        values = expectations(vector, hermitian_B).real
        return np.r_[
            values,
            np.vdot(vector, vector).real - 1,
        ]

    best_result = None

    for _ in range(restarts):
        initial = (
            rng.normal(size=dimension)
            + 1j * rng.normal(size=dimension)
        )
        initial /= np.linalg.norm(initial)
        initial_parameters = np.r_[initial.real, initial.imag]

        result = scipy.optimize.minimize(
            objective,
            initial_parameters,
            method="SLSQP",
            constraints={"type": "eq", "fun": constraints},
            options={"ftol": 1e-12, "maxiter": 2000},
        )

        if result.success and (
            best_result is None or result.fun < best_result.fun
        ):
            best_result = result

    if best_result is None:
        raise RuntimeError("No constrained optimization succeeded.")

    vector = unpack(best_result.x)
    vector /= np.linalg.norm(vector)
    phase_index = np.argmax(np.abs(vector))
    return vector * np.exp(-1j * np.angle(vector[phase_index]))



class IrrepData:
    def __init__(self, filename, dimension, class_sizes, characters, matrices, generator_indices=None, generator_matrices=None):
        self.filename = filename
        self.dimension = dimension
        self.class_sizes = class_sizes
        self.characters = characters
        self.matrices = matrices
        self.generator_indices = [] if generator_indices is None else generator_indices
        self.generator_matrices = [] if generator_matrices is None else generator_matrices
        self.generators = self.generator_matrices

    @staticmethod
    def E(n): 
        return cmath.exp(2j * math.pi / n)

    @staticmethod
    def eval_cyclo(s):
        """Turn GAP's 'E(n)^k' syntax into a Python complex number."""
        if isinstance(s, (int, float, complex)): 
            return s
        s = s.replace("^", "**")
        return eval(s, {"E": IrrepData.E, "cmath": cmath, "math": math})

    @staticmethod
    def unitarize_rep(mats):
        # 1) Build invariant Hermitian form
        H = sum(M.conj().T @ M for M in mats) / len(mats)
        H = 0.5 * (H + H.conj().T)  # symmetrize

        # 2) Cholesky: H = L L†  (L lower-triangular)
        L = np.linalg.cholesky(H)
        C = np.linalg.inv(L.conj().T)   # C = (L†)^(-1)
        Cinv = L.conj().T               # C^{-1} = L†

        # 3) Transform representation
        U_mats = [Cinv @ M @ C for M in mats]
        return U_mats, C

    @classmethod
    def load(cls, filename):
        """Read GAP irrep JSON (with 'characters' and flattened 'matrices') and convert to numbers."""
        with open(filename) as f:
            data = json.load(f)

        dim = data["dimension"]

        # characters: list of strings -> complex
        chars = [cls.eval_cyclo(x) for x in data["characters"]]

        # matrices: each entry is a flat list of length dim^2 (row-major)
        mats = []
        for flat in data["matrices"]:
            vals = [cls.eval_cyclo(s) for s in flat]
            M = np.array(vals, dtype=complex).reshape((dim, dim))
            mats.append(M)

        gen_mats = []
        for flat in data.get("generator_matrices", []):
            vals = [cls.eval_cyclo(s) for s in flat]
            M = np.array(vals, dtype=complex).reshape((dim, dim))
            gen_mats.append(M)

        # GAP's Position(...) is 1-based; store indices in Python's 0-based convention.
        generator_indices = []
        for idx in data.get("generator_indices", []):
            py_idx = int(idx) - 1
            if py_idx < 0 or py_idx >= len(mats):
                raise ValueError(f"Generator index {idx} is out of range for {filename}")
            generator_indices.append(py_idx)

        # optional: unitarize if you want
        mats, C = cls.unitarize_rep(mats)  # your existing function
        Cinv = np.linalg.inv(C)
        gen_mats = [Cinv @ M @ C for M in gen_mats]
        

        return cls(
            filename=filename,
            dimension=dim,
            class_sizes=data["class_sizes"],
            characters=chars,
            matrices=mats,
            generator_indices=generator_indices,
            generator_matrices=gen_mats,
        )

    def as_dict(self):
        return {
            "dimension": self.dimension,
            "class_sizes": self.class_sizes,
            "characters": self.characters,
            "matrices": self.matrices,
            "generator_indices": self.generator_indices,
            "generator_matrices": self.generator_matrices,
        }

    def frobenius_schur_indicator(self, tol=1e-8):
        """Return the Frobenius--Schur indicator of this representation."""
        if not self.matrices:
            raise ValueError("The representation must contain at least one group element")

        indicator = sum(
            np.trace(matrix @ matrix) for matrix in self.matrices
        ) / len(self.matrices)
        rounded_indicator = int(round(indicator.real))

        if not np.isclose(indicator, rounded_indicator, atol=tol, rtol=0):
            raise ValueError(
                "The computed Frobenius--Schur indicator is not an integer "
                f"within tolerance: {indicator}"
            )

        return rounded_indicator

    @classmethod
    def direct_sum(cls, *irreps):
        """Return the block-diagonal direct sum of compatible representations."""
        if len(irreps) == 1 and isinstance(irreps[0], (list, tuple)):
            irreps = tuple(irreps[0])
        if len(irreps) < 2:
            raise ValueError("A direct sum requires at least two representations")
        if not all(isinstance(irrep, cls) for irrep in irreps):
            raise TypeError("Each summand must be an IrrepData object")

        reference = irreps[0]
        for irrep in irreps[1:]:
            if len(irrep.matrices) != len(reference.matrices):
                raise ValueError("Direct-sum representations must use the same group elements")
            if irrep.class_sizes != reference.class_sizes:
                raise ValueError("Direct-sum representations must use the same conjugacy-class ordering")
            if irrep.generator_indices != reference.generator_indices:
                raise ValueError("Direct-sum representations must use the same generators")

        matrices = [
            scipy.linalg.block_diag(*(irrep.matrices[i] for irrep in irreps))
            for i in range(len(reference.matrices))
        ]
        characters = [
            sum(irrep.characters[i] for irrep in irreps)
            for i in range(len(reference.characters))
        ]
        filename = "direct_sum(" + ",".join(
            os.path.basename(irrep.filename) for irrep in irreps
        ) + ")"

        return cls(
            filename=filename,
            dimension=sum(irrep.dimension for irrep in irreps),
            class_sizes=list(reference.class_sizes),
            characters=characters,
            matrices=matrices,
            generator_indices=list(reference.generator_indices),
            generator_matrices=[
                matrices[index] for index in reference.generator_indices
            ],
        )


class FiniteGroupData:
    def __init__(self, folder):
        self.folder = folder
        self.irrep_files = self.list_irrep_files(folder)
        self.irreps = {}
        for filename in self.irrep_files:
            self.irreps[os.path.basename(filename)] = IrrepData.load(filename)

    @staticmethod
    def list_irrep_files(folder):
        return sorted(
            os.path.join(folder, name)
            for name in os.listdir(folder)
            if name.endswith(".json") and name.startswith("irrep_dim")
        )

    @staticmethod
    def character_inner_product(char1, char2):
        # compute the inner product of two characters
        if len(char1) != len(char2):
            raise ValueError("Characters must have the same length")
        
        size = len(char1)
        inner_prod = sum(c1 * np.conj(c2) for c1, c2 in zip(char1, char2)) / size
        return np.round(np.real_if_close(inner_prod)).astype(int)

    @staticmethod
    def end_representation_character(char):
        # given a character, compute the character of the end representation
        return [abs(c)**2 for c in char]

    @staticmethod
    def find_irrep_multiplicity_from_characters(irrep_char, representation_char):
        # given the character of an irrep, and the character of a reducible representation, find the multiplicity     

        return FiniteGroupData.character_inner_product(irrep_char, representation_char)

    @staticmethod
    def find_irrep_decomposition_from_characters(irrep_chars, representation_char):
        # given a list of irrep characters, and the character of a reducible representation, find the multiplicity of each irrep in the decomposition
        multiplicities = {}
        for name, irrep_char in irrep_chars.items():
            m = FiniteGroupData.find_irrep_multiplicity_from_characters(irrep_char, representation_char)
            multiplicities[name] = m

        return multiplicities

    def get_all_characters(self):
        all_irrep_chars = {}
        for filename, data in self.irreps.items():
            all_irrep_chars[filename] = data.characters

        return all_irrep_chars

    def get_irrep(self, irrep):
        if isinstance(irrep, IrrepData):
            return irrep

        name = os.path.basename(irrep)
        if name not in self.irreps:
            raise KeyError(f"Irrep {irrep} was not found in {self.folder}")

        return self.irreps[name]


class SpinSector:
    def __init__(self, n):
        self.n = n
        self.dim = int(n + 1)
        self.j = n / 2
        self.Jz = self.J_z()
        self.Jp = self.J_plus()
        self.Jm = self.J_minus()
        self.Jx = self.J_x()
        self.Jy = self.J_y()
        self._irreducible_tensors = {}

    # Generate the J_z rotation matrix in the Schwinger model
    # This is clearly diagonal
    def Jz_rotation_matrix(self, theta):
        n = self.n
        matrix = np.zeros((self.dim, self.dim), dtype=complex)

        for k in range(n+1):
            phase = np.exp(-1j*(theta/2)*(2*k - n))
            matrix[k, k] = phase

        return matrix

    def J_z(self):
        n = self.n
        matrix = np.zeros((self.dim, self.dim), dtype=complex)
        for k in range(n+1):
            diag =k - n/2
            matrix[k,k] = diag

        return matrix 

    def J_plus(self):
        n = self.n
        dim = self.dim

        Jp = np.zeros((dim, dim), dtype=complex)
        for k in range(n):
            Jp[k+1, k] = np.sqrt((k+1)*(n-k))


        return Jp

    def J_minus(self):
        n = self.n
        dim = self.dim

        Jm = np.zeros((dim, dim), dtype=complex)
        for k in range(1,n+1):
            Jm[k-1, k] = np.sqrt(k*(n-k+1))

        return Jm

    def J_x(self):
        return 0.5*(self.J_plus() + self.J_minus())

    def J_y(self):
        return -0.5j*(self.J_plus() - self.J_minus())

    def sym_n_rep_from_2d(self, U):
        """Compute the induced unitary action on the fixed-n photon sector.

        Convention used throughout this file and faster_gen_dual_rail.py:
            basis index k represents |k, n-k> = |n_a=k, n_b=n-k>.
            U is expressed in the ordered one-photon mode basis (|a>, |b>).

        Thus D[l, k] = <l, n-l| Sym^n(U) |k, n-k>.
        """
        U = np.asarray(U, dtype=complex)
        if U.shape != (2, 2):
            raise ValueError("U must have shape (2, 2)")

        unitary_residual = np.linalg.norm(
            U.conj().T @ U - np.eye(2),
            ord=2,
        )
        if unitary_residual > 1e-10:
            raise ValueError(
                "The symmetric-power construction requires a unitary "
                f"2x2 matrix; residual={unitary_residual:.3e}"
            )

        generator = scipy.linalg.logm(U)
        generator = 0.5 * (generator - generator.conj().T)

        reconstruction_residual = np.linalg.norm(
            scipy.linalg.expm(generator) - U,
            ord=2,
        )
        if reconstruction_residual > 1e-10:
            raise RuntimeError(
                "Could not construct a reliable anti-Hermitian logarithm "
                f"of U; residual={reconstruction_residual:.3e}"
            )

        identity = np.eye(self.dim, dtype=complex)
        number_a = (self.n / 2) * identity + self.Jz
        number_b = (self.n / 2) * identity - self.Jz

        # Use fresh ladder matrices so this construction does not depend on
        # callers preserving the cached matrices in their original scale.
        Jp = self.J_plus()
        Jm = self.J_minus()

        physical_generator = (
            generator[0, 0] * number_a
            + generator[1, 1] * number_b
            + generator[0, 1] * Jp
            + generator[1, 0] * Jm
        )
        return scipy.linalg.expm(physical_generator)

    def full_symmetric_n_rep_from_2d(self, matrices, T=None):
        """Compute Sym^n(U) for every fundamental representation matrix U.

        The optional T argument is accepted for compatibility with the previous
        tensor-product implementation and is intentionally unused.
        """
        return [self.sym_n_rep_from_2d(U) for U in matrices]
    
    def irreducible_tensors(self, rank=1):
        """Generate the trace-orthonormal irreducible tensors of a given rank."""
        if not isinstance(rank, (int, np.integer)):
            raise TypeError("rank must be an integer")
        if rank < 0 or rank > self.n:
            raise ValueError("rank must satisfy 0 <= rank <= n")
        if rank in self._irreducible_tensors:
            return self._irreducible_tensors[rank]

        j = sympy.Rational(self.n, 2)
        rank_tensors = {}

        for q in range(rank, -rank - 1, -1):
            tensor = np.zeros((self.dim, self.dim), dtype=complex)

            for row in range(self.dim):
                col = row - q
                if not 0 <= col < self.dim:
                    continue

                m_row = -j + row
                m_col = -j + col
                phase = (-1) ** int(j - m_row)
                tensor[row, col] = (
                    np.sqrt(2 * rank + 1)
                    * phase
                    * complex(
                        spw.wigner_3j(
                            j,
                            rank,
                            j,
                            -m_row,
                            q,
                            m_col,
                        )
                    )
                )

            tensor /= np.linalg.norm(tensor)
            rank_tensors[q] = tensor

        # Match the existing convention in which the highest-weight tensor is
        # the normalized positive multiple of J_+^rank.
        highest_weight = np.linalg.matrix_power(self.J_plus(), rank)
        overlap = np.vdot(rank_tensors[rank], highest_weight)
        common_phase = overlap / abs(overlap)
        rank_tensors = {
            q: common_phase * tensor
            for q, tensor in rank_tensors.items()
        }

        self._irreducible_tensors[rank] = rank_tensors
        return rank_tensors


class IsotypicComponent:
    def __init__(self, n, fundamental_irrep, target_irrep, tol=1e-6, spin_sector=None, group_representatives=None):
        self.n = n
        self.tol = tol
        self.spin_sector = SpinSector(n) if spin_sector is None else spin_sector
        self.fundamental_irrep = self._coerce_irrep(fundamental_irrep)
        self.target_irrep = self._coerce_irrep(target_irrep)
        self.dimension = self.target_irrep.dimension
        self.irrep_dimension = self.dimension
        if group_representatives is None:
            self.group_representatives = self.symmetric_representations()
        else:
            self.group_representatives = group_representatives
            if len(self.group_representatives) != len(self.target_irrep.matrices):
                raise ValueError("Physical representation and target irrep must use the same group element ordering")
        self.representations = self.group_representatives
        self.multiplicity_raw = self.irrep_multiplicity(
            self.target_irrep.characters,
            self.group_representatives,
        )
        self.multiplicity = int(round(float(np.real(np.real_if_close(self.multiplicity_raw)))))
        self.rank = self.dimension * self.multiplicity
        
        self.projector_matrix = self.projector_from_irrep(
            self.target_irrep.characters,
            self.group_representatives,
            self.target_irrep.dimension,
        )
        self.projector_eigvals, self.projector_eigvecs = np.linalg.eigh(self.projector_matrix)
        self.basis = self.projector_eigvecs[:, self.projector_eigvals > 1 - self.tol]
        self.projector_rank = self.basis.shape[1]
        self.embedding_operator = self.isotypic_embedding()



    @staticmethod
    def _coerce_irrep(irrep):
        if isinstance(irrep, IrrepData):
            return irrep
        return IrrepData.load(irrep)

    @staticmethod
    def irrep_multiplicity(irrep_char, representation):
        # given the character of an irrep, and a reducible representation, find the multiplicity 
        m = 0
        size = len(irrep_char)

        for i in range(size):
            char1 = irrep_char[i]
            char2 = np.trace(representation[i])
            m += (1 / size) * np.conj(char1) * char2

        return m

    @staticmethod
    def projector_from_irrep(characters, representations, dim):
        group_size = len(characters)
        proj = np.zeros_like(representations[0], dtype=complex)
        for i in range(group_size):
            proj += (dim/group_size)*np.conj(characters[i])*representations[i]

        return proj
    
    @staticmethod
    def tr_inner_product(A,B):
        dim = A.shape[1]
        return 1/dim*np.linalg.trace(A.T.conj()@B)

    @staticmethod
    def devectorize_operator_columns(vectorized_solutions, target_dimensions):
        rows, cols = target_dimensions
        if vectorized_solutions.ndim != 2:
            raise ValueError("vectorized_solutions must be a matrix")
        if vectorized_solutions.shape[0] != rows * cols:
            raise ValueError("Vectorized solution length does not match target dimensions")

        operators = []
        for i in range(vectorized_solutions.shape[1]):
            vector = vectorized_solutions[:, i]
            operator = vector.reshape((rows, cols), order='F')
            operators.append(operator)

        return operators

    @staticmethod
    def orthonormalize_operators(operators):
        # given a list of operators, 
        # make them orthonormal wrt 1/d Tr(A^\dagger B)
        size = len(operators)
        Gamma_matrix = np.zeros((size, size), dtype=complex)

        for a in range(size):
            for b in range(size):
                opA = operators[a]
                opB = operators[b]
                Gamma_matrix[a][b] = IsotypicComponent.tr_inner_product(opA, opB)


        Gamma_matrix = 0.5 * (Gamma_matrix + Gamma_matrix.conj().T)
        eigenvalues, eigenvectors = np.linalg.eigh(Gamma_matrix)
        threshold = (
            np.finfo(float).eps
            * Gamma_matrix.shape[0]
            * max(1.0, np.max(np.abs(eigenvalues)))
        )
        if np.min(eigenvalues) <= threshold:
            raise ValueError("Operator Gram matrix is not positive definite")

        Gamma_inv_sqrt = (
            eigenvectors
            @ np.diag(1 / np.sqrt(eigenvalues))
            @ eigenvectors.conj().T
        )

        new_operators = []

        for a in range(size):
            new_op_A = np.zeros(operators[0].shape, dtype=complex)
            for b in range(size):
                opB = operators[b]
                new_op_A += opB*Gamma_inv_sqrt[b][a]
        
            new_operators.append(new_op_A)

        return new_operators


    def symmetric_representations(self):
        return self.spin_sector.full_symmetric_n_rep_from_2d(self.fundamental_irrep.matrices)

    def projector(self):
        return self.projector_matrix

    def optimize_isotypic_component(self, d, eta, rounds, file_output, tol=1e-6):
        return optimize_isotypic_component(self.projector_matrix, n=self.n, d=d, eta=eta, rounds=rounds, file_output=file_output, tol=tol)

    def isotypic_embedding(self):
        # The goal is to construct the operator S_C derived in the notes 
        equation_matrices = []
        for generator_index in self.target_irrep.generator_indices:
            elementary_generator = self.target_irrep.matrices[generator_index]
            reducible_generator = self.group_representatives[generator_index]

            equation_matrix = np.kron(np.identity(self.dimension), reducible_generator) - np.kron(elementary_generator.T, np.identity(self.n + 1))
            equation_matrices.append(equation_matrix)

        final_vector_equation =  np.vstack(equation_matrices)
        kernel = scipy.linalg.null_space(final_vector_equation, rcond=self.tol)

        operators = IsotypicComponent.devectorize_operator_columns(
            kernel,
            (self.n + 1, self.dimension),
        )

        operators = IsotypicComponent.orthonormalize_operators(operators)
        return np.hstack(operators)


class GroupQECModel:
    DEFAULT_FUNDAMENTAL_IRREPS = {
        "2I": "irrep_dim2_2I_2.json",
        "2O": "irrep_dim2_2O_4.json",
        "2T": "irrep_dim2_2T_4.json",
        "C4": ("irrep_dim1_C4_3.json", "irrep_dim1_C4_4.json"),
        "C6": ("irrep_dim1_C6_4.json", "irrep_dim1_C6_6.json"),
        "C8": ("irrep_dim1_C8_5.json", "irrep_dim1_C8_8.json"),
        "C10": ("irrep_dim1_C10_4.json", "irrep_dim1_C10_10.json"),
        '2D2': 'irrep_dim2_2D2_5.json',
        '2D3': 'irrep_dim2_2D3_6.json',
        '2D4': 'irrep_dim2_2D4_6.json',
        '2D5': 'irrep_dim2_2D5_6.json',
        '2D6': 'irrep_dim2_2D6_8.json'
    }

    def __init__(self, group_folder, n, fundamental_irrep=None, tol=1e-6):
        self.group_folder = group_folder
        self.group_name = os.path.basename(os.path.normpath(group_folder))
        self.n = n
        self.tol = tol

        self.group = FiniteGroupData(group_folder)
        self.irreps = self.group.irreps
        self.spin_sector = SpinSector(n)
        self.fundamental_irrep = self._choose_fundamental_irrep(fundamental_irrep)
        self.generator_indices = self.fundamental_irrep.generator_indices

        self.physical_representations = self.spin_sector.full_symmetric_n_rep_from_2d(
            self.fundamental_irrep.matrices
        )
        self.group_representatives = self.physical_representations
        self.physical_character = [np.trace(M) for M in self.physical_representations]
        self.decomposition = self.decompose()

        self.components = {}
        self.interwiners = {}


    def _choose_fundamental_irrep(self, fundamental_irrep):
        if fundamental_irrep is None and self.group_name in self.DEFAULT_FUNDAMENTAL_IRREPS:
            fundamental_irrep = self.DEFAULT_FUNDAMENTAL_IRREPS[self.group_name]

        if isinstance(fundamental_irrep, (list, tuple)):
            if len(fundamental_irrep) != 2:
                raise ValueError("The physical fundamental representation must use exactly two 1D irreps")
            summands = [self.get_irrep(irrep) for irrep in fundamental_irrep]
            if any(irrep.dimension != 1 for irrep in summands):
                raise ValueError("A fundamental-irrep pair must contain two 1D irreps")
            irrep = IrrepData.direct_sum(summands)
        elif fundamental_irrep is not None:
            irrep = self.get_irrep(fundamental_irrep)
        else:
            dim_two_irreps = [
                irrep
                for irrep in self.irreps.values()
                if irrep.dimension == 2
            ]
            if len(dim_two_irreps) != 1:
                raise ValueError("Specify fundamental_irrep when the group does not have a unique 2D irrep")
            irrep = dim_two_irreps[0]

        if irrep.dimension != 2:
            raise ValueError("The physical fundamental irrep must be 2-dimensional")

        return irrep

    def _irrep_key(self, irrep):
        if isinstance(irrep, IrrepData):
            return os.path.basename(irrep.filename)
        return os.path.basename(irrep)

    def get_irrep(self, irrep):
        return self.group.get_irrep(irrep)

    def irrep_names(self, dimension=None):
        names = sorted(self.irreps)
        if dimension is None:
            return names

        return [
            name
            for name in names
            if self.irreps[name].dimension == dimension
        ]

    def irrep_multiplicity(self, target_irrep):
        target = self.get_irrep(target_irrep)
        return IsotypicComponent.irrep_multiplicity(
            target.characters,
            self.physical_representations,
        )

    def decompose(self, tol=None):
        if tol is None:
            tol = self.tol

        decomposition = {}
        for name, irrep in self.irreps.items():
            multiplicity = self.irrep_multiplicity(irrep)
            if abs(multiplicity) <= tol:
                continue

            multiplicity = np.real_if_close(multiplicity)
            if np.isclose(multiplicity, round(float(np.real(multiplicity))), atol=tol):
                multiplicity = int(round(float(np.real(multiplicity))))

            decomposition[name] = multiplicity

        return decomposition

    def isotypic_component(self, target_irrep):
        target = self.get_irrep(target_irrep)
        key = self._irrep_key(target)
        if key not in self.components:
            self.components[key] = IsotypicComponent(
                n=self.n,
                fundamental_irrep=self.fundamental_irrep,
                target_irrep=target,
                tol=self.tol,
                spin_sector=self.spin_sector,
                group_representatives=self.physical_representations,
            )

        return self.components[key]

    def projector(self, target_irrep):
        return self.isotypic_component(target_irrep).projector()

    def embedding(self, target_irrep):
        return self.isotypic_component(target_irrep).embedding_operator

    def code_from_multiplicity_vector(self, target_irrep, multiplicity_vector):
        """Build the fixed-n code selected by a normalized multiplicity vector."""
        component = self.isotypic_component(target_irrep)
        vector = np.asarray(multiplicity_vector, dtype=complex).reshape(-1)
        if vector.size != component.multiplicity:
            raise ValueError(
                f"Multiplicity vector must have length {component.multiplicity}"
            )
        if not np.isclose(np.linalg.norm(vector), 1, atol=self.tol):
            raise ValueError("Multiplicity vector must be normalized")

        isometry = component.embedding_operator @ np.kron(
            vector[:, None], np.eye(component.dimension)
        )
        return RandomConstantNCode(n=self.n, d=component.dimension, U=isometry)

    def maximum_kl_residual(self, code, maximum_rank):
        """Return the maximum normalized KL residual through ``maximum_rank``.

        Each irreducible tensor is Hilbert--Schmidt normalized, and the residual
        is the Frobenius norm of the traceless part of its compression to the
        code space.
        """
        if not isinstance(maximum_rank, (int, np.integer)):
            raise TypeError("maximum_rank must be an integer")
        if maximum_rank < 0 or maximum_rank > self.n:
            raise ValueError("maximum_rank must satisfy 0 <= maximum_rank <= n")

        encoding = np.asarray(code.U, dtype=np.complex128)
        if encoding.ndim != 2 or encoding.shape[0] != self.n + 1:
            raise ValueError("Code isometry is incompatible with this spin sector")

        logical_dimension = encoding.shape[1]
        logical_identity = np.eye(logical_dimension, dtype=np.complex128)
        maximum_residual = 0.0
        for rank in range(maximum_rank + 1):
            for tensor in self.spin_sector.irreducible_tensors(rank=rank).values():
                compressed_tensor = encoding.conj().T @ tensor @ encoding
                traceless_part = compressed_tensor - (
                    np.trace(compressed_tensor)
                    * logical_identity
                    / logical_dimension
                )
                maximum_residual = max(
                    maximum_residual,
                    float(np.linalg.norm(traceless_part)),
                )

        return maximum_residual

    def tensor_product_decomposition(self, irrep_i, irrep_j):
        """Decompose ``rho_i tensor rho_j^*`` into irreducible representations.

        The irreps may be specified either by name/path or by ``IrrepData``
        objects.  The returned dictionary maps irrep filenames to their
        (positive) multiplicities, matching :meth:`end_decomposition`.
        """
        rho_i = self.get_irrep(irrep_i)
        rho_j = self.get_irrep(irrep_j)
        product_char = [
            char_i * np.conj(char_j)
            for char_i, char_j in zip(rho_i.characters, rho_j.characters)
        ]
        return {
            name: mult
            for name, rep in self.irreps.items()
            if (
                mult := FiniteGroupData.character_inner_product(
                    product_char,
                    rep.characters,
                )
            ) > 0
        }
    
    def end_decomposition(self, end_irrep):
        return self.tensor_product_decomposition(end_irrep, end_irrep)

    def rank_seeds(self, rank, random_state=0, max_attempts=10_000):
        """Generate the rank-``rank`` seed operators from Theorem 11.

        The returned list contains the minimum number

            ceil((2 * rank + 1) / |G_tilde|)

        of operators whose G-orbits span the irreducible tensor space of the
        requested rank.  The rank-zero seed is defined to be the identity.
        """
        if not isinstance(rank, (int, np.integer)):
            raise TypeError("rank must be an integer")
        if rank < 0 or rank > self.n:
            raise ValueError("rank must satisfy 0 <= rank <= n")
        if not isinstance(max_attempts, (int, np.integer)) or max_attempts <= 0:
            raise ValueError("max_attempts must be a positive integer")

        if rank == 0:
            return [np.eye(self.n + 1, dtype=complex)]

        def same_projective_vector(vector_1, vector_2):
            overlap = abs(np.vdot(vector_1, vector_2))
            overlap /= np.linalg.norm(vector_1) * np.linalg.norm(vector_2)
            return overlap >= 1 - self.tol

        def same_projective_matrix(matrix_1, matrix_2):
            overlap = abs(np.vdot(matrix_1, matrix_2))
            overlap /= np.linalg.norm(matrix_1) * np.linalg.norm(matrix_2)
            return overlap >= 1 - self.tol

        projective_group = []
        for group_element in self.fundamental_irrep.matrices:
            already_present = any(
                same_projective_matrix(group_element, representative)
                for representative in projective_group
            )
            if not already_present:
                projective_group.append(group_element)

        projective_order = len(projective_group)
        number_of_seeds = math.ceil((2 * rank + 1) / projective_order)

        rng = np.random.default_rng(random_state)
        seed_points = []
        occupied_orbit_points = []

        for _ in range(max_attempts):
            point = rng.normal(size=2) + 1j * rng.normal(size=2)
            point /= np.linalg.norm(point)

            orbit = []
            for group_element in projective_group:
                orbit_point = group_element @ point
                orbit_point /= np.linalg.norm(orbit_point)

                if not any(
                    same_projective_vector(orbit_point, previous)
                    for previous in orbit
                ):
                    orbit.append(orbit_point)

            if len(orbit) != projective_order:
                continue

            intersects_previous_orbit = any(
                same_projective_vector(orbit_point, previous)
                for orbit_point in orbit
                for previous in occupied_orbit_points
            )
            if intersects_previous_orbit:
                continue

            seed_points.append(point)
            occupied_orbit_points.extend(orbit)

            if len(seed_points) == number_of_seeds:
                break

        if len(seed_points) != number_of_seeds:
            raise RuntimeError(
                "Could not find enough free, pairwise disjoint projective orbits"
            )

        highest_weight = np.linalg.matrix_power(
            self.spin_sector.Jp,
            rank,
        ).copy()
        highest_weight /= np.linalg.norm(highest_weight)
        seeds = []

        for a, b in seed_points:
            rotation = np.array(
                [
                    [a, -np.conj(b)],
                    [b, np.conj(a)],
                ],
                dtype=complex,
            )
            physical_rotation = self.spin_sector.sym_n_rep_from_2d(rotation)
            seed = (
                physical_rotation
                @ highest_weight
                @ physical_rotation.conj().T
            )
            seeds.append(seed)

        orbit_columns = []
        for seed in seeds:
            for group_element in self.physical_representations:
                orbit_operator = (
                    group_element
                    @ seed
                    @ group_element.conj().T
                )
                orbit_columns.append(orbit_operator.reshape(-1))

        orbit_matrix = np.column_stack(orbit_columns)
        singular_values = np.linalg.svd(orbit_matrix, compute_uv=False)
        numerical_rank = np.count_nonzero(
            singular_values > self.tol * singular_values[0]
        )
        expected_rank = 2 * rank + 1

        if numerical_rank != expected_rank:
            raise RuntimeError(
                "The generated seed orbits have numerical rank "
                f"{numerical_rank}; expected {expected_rank}"
            )

        return seeds

    def optimize_multiplicity_code(
        self,
        code_irreps,
        restarts=50,
        tol=1e-8,
    ):
        """Construct a code from rank-one projectors on multiplicity spaces.

        Minimal orbit seeds are used to impose every exactly feasible tensor
        rank. At the first infeasible rank, the residual is minimized over the
        complete trace-orthonormal tensor basis at that rank.
        """
        code_irreps = [
            self._irrep_key(self.get_irrep(C))
            for C in code_irreps
        ]
        components = {
            C: self.isotypic_component(C)
            for C in code_irreps
        }
        mult = {C: components[C].multiplicity for C in code_irreps}
        dim = {C: components[C].dimension for C in code_irreps}

        # Locate each M_C inside the direct sum of multiplicity spaces.
        slices = {}
        start = 0
        for C in code_irreps:
            slices[C] = slice(start, start + mult[C])
            start += mult[C]
        total_mult = start

        # Embed B : M_D -> M_C into its (C,D) block.
        def embed(B, C, D):
            result = np.zeros((total_mult, total_mult), dtype=complex)
            result[slices[C], slices[D]] = B
            return result

        trivial = next(
            name
            for name, irrep in self.irreps.items()
            if irrep.dimension == 1
            and np.allclose(irrep.characters, 1, atol=tol)
        )

        # Convert the B-matrices of a set of seeds into global constraints.
        def constraints_from_seeds(seeds):
            constraints = []
            reference = code_irreps[0]

            for operator in seeds:
                trivial_B = {}

                for C in code_irreps:
                    for D in code_irreps:
                        for U in self.tensor_product_decomposition(C, D):
                            family = self.get_multiplicity_matrices(
                                operator=operator,
                                input_irrep=U,
                                output_irrep=C,
                                source_irrep=D,
                            )

                            if C == D and U == trivial:
                                F = self.intertwiner_basis(trivial, C, C)[0, 0]
                                kappa = np.trace(F) / dim[C]
                                trivial_B[C] = kappa * family[0, 0]
                            else:
                                for B in family.reshape(-1, mult[C], mult[D]):
                                    constraints.append(embed(B, C, D))

                # Enforce the same trivial-sector KL constant for every C.
                for C in code_irreps[1:]:
                    constraints.append(
                        embed(trivial_B[C], C, C)
                        - embed(trivial_B[reference], reference, reference)
                    )

            return np.asarray(constraints).reshape(
                -1, total_mult, total_mult
            )

        reference = code_irreps[0]

        # Ensure that every block of the combined vector has the same norm.
        exact_B = [
            embed(np.eye(mult[C]), C, C)
            - embed(np.eye(mult[reference]), reference, reference)
            for C in code_irreps[1:]
        ]
        if not exact_B:
            exact_B = [np.zeros((total_mult, total_mult), dtype=complex)]
        exact_B = np.asarray(exact_B)

        exact_rank = 0
        vector = None

        # Use the existing seed construction and solve rank by rank.
        for rank in range(1, self.n + 1):
            rank_B = constraints_from_seeds(self.rank_seeds(rank))
            trial_B = row_reduce_multiplicity_matrices(
                exact_B, rank_B, tol=tol
            )

            try:
                vector = find_multiplicity_vector(
                    trial_B,
                    restarts=restarts,
                )
            except RuntimeError:
                full_rank_B = constraints_from_seeds(
                    list(
                        self.spin_sector.irreducible_tensors(
                            rank=rank,
                        ).values()
                    )
                )
                vector = find_multiplicity_vector(
                    exact_B,
                    full_rank_B,
                    restarts=restarts,
                )
                break

            exact_B = trial_B
            exact_rank = rank

        # Split the combined vector into one normalized vector on each M_C.
        vectors = {}
        for C in code_irreps:
            vectors[C] = vector[slices[C]]
            vectors[C] /= np.linalg.norm(vectors[C])

        # Embed one selected copy of each requested irrep into physical space.
        encoding = np.hstack([
            components[C].embedding_operator
            @ np.kron(vectors[C][:, None], np.eye(dim[C]))
            for C in code_irreps
        ])

        code = RandomConstantNCode(
            n=self.n,
            d=encoding.shape[1],
            U=encoding,
        )
        maximum_residual = self.maximum_kl_residual(code, exact_rank)
        print(
            f"KL verification for n={self.n} {self.group_name} code through "
            f"rank {exact_rank}: maximum normalized residual = "
            f"{maximum_residual:.3e}"
        )
        return code, vectors, exact_rank


    

    def intertwiner_basis(self, input_irrep, output_irrep, source_irrep):
        """Construct the F_{j, mu}^{U; C,D} basis from the notes.

        Here ``input_irrep`` is U, ``output_irrep`` is C, and
        ``source_irrep`` is D, so each F_{j, mu}^{U; C,D} lies in Hom(D, C).
        """
        input_irrep = self.get_irrep(input_irrep)
        output_irrep = self.get_irrep(output_irrep)
        source_irrep = self.get_irrep(source_irrep)

        key = (
            self._irrep_key(input_irrep),
            self._irrep_key(output_irrep),
            self._irrep_key(source_irrep),
        )
        if key in self.interwiners:
            return self.interwiners[key]

        decomposition = self.tensor_product_decomposition(
            output_irrep,
            source_irrep,
        )
        if self._irrep_key(input_irrep) not in decomposition:
            return None

        matrix_equations = []
        for index in self.generator_indices:
            input_representative = input_irrep.matrices[index]
            output_representative = output_irrep.matrices[index]
            source_representative = source_irrep.matrices[index]

            hom_representative = np.kron(
                source_representative.conj(),
                output_representative,
            )
            matrix = (
                np.kron(np.identity(input_irrep.dimension), hom_representative)
                - np.kron(
                    input_representative.T,
                    np.identity(
                        output_irrep.dimension * source_irrep.dimension
                    ),
                )
            )
            matrix_equations.append(matrix)

        full_matrix = np.vstack(matrix_equations)
        _, singular_values, vh = scipy.linalg.svd(full_matrix, full_matrices=True)
        numerical_rank = np.count_nonzero(singular_values > self.tol)
        kernel = vh[numerical_rank:].conj().T

        operators = IsotypicComponent.devectorize_operator_columns(
            kernel,
            (
                output_irrep.dimension * source_irrep.dimension,
                input_irrep.dimension,
            ),
        )
        operators = IsotypicComponent.orthonormalize_operators(operators)

        F_basis = np.array([
            IsotypicComponent.devectorize_operator_columns(
                intertwiner,
                (output_irrep.dimension, source_irrep.dimension),
            )
            for intertwiner in operators
        ])

        self.interwiners[key] = F_basis
        return F_basis
    

        
    def get_multiplicity_matrices(
        self,
        operator,
        input_irrep,
        output_irrep,
        source_irrep=None,
    ):
        """Return the matrices B_{j,mu}^{U;C,D}(O).

        Here ``input_irrep`` is U, ``output_irrep`` is C, and
        ``source_irrep`` is D. If ``source_irrep`` is omitted, D is taken
        to be C.
        """
        input_irrep = self.get_irrep(input_irrep)
        output_irrep = self.get_irrep(output_irrep)
        if source_irrep is None:
            source_irrep = output_irrep
        else:
            source_irrep = self.get_irrep(source_irrep)

        for name, irrep in [
            ("output", output_irrep),
            ("source", source_irrep),
        ]:
            key = self._irrep_key(irrep)
            if key not in self.decomposition or self.decomposition[key] == 0:
                raise ValueError(
                    f"The {name} irrep does not occur in the physical "
                    "representation"
                )

        F_basis = self.intertwiner_basis(
            input_irrep=input_irrep,
            output_irrep=output_irrep,
            source_irrep=source_irrep,
        )
        if F_basis is None:
            raise ValueError(
                "The input irrep does not occur in "
                "output_irrep tensor source_irrep^*"
            )

        output_component = self.isotypic_component(output_irrep)
        source_component = self.isotypic_component(source_irrep)

        operator = np.asarray(operator, dtype=complex)
        expected_shape = (self.n + 1, self.n + 1)
        if operator.shape != expected_shape:
            raise ValueError(f"operator must have shape {expected_shape}")

        compressed_operator = (
            output_component.embedding_operator.conj().T
            @ operator
            @ source_component.embedding_operator
        )
        blocks = compressed_operator.reshape(
            output_component.multiplicity,
            output_component.dimension,
            source_component.multiplicity,
            source_component.dimension,
        )

        return np.einsum(
            "jmrs,arbs->jmab",
            F_basis.conj(),
            blocks,
        )
        

def get_code_performance(code, eta,file=None):
    """Evaluate a fixed RandomConstantNCode using its optimal decoder."""
    loss = TruncatedLoss(code.n, eta)
    optimizer = SDPCodeOptimizer(
        loss, T_E_init=code.superoperator, n=code.n, d=code.d, eps=1e-9, file_output=file
    )
    optimizer.solve_decoder_problem()
    optimizer.calculate_performance()
    if file is not None:
        optimizer.save()

    return optimizer.Fe, optimizer.I_c


def sym_n_rep_from_2d(U, n):
    return SpinSector(n).sym_n_rep_from_2d(U)

def full_symmetric_n_rep_from_2d(matrices, n, T=None):
    return SpinSector(n).full_symmetric_n_rep_from_2d(matrices, T=T)


def projector_from_irrep_file(n, fundamental_irrep, target_irrep ):
    component = IsotypicComponent(n, fundamental_irrep, target_irrep)
    return component.projector()




def projector(logical_0, logical_1=None):
    # Accept either two logical vectors or a list of logical vectors
    if logical_1 is None and isinstance(logical_0, (list, tuple)):
        logicals = list(logical_0)
    elif logical_1 is None:
        raise ValueError('projector requires either two logical vectors or a list of logical vectors')
    else:
        logicals = [logical_0, logical_1]

    P = sum(np.outer(l, l.conj()) for l in logicals)
    return P






def optimize_isotypic_component(projector, n, d, eta,rounds,  file_output, tol=1e-6):
    eigvals, eigvecs= np.linalg.eigh(projector)
    keep = eigvals > 1-tol
    V = eigvecs[:, keep]
    
    # print dimension of space 
    print(f"Projector has rank {V.shape[1]}")

    myCode = RandomConstantNCode(n=n, d=d)
    T_E_init = myCode.superoperator    

    LossOp = TruncatedLoss(n, eta)
    optimizer = SDPCodeOptimizer(LossOp, T_E_init=T_E_init, n=n, d=d, encoding_subspace=V, file_output=file_output, save_progress=True)
    optimizer.run(rounds=rounds)
    optimizer.calculate_performance()
    optimizer.save()

    return optimizer.Fe, optimizer.I_c






def get_optimal_isotypic_code(fundamental_irrep_file, code_irrep_file,file_output=None, n = 10, d= 2, eta = 0.9, rounds=5):
    component = IsotypicComponent(n, fundamental_irrep_file, code_irrep_file)
    proj = component.projector()
    if file_output is None:
        group = fundamental_irrep_file.split("_")[2]
        irrep = code_irrep_file.split("_")[3][0]
        file_output = f"{group}_irrep_{irrep}_code_optimization_results_d{d}_n{n}_{rounds}.pkl"

    print(file_output)
    optimize_isotypic_component(proj, n=n, d=d, eta=eta, rounds=rounds, file_output=file_output)





def __main__():
    parser = argparse.ArgumentParser(description='Optimize isotypic component for quantum error correction.')
    parser.add_argument('--fundamental_irrep_file', type=str, required=True, help='Path to the JSON file containing the fundamental irrep data.')
    parser.add_argument('--code_irrep_file', type=str, required=True, help='Path to the JSON file containing the code irrep data.')
    parser.add_argument('--n', type=int, default=10, help='Number of photons (default: 10)')
    parser.add_argument('--d', type=int, default=2, help='Dimension of the code space (default: 2)')
    parser.add_argument('--eta', type=float, default=0.9, help='Loss parameter (default: 0.9)')
    parser.add_argument('--rounds', type=int, default=200, help='Number of optimization rounds (default: 200)')
    parser.add_argument('--file_output', type=str, help='Filename for saving optimization results')

    args = parser.parse_args()
    get_optimal_isotypic_code(fundamental_irrep_file=args.fundamental_irrep_file, code_irrep_file=args.code_irrep_file, n=args.n, d=args.d, eta=args.eta, rounds=args.rounds, file_output=args.file_output)




if __name__ == "__main__":
    n = 15
    group = '2T'
    #code_irrep  = '2T/irrep_dim2_2T_4.json'
    
    groupModel = GroupQECModel(group, n = n)
    print(groupModel.decomposition)

    code, vectors, rank = groupModel.optimize_multiplicity_code(
    [
       "irrep_dim2_2T_5.json", 
    ],
    restarts=50,
)

    fe, ic = get_code_performance(code, eta=0.95, file='n15/n15_2T_code.pkl')
    print(fe, ic, rank)
