import os
import cvxpy as cp
import pickle
import numpy as np
import math 
from scipy.stats import unitary_group
from scipy.special import factorial, gammaln

from scipy.sparse import csr_matrix, kron as spkron, csc_matrix
from scipy.sparse import vstack as spvstack 
from scipy.special import gammaln
import argparse
import time
from concurrent.futures import ThreadPoolExecutor

import warnings

warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=".*encountered in matmul.*",
)

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=".*Defaulting to the SCIPY backend for canonicalization.*",
)



# functions that gives the superoperator corresponding to a choi matrix with given input/output dimensions
# Remember that superoperator has indices jj', ii'
# Choi matrix has indices ij, i'j'
def superoperator_to_choi(superoperator, d_in, d_out):
    new_op = superoperator.reshape(d_out, d_out, d_in, d_in)
    choi = np.transpose(new_op, (2,0, 3,1)).reshape(d_in*d_out, d_in*d_out)
    return choi

# Inverse function
def choi_to_superoperator(choi, d_in, d_out):
    new_op = choi.reshape(d_in, d_out, d_in, d_out)
    superoperator = np.transpose(new_op,(1,3,0,2)).reshape(d_out*d_out, d_in*d_in)
    return superoperator


# Take partial trace wrt one of the system
# This assumes a choi matrix representation of the form X[ij, i'j']
def partial_trace(operator, dimA, dimB, to_trace=1):
    # if you want to trace out system A enter to_trace=0
    axis1 = to_trace + 0
    axis2 = to_trace+ 2

    reshaped_op = operator.reshape(dimA, dimB, dimA, dimB)

    return np.trace(reshaped_op, axis1 = axis1, axis2=axis2)

    

class RandomConstantNCode:
    def __init__(self, n= 5, d = 2, U= None):
        # The number of photons, and dimension of logical qubit
        self.n = n
        self.d = d

        # We label all possible states n_a + n_b = n with n_a
        self.target_space_basis = [np.eye(n+1)[i].reshape(n+1,1) for i in range(n+1)]

        # Generate the logical basis states
        self.logical_space_basis = [np.eye(d)[i].reshape(d,1) for i in range(d)]

        # If no unitary is provided, generate a random one
        if U is not None:
            self.U = U
        else:
            self.generate_isometry()

        # Get the superoperator matrix
        self.generate_superoperator()

        


    def generate_isometry(self):
        # generate an isometry mapping from C^d to C^{n+1}
        U = unitary_group.rvs(self.n+1)
        self.U = U[:, :self.d]

    def encode(self, density_matrix):
        return self.U@ density_matrix @ self.U.conj().T


    

    def generate_superoperator(self):
        # Define a superoperator representing the encoding map
        # We map (T_A)_{jj', ii'} = <j|A(|i><i'|)|j'>

        self.superoperator = np.kron(self.U, self.U.conj())  # shape: ((n+1)^2, d^2)
        self.choi = superoperator_to_choi(self.superoperator, self.d, self.n+1)



    def check_isometry(self):
        return self.U.conj().T@self.U






# Define the action of a_1^(l1) a_2^(l2) on |x,y>
def a(l1, l2, x,y):
    if l1 > x or l2 > y:
        return 0,x, y
    else:
        x_new = x-l1
        y_new = y-l2

        coefficient = np.sqrt((factorial(x)*factorial(y))/(factorial(x_new)*factorial(y_new)))

        return coefficient, x_new, y_new



# Find coefficient of \eta^{n_1/2} \eta^{n_2/2} acting on |x,y>
def eta_power(eta, x,y):
    return eta**(x/2 + y/2)

    

def von_neumann_entropy(rho, base=2):
    """
    Compute the von Neumann entropy S(rho) = -Tr(rho log rho).
    rho: density matrix (numpy array).
    base: logarithm base (default = 2 for bits).
    """
    # Diagonalize
    eigvals = np.linalg.eigvalsh(rho)  # safer for Hermitian matrices
    
    # Keep only nonzero eigenvalues (avoid log(0))
    eigvals = eigvals[eigvals > 1e-15]
    
    # Compute entropy
    entropy = -np.sum(eigvals * np.log(eigvals)) / np.log(base)
    return np.real_if_close(entropy)


# Note that we use C-ordering everywhere in this file
# Including the superoperator generation
def vectortize(matrix):
    return np.asarray(matrix).reshape(-1, 1, order='C')


def devectortize(vector, rows, cols=None):
    if cols is None:
        cols = rows
    return np.asarray(vector).reshape(rows, cols, order='C')


def evolve_under_superoperator(superoperator, rho, dim_out):
    evolved_vectorized_state = superoperator @ vectortize(rho)
    return devectortize(evolved_vectorized_state, dim_out)



# We represent the relevant photon loss operators for a system of two cavities
# with total photon number n and loss eta

class TruncatedLoss:
    def __init__(self, n, eta):
        self.n = n
        self.eta = eta


        # We have a total of (n+1)(n+2)/2 states with <= n photons
        self.num_states =(n+1)*(n+2)//2
        self.dim_in = n+1
        self.dim_out = self.num_states


        self.fac = factorial(np.arange(n+2))

        # Precompute triangular index map idx[x,y] -> k, with -1 for invalid
        self.idx = -np.ones((n+1, n+1), dtype=int)
        k = 0
        for x in range(n+1):
            for y in range(n+1-x):
                self.idx[x, y] = k
                k += 1

        # Find the loss channel operators in an efficient matrix representation
        self.loss_operator_matrices()
        self.generate_superoperator()
        self.generate_inclusion_isometries()
        self.generate_normalized_subchannels()


    def generate_inclusion_isometries(self):
        # Defines the canonical inclusions V_k: H_k --> S 
        # Note here k denotes the total number of photons, not the number of photons lost
        # Map the vector |x, k-x > to self.idx[x, k-x]

        # We also return the superoperator corresponding to V_k^\dagger
        # Remember that with row-first indexing, that is K \otimes K^* for K = V^\dagger

        self.inclusion_isometries = {}
        self.compression_superoperators = {}
        for k in range(self.n+1):
            rows = [self.idx[x, k-x] for x in range(k+1)]
            columns = [x for x in range(k+1)]

            isometry = csc_matrix((np.ones(k+1), (rows, columns)), shape=(self.dim_out, k+1),dtype=np.complex128 )
            self.inclusion_isometries[k] = isometry
            self.compression_superoperators[k] = spkron(isometry.T.conj(), isometry.T, format='csr')

        return self.inclusion_isometries, self.compression_superoperators
    



    # Find action of full loss operator on two modes, with l1, l2 losses for each mode on state |x,y>

    def a_fast(self, l1, l2, x,y):
        if l1 > x or l2 > y:
            return 0,x, y
        else:
            x_new = x-l1
            y_new = y-l2

            coefficient = np.sqrt((self.fac[x]*self.fac[y])/(self.fac[x_new]*self.fac[y_new]))

            return coefficient, x_new, y_new

    def loss_operator_action(self, l1, l2, x,y, backaction=True):
        coeff1, x_new, y_new = self.a_fast(l1, l2, x,y)
        coeff2 = eta_power(self.eta, x_new, y_new)

        coeff3 = np.sqrt(((1-self.eta)**l2/self.fac[l2])*((1-self.eta)**l1/self.fac[l1]))

        if backaction:
            return coeff1*coeff2*coeff3, x_new, y_new
        else:
            return coeff1, x_new, y_new
    
    
    def probability_of_k_losses(self, k):
        # Probability of losing k photons in total
        binomial_coeff = self.fac[self.n] / (self.fac[k] * self.fac[self.n - k])
        prob = binomial_coeff * (self.eta ** (self.n - k)) * ((1 - self.eta) ** k)

        return prob


    # Find the matrix elementrs of a given loss operator characterized by l1, l2 losses
    def loss_operator_matrix(self, l1, l2, backaction=True):
        rows = []
        cols = []
        data = []
        
        for x in range(self.n+1):
            y = self.n - x
            coeff, x_new, y_new = self.loss_operator_action(l1, l2, x, y, backaction=backaction)

            if coeff!=0:
                row = self.idx[x_new, y_new]
                rows.append(row)
                cols.append(x)

                data.append(coeff)


        if not data:
            return csc_matrix((self.dim_out, self.dim_in), dtype=np.complex128)
        return csc_matrix((data, (rows, cols)), shape=(self.dim_out, self.dim_in), dtype=np.complex128)



    # Find all loss operators 
    def loss_operator_matrices(self):
        self.loss_operators = []
        self.loss_operators_indexed_no_backaction = {}

        for l1 in range(self.n+1):
            for l2 in range(self.n+1-l1):
                op = self.loss_operator_matrix(l1, l2)
                op_no_backaction = self.loss_operator_matrix(l1, l2, backaction=False)
                if op.nnz:
                    self.loss_operators.append(op)

                if op_no_backaction.nnz:
                    self.loss_operators_indexed_no_backaction[(l1,l2)] = op_no_backaction

   
                    
    def generate_superoperator(self):
        # Define a superoperator representing the encoding map
        # We map (T_A)_{jj', ii'} = <j|A(|i><i'|)|j'>
        # Build the list of Kronecker products first and sum once to avoid
        # many incremental sparse matrix additions which are expensive.
        kron_list = [spkron(K, K.conj(), format='csc') for K in self.loss_operators if K.nnz]
        if kron_list:
            # summing a list of sparse matrices is reasonably efficient
            self.superoperator = sum(kron_list)
        else:
            self.superoperator = csc_matrix((self.dim_out**2, self.dim_in**2), dtype=np.complex128)



    def generate_normalized_subchannels(self):
        # generate the \overline{N}_q subchannels from the main.tex
        # These contain no noise dependence
        self.normalized_subchannel_kraus = {}
        self.normalized_complementary_subchannel_superoperators = {}
        self.normalized_restricted_subchannel_superoperators = {}

        for r in range(self.n+1):
            q = self.n - r
            kraus_q = []
            binomial_coefficient = (
                self.fac[self.n] / (self.fac[q] * self.fac[self.n-q])
            )

            for l1 in range(r+1):
                l2 = r - l1
                no_backaction_operator = self.loss_operators_indexed_no_backaction[(l1, l2)]
                normalization = np.sqrt(
                    binomial_coefficient * self.fac[l1] * self.fac[l2]
                )

                reduced_kraus = (
                    self.inclusion_isometries[q].conj().T @ no_backaction_operator
                ) / normalization
                kraus_q.append(reduced_kraus.tocsr())

            self.normalized_subchannel_kraus[q] = kraus_q
            self.normalized_restricted_subchannel_superoperators[q] = sum(
                spkron(K, K.conj(), format='csc') for K in kraus_q
            )

            complementary_kraus_q = [spvstack([K.getrow(a) for K in kraus_q], format='csr') for a in range(q+1)]
            self.normalized_complementary_subchannel_superoperators[q] = sum(spkron(M, M.conj(), format='csr') for M in complementary_kraus_q)




    
    def environment_normalized_sub_evolution(self, rho, q):
        # Evolve rho under the complement of the q-th sub-channel
        # Dim out here is r + 1 = n-q + 1 
        return evolve_under_superoperator(self.normalized_complementary_subchannel_superoperators[q], rho, self.n - q + 1 )
    

    def system_normalized_sub_evolution(self, rho, q):
        return evolve_under_superoperator(self.normalized_restricted_subchannel_superoperators[q], rho, q+1)
    

    def normalized_subsystem_coherent_info(self, rho, q):
        evolved_sys = self.system_normalized_sub_evolution(rho, q)
        evolved_env = self.environment_normalized_sub_evolution(rho, q)

        entropy_sys = von_neumann_entropy(evolved_sys)
        entropy_env = von_neumann_entropy(evolved_env)

        return entropy_sys - entropy_env
    
    def coherent_information(self, rho):
        return sum([self.normalized_subsystem_coherent_info(rho, q)*self.probability_of_k_losses(self.n-q) for q in range(self.n + 1)])
    
    


    def sanity_check(self, tol=1e-10):
        acc = np.zeros((self.dim_in, self.dim_in), dtype=np.complex128)
        for K in self.loss_operators:
            acc += K.conj().T @ K
        err = np.linalg.norm(acc - np.eye(self.dim_in))
        return err <= tol, err

        



def superoperator_to_choi_cvx(T, d_in, d_out):
    """
    CVXPY version of T[j j', i i'] -> X[i j, i' j'] (C-order).
    """
    T4 = cp.reshape(T, (d_out, d_out, d_in, d_in), order='C')            # (j, j', i, i')
    X  = cp.reshape(cp.transpose(T4, (2, 0, 3, 1)),
                    (d_in*d_out, d_in*d_out), order='C')
    return X


def choi_to_superoperator_cvx(choi, d_in, d_out):
    new_op = cp.reshape(choi,(d_in, d_out, d_in, d_out), order='C')
    superoperator = cp.reshape(cp.transpose(new_op, (1,3,0,2)), (d_out*d_out, d_in*d_in), order='C')
    return superoperator


def ptrace_out_cvx(X, d_in, d_out):
    """
    Tr_out over the output system (size d_out).
    X is a CVXPY (d_in*d_out, d_in*d_out) expression.
    Returns a (d_in, d_in) CVXPY expression.
    """
    X4 = cp.reshape(X, (d_in, d_out, d_in, d_out), order='C')  # (i, j, i', j')
    # Sum the diagonal blocks over the output index j
    return sum(X4[:, j, :, j] for j in range(d_out))




class SDPCodeOptimizer:
    def __init__(self, loss_op, T_E_init, n, d, eps=1e-5, save_progress= False, file_output = None, encoding_subspace=None, save_every=50, sector_decoding=True, decoder_workers=1, max_solver_retries=5):
        self.T_N = loss_op.superoperator
        self.loss_op = loss_op
        self.loss_sector_superoperators = loss_op.normalized_restricted_subchannel_superoperators

        self.n =n 
        self.dim_code= n + 1 # dimension of fixed photon number space                         
        self.dim_logical= d # logical dimension                         
        self.dim_physical = (n + 1) * (n + 2) // 2 # states with \leq n photons

        self.T_E = T_E_init.astype(complex)         # current encoder
        self.T_D = None
        self.M_eff = None
        self.save_every = save_every
        self.decoder_workers = decoder_workers
        self.max_solver_retries = max_solver_retries

        # add a constraint that the encoding superoperator maps into a given subspace if provided. This can be used to speed up optimization by reducing the search space, but is optional.
        # The columns of encoding_subspace should be orthogonal and span the subspace
        self.encoding_subspace = encoding_subspace
        subspace_shape_x, subspace_shape_y = encoding_subspace.shape if encoding_subspace is not None else (None, None)
        assert subspace_shape_x is None or subspace_shape_x == self.dim_code, "Encoding subspace must have dimension matching the constant-excitation space dimension."
        assert subspace_shape_y is None or subspace_shape_y >= self.dim_logical, "Encoding subspace must be at least as large as logical space."
        self.subspace_dim = subspace_shape_y

        # make sure columns are orthonormal
        subspace_row = encoding_subspace.T if encoding_subspace is not None else None
        if encoding_subspace is not None:
            for i in range(self.subspace_dim):
                for j in range(i, self.subspace_dim):
                    inner_prod = np.vdot(subspace_row[i], subspace_row[j])
                    if i == j:
                        assert np.isclose(inner_prod, 1), f"Encoding subspace vectors must be normalized. Inner product of vector {i} with itself is {inner_prod}."
                    else:
                        assert np.isclose(inner_prod, 0), f"Encoding subspace vectors must be orthogonal. Inner product of vector {i} with vector {j} is {inner_prod}."


        if encoding_subspace is not None:
            self.subspace_inclusion_superoperator = spkron(
                csc_matrix(encoding_subspace),
                csc_matrix(encoding_subspace.conj()),
                format='csc',
            )
        else:
            self.subspace_inclusion_superoperator = None
        
        # Solver tuning options
        self.solver = cp.SCS
        self.target_gap = 1e-5
        self.target_residual = 1e-4
        self.target_physical = 1e-6
        self.strict_solver_mode = False
        solver_info = {'eps': eps, 'max_iters': 10000, 'primal_residual': np.inf,
                       'dual_residual': np.inf, 'gap': np.inf, 'iterations': 0,
                       'R': np.inf, 'cap_count': 0, 'new_result': False}

        # Build the CVX problems
        self.sector_decoding = sector_decoding
        if sector_decoding:
            self.decoder_solver_info = {q: solver_info.copy() for q in range(self.dim_code)}
            for q, info in self.decoder_solver_info.items():
                info['weight'] = self.loss_op.probability_of_k_losses(self.n - q)
            self.create_sector_decoder_problems()
        else:
            self.decoder_solver_info = {'global': solver_info.copy()}
            self.create_decoder_problem()

        self.encoder_solver_info = solver_info.copy()
        self.create_encoder_problem()
        self.history = []

        self.save_progress = save_progress    
        if file_output is None:
            self.file_output = f"n{n}_eta{float(self.loss_op.eta):.6f}_d{d}_PROGRESS_s{np.random.randint(low=0, high=10000)}.pkl"
        else:
            self.file_output = file_output


        



    def create_decoder_problem(self):
        # create a variable for the choi matrix/superoperator of the decoder
        self.X_D_var = cp.Variable((self.dim_logical*self.dim_physical, self.dim_physical*self.dim_logical), complex=True, hermitian=True)
        self.T_D_var = choi_to_superoperator_cvx(self.X_D_var, self.dim_physical, self.dim_logical)

        # impose positive semidefiniteness and TP
        constraints = [self.X_D_var>>0, ptrace_out_cvx(self.X_D_var, self.dim_physical, self.dim_logical)==np.eye(self.dim_physical, dtype=complex)]

        # Use Parameter to avoid rebuilding problem each iteration
        self.T_E_param = cp.Parameter((self.dim_code**2, self.dim_logical**2), complex=True)
        self.T_E_param.value = self.T_E
        M_expr = self.T_N @ self.T_E_param

        # proportional to entanglement fidelity
        objective = cp.Maximize(cp.real(cp.trace(self.T_D_var @ M_expr)))

        self.dec_prob = cp.Problem(objective, constraints)



    def create_sector_decoder_problems(self):
        self.decoder_sector_problems = {}
        self.sector_X_D_vars = {}
        self.sector_T_D_vars = {}

        # Use Parameter to avoid rebuilding problem each iteration
        self.T_E_param = cp.Parameter((self.dim_code**2, self.dim_logical**2), complex=True)
        self.T_E_param.value = self.T_E

        for q in range(self.dim_code):
            dim_phys = q + 1
            new_var = cp.Variable((self.dim_logical*dim_phys, self.dim_logical*dim_phys), complex=True, hermitian=True)
            new_var_T_D = choi_to_superoperator_cvx(new_var, dim_phys, self.dim_logical)
            self.sector_X_D_vars[q] = new_var 
            self.sector_T_D_vars[q] = new_var_T_D

            constraints = [new_var>>0, ptrace_out_cvx(new_var, dim_phys, self.dim_logical)==np.eye(dim_phys, dtype=complex)]
            M_expr = self.loss_sector_superoperators[q] @ self.T_E_param

            objective = cp.Maximize(cp.real(cp.trace(new_var_T_D@M_expr)))
            self.decoder_sector_problems[q] = cp.Problem(objective, constraints)


    def get_full_decoder(self):
        if self.T_D is not None:
            return self.T_D

        if not self.sector_decoding:
            raise RuntimeError("The full decoder has not been solved yet.")

        if any(
            self.sector_T_D_vars[q].value is None
            for q in range(self.dim_code)
        ):
            raise RuntimeError("The sector decoders must be solved first.")

        self.T_D = sum(
            csr_matrix(self.sector_T_D_vars[q].value)
            @ self.loss_op.compression_superoperators[q]
            for q in range(self.dim_code)
        )

        return self.T_D


    def create_encoder_problem(self):
        # If we're given an encoder subspace, we can add a constraint that the encoder superoperator maps into that subspace. This is optional and can be used to speed up optimization by reducing the search space.
        if self.encoding_subspace is not None:
            encoding_output_dim = self.subspace_dim
            self.X_E_var = cp.Variable((encoding_output_dim*self.dim_logical, encoding_output_dim*self.dim_logical), complex=True, hermitian=True)
            self.T_E_var = self.subspace_inclusion_superoperator @ choi_to_superoperator_cvx(self.X_E_var, self.dim_logical, encoding_output_dim)


        else:# create a variable for the choi/superoperator matrix of encoder
            encoding_output_dim = self.dim_code
            self.X_E_var = cp.Variable((encoding_output_dim*self.dim_logical, encoding_output_dim*self.dim_logical), complex=True, hermitian=True)
            self.T_E_var = choi_to_superoperator_cvx(self.X_E_var, self.dim_logical, encoding_output_dim)

        # Enforce positive definite and TP
        constraints = [self.X_E_var>>0, ptrace_out_cvx(self.X_E_var, self.dim_logical, encoding_output_dim)==np.eye(self.dim_logical,dtype=complex)]

        # Use parameter to speed up. Delay assigning value until available.
        self.M_eff_param = cp.Parameter((self.dim_logical**2, self.dim_code**2), complex=True)
        if self.M_eff is not None:
            self.M_eff_param.value = self.M_eff

        objective = cp.Maximize(cp.real(cp.trace(self.M_eff_param @ self.T_E_var)))
        self.enc_prob = cp.Problem(objective, constraints)


    def encoded_coherent_information(self):
        # Compute the coherent information of the full loss channel for the current encoded state
        # from the current encoder superoperator T_E

        input_state = np.eye(self.dim_logical, dtype=complex) / self.dim_logical
        rho_encoded = evolve_under_superoperator(
            self.T_E, input_state, self.dim_code
        )

        I_c = self.loss_op.coherent_information(rho_encoded)
        return I_c

    

    def _solve_sector_decoder(self, q):
        info = self.decoder_solver_info[q]
        problem = self.decoder_sector_problems[q]
        dim_phys = q + 1
        fallback = np.kron(
            np.eye(dim_phys),
            np.eye(self.dim_logical) / self.dim_logical,
        )
        return self._solve_validated(
            problem, self.sector_X_D_vars[q], dim_phys,
            self.dim_logical, info, fallback=fallback,
        )

    def solve_decoder_problem(self):
        if self.sector_decoding:
            sectors = range(self.dim_code - 1, -1, -1)
            if self.decoder_workers == 1:
                accepted = [self._solve_sector_decoder(q) for q in sectors]
            else:
                with ThreadPoolExecutor(max_workers=self.decoder_workers) as pool:
                    accepted = list(pool.map(self._solve_sector_decoder, sectors))

            self.T_D = None
            self.M_eff = sum(self.loss_op.probability_of_k_losses(self.n - q)*self.sector_T_D_vars[q].value@self.loss_sector_superoperators[q] for q in range(self.dim_code))
            self.M_eff_param.value = self.M_eff
            return all(accepted)

        else:
            info = self.decoder_solver_info['global']
            fallback = np.kron(
                np.eye(self.dim_physical),
                np.eye(self.dim_logical) / self.dim_logical,
            )
            accepted = self._solve_validated(
                self.dec_prob, self.X_D_var, self.dim_physical,
                self.dim_logical, info, fallback=fallback,
            )
            self.T_D = self.T_D_var.value
            self.M_eff = self.T_D_var.value@self.T_N
            self.M_eff_param.value = self.M_eff
            return accepted

    

    def solve_encoder_problem(self):
        info = self.encoder_solver_info
        dim_output = self.subspace_dim if self.encoding_subspace is not None else self.dim_code
        accepted = self._solve_validated(
            self.enc_prob, self.X_E_var, self.dim_logical,
            dim_output, info,
        )
        if accepted:
            self.T_E_param.value = self.T_E_var.value
            self.T_E = self.T_E_var.value
        return accepted

    @staticmethod
    def _record_solver_info(problem, info):
        stats = problem.solver_stats
        scs_info = (stats.extra_stats or {}).get('info', {})
        info.update(primal_residual=abs(float(scs_info.get('res_pri', np.inf))),
                    dual_residual=abs(float(scs_info.get('res_dual', np.inf))),
                    gap=abs(float(scs_info.get('gap', np.inf))),
                    iterations=int(scs_info.get('iter', stats.num_iters or 0)), new_result=True)

    @staticmethod
    def _choi_errors(X, d_in, d_out):
        if X is None or not np.all(np.isfinite(X)):
            return np.inf, np.inf

        X = (X + X.conj().T) / 2
        try:
            cp_error = np.maximum(-np.linalg.eigvalsh(X), 0).sum() / d_in
        except np.linalg.LinAlgError:
            return np.inf, np.inf
        tp_error = np.linalg.norm(
            partial_trace(X, d_in, d_out) - np.eye(d_in), 'fro'
        ) / np.sqrt(d_in)
        return cp_error, tp_error

    def _accuracy_ratio(self, info):
        weight = info.get('weight', 1.0)
        R = max(
            weight * info['gap'] / self.dim_logical**2 / self.target_gap,
            info['primal_residual'] / self.target_residual,
            weight * info['dual_residual'] / self.target_residual,
        )
        return R if np.isfinite(R) else np.inf

    def _solve_validated(self, problem, X, d_in, d_out, info, fallback=None):
        old_value = None if X.value is None else X.value.copy()
        ratio_limit = 2 if self.strict_solver_mode else 10
        physical_limit = max(self.target_physical, 2e-9 * d_out)

        for attempt in range(self.max_solver_retries + 1):
            try:
                problem.solve(
                    solver=self.solver, warm_start=True, eps=info['eps'],
                    max_iters=info['max_iters'], verbose=False,
                )
            except cp.error.SolverError:
                severity = np.inf
            else:
                self._record_solver_info(problem, info)
                cp_error, tp_error = self._choi_errors(X.value, d_in, d_out)
                R = self._accuracy_ratio(info)
                info.update(cp_error=cp_error, tp_error=tp_error, R=R)

                if (problem.status in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE)
                        and cp_error <= physical_limit
                        and tp_error <= physical_limit
                        and R <= ratio_limit):
                    return True

                severity = max(
                    cp_error / physical_limit,
                    tp_error / physical_limit,
                    R / ratio_limit,
                )

            if attempt == self.max_solver_retries:
                break

            if np.isfinite(severity):
                decades = int(np.clip(
                    np.ceil(np.log10(max(severity, 1.01))), 1, 3
                ))
            else:
                decades = 1
            info['eps'] = max(1e-9, info['eps'] / 10**decades)
            info['max_iters'] = min(
                100000, info['max_iters'] * 2**min(decades, 2)
            )

        X.value = old_value if old_value is not None else fallback
        warnings.warn(
            "SDP remained insufficiently accurate after retries; "
            "keeping the previous valid iterate."
        )
        return False

    def update_solver_settings(self):
        infidelity = max(0.0, 1.0 - min(self.Fe, 1.0))
        self.target_gap = min(
            self.target_gap, np.clip(0.01 * infidelity, 1e-7, 1e-5)
        )
        self.target_residual = 10 * self.target_gap
        self.target_physical = min(
            self.target_physical, np.clip(0.01 * infidelity, 1e-7, 1e-6)
        )
        for info in [self.encoder_solver_info, *self.decoder_solver_info.values()]:
            if not info['new_result']:
                continue
            info['new_result'] = False
            R = self._accuracy_ratio(info)
            info['R'] = R
            if R > 2:
                info['eps'] = max(
                    1e-9,
                    min(info['eps'], self.target_gap)
                    / max(1.0, np.sqrt(R / 2)),
                )
            else:
                info['eps'] = min(info['eps'], self.target_gap)

            info['cap_count'] = info['cap_count'] + 1 if info['iterations'] >= info['max_iters'] and R > 1 else 0
            if info['cap_count'] >= 2:
                info['max_iters'] = min(2*info['max_iters'], 100000)
                info['cap_count'] = 0

    def calculate_performance(self,step=0):
        # Calculate entanglement fidelity of full channel
        with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
            full_channel = self.M_eff@ self.T_E
            
        self.Fe = np.real(np.trace(full_channel)) / (self.dim_logical**2)
        self.update_solver_settings()

        # Calculate coherent information of full channel
        if step % self.save_every == 0:
            self.I_c = self.encoded_coherent_information()

    def save(self):
        decoder = self.get_full_decoder()
        n, eta, d, rounds, seed, dec, enc, fe, I_c = self.dim_code-1, self.loss_op.eta, self.dim_logical, len(self.history), None, decoder, self.T_E, self.Fe, self.I_c
        results = (n, eta, d, rounds, seed, dec, enc, fe, I_c)


        with open(self.file_output, "wb") as f:
            pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)
        
        

    def step(self, step_num):
        # Optimize decoder and pass value as parameter to the encoding problem
        # Solve decoder with configured solver options
        
        if not self.solve_decoder_problem():
            warnings.warn("Skipping this iteration because the decoder did not pass validation.")
            return False
        # Solve encoder and pass value as parameter to decoding probelm
        # Ensure parameter is set before solving encoder
        if not self.solve_encoder_problem():
            warnings.warn("Skipping this iteration because the encoder did not pass validation.")
            return False

        self.calculate_performance(step=step_num)

        self.history.append((self.Fe, self.I_c))

        if self.save_progress and step_num % self.save_every == 0:
            self.save()

        return True

    def _final_polish(self):
        self.strict_solver_mode = True
        try:
            self.target_gap = 1e-8
            self.target_residual = 1e-7
            self.target_physical = 1e-8
            for info in [self.encoder_solver_info, *self.decoder_solver_info.values()]:
                info['eps'] = min(info['eps'], 1e-8)
                info['max_iters'] = max(info['max_iters'], 50000)

            previous_Fe = getattr(self, 'Fe', None)
            converged = False
            strict_update_accepted = False
            for _ in range(3):
                if not self.solve_decoder_problem():
                    continue
                if not self.solve_encoder_problem():
                    continue
                strict_update_accepted = True
                self.calculate_performance(step=1)

                tolerance = max(2e-8, 1e-4 * max(1.0 - self.Fe, 0.0))
                if previous_Fe is not None and abs(self.Fe - previous_Fe) <= tolerance:
                    converged = True
                    break
                previous_Fe = self.Fe

            decoder_accepted = self.solve_decoder_problem()
            self.get_full_decoder()
            self.calculate_performance(step=0)
            if not decoder_accepted:
                warnings.warn("The final decoder retained a previous valid iterate.")
            if not strict_update_accepted:
                warnings.warn(
                    "Strict final polishing did not complete; the output "
                    "retains an earlier valid encoder."
                )
            elif not converged:
                warnings.warn(
                    "The strict solves passed, but the alternating "
                    "optimization may benefit from more rounds."
                )
        finally:
            self.strict_solver_mode = False

        

    def run(self, rounds=50):
        for i in range(rounds):
            self.step(i)

        self._final_polish()

    def run_timed(self, rounds=50, verbose=False):
        import time
        times = []
        for i in range(rounds):
            t0 = time.time()
            self.step(i)
            times.append(time.time()-t0)

        self._final_polish()

        if verbose:
            print(f'Average time per round: {sum(times)/len(times):.3f}s over {len(times)} rounds')
        return 
    

    
        

def optimize(n, eta, d=2, rounds=50, save_progress=False, outpath=None, input_encoder=None, sector_decoding=True, decoder_workers=1):
    # Start with a random code
    if input_encoder is  None:
        myCode = RandomConstantNCode(n=n, d=d)
        T_E_init = myCode.superoperator
    else:
        T_E_init = input_encoder
    
    # Get a loss superoperato
    herLoss = TruncatedLoss(n, eta)

    # Now optimize
    myOptimizer = SDPCodeOptimizer(herLoss, T_E_init, n, d, save_progress=save_progress, file_output=outpath, sector_decoding=sector_decoding, decoder_workers=decoder_workers)
    # Use timed run for quick profiling
    times = myOptimizer.run_timed(rounds=rounds, verbose=False)
    return myOptimizer.T_D, myOptimizer.T_E, myOptimizer.Fe, myOptimizer.I_c






def load_optimization_results(pickle_path, as_dense=False):
    """Load results saved by parallel_optimize (single pickle file).

    Returns a list of tuples (eta, decoder, encoder, fidelity, I_c).
    If as_dense=True, attempts to convert sparse matrices (scipy) to dense numpy arrays.
    """
    
    if not os.path.exists(pickle_path):
        raise FileNotFoundError(pickle_path)
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f)

    if as_dense:
        try:
            from scipy.sparse import issparse
        except Exception:
            issparse = lambda x: False

        new_data = []
        for item in data:
            eta, dec, enc, fe, Ic = item
            if issparse(dec):
                dec = dec.toarray()
            if issparse(enc):
                enc = enc.toarray()
            new_data.append((eta, dec, enc, fe, Ic))
        return new_data

    return data




def run_single_optimize(n, eta, d=2, rounds=50, outpath=None, seed=None, input_file = None, folder=None, decoder_workers=1):
    """Run a single optimize() call and save the results to a pickle.

    Returns the tuple saved: (n, eta, d, rounds, seed, decoder, encoder, fidelity, I_c)
    """
    if seed is not None and input_file is not None:
        raise ValueError(
            'seed and input_file are mutually exclusive: use seed for a '
            'random initialization or input_file for a warm start'
        )

    print(f"Running optimize with n={n}, eta={eta}, d={d}, rounds={rounds}, seed={seed}")

    if seed is not None:
        # seed numpy RNG for reproducible RandomConstantNCode initialization
        np.random.seed(seed)

    T_E_init = None
    if input_file is not None:
        # Load an existing encoding from file
        with open(input_file, "rb") as f:
            data = pickle.load(f)
            # Expecting (n, eta, d, rounds, seed, dec, enc, fe, I_c)
            if len(data) != 9:
                raise ValueError(f"Input file {input_file} does not contain expected data format.")
            n_org, _, d_org, _, _, _, enc, _, _ = data
            if n_org != n or d_org != d:
                raise ValueError(f"Input file {input_file} has incompatible n or d (got n={n_org}, d={d_org}, expected n={n}, d={d}).")
            print(f"Loaded encoder from {input_file} for initialization.")

        T_E_init = enc
    
    
    if outpath is None:
        base_tag = (
            f"n{n}_eta{float(eta):.6f}_d{d}_r{rounds}_updated"
        )
        if input_file is not None:
            input_stem = os.path.splitext(
                os.path.basename(input_file)
            )[0]
            tag = f"{base_tag}_warm_from_{input_stem}"
        else:
            tag = f"{base_tag}_s{(seed if seed is not None else 0)}"
        outpath = os.path.abspath(tag + '.pkl')

    if folder is not None:
        os.makedirs(folder, exist_ok=True)
        outpath = os.path.join(folder, os.path.basename(outpath))


    dec, enc, fe, I_c = optimize(n, eta, d=d, rounds=rounds, save_progress=True, outpath=outpath, input_encoder= T_E_init, decoder_workers=decoder_workers)

    results = (n, eta, d, rounds, seed, dec, enc, fe, I_c)
    

    with open(outpath, "wb") as f:
        pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)




def main():
    parser = argparse.ArgumentParser(description='Run optimize() for a single (n,eta,d,rounds) and save results.')
    parser.add_argument('--n', type=int, required=True, help='total photon number n')
    parser.add_argument('--eta', type=float, required=True, help='channel parameter eta')
    parser.add_argument('--d', type=int, default=2, help='logical dimension d (default: 2)')
    parser.add_argument('--rounds', type=int, default=50, help='number of optimization rounds (default: 50)')
    initialization_group = parser.add_mutually_exclusive_group()
    initialization_group.add_argument('--seed', type=int, default=None, help='RNG seed for a reproducible random initialization')
    initialization_group.add_argument('--input_file', type=str, default=None, help='existing result file to use as a warm-start encoder')
    parser.add_argument('--out', type=str, default=None, help='optional output file path for pickle')
    parser.add_argument('--folder' , type=str, default=None, help='optional output folder to save with auto-generated filename')
    parser.add_argument('--workers', type=int, default=1, help='number of sector-decoder worker threads (default: 1)')

    args = parser.parse_args()

    run_single_optimize(args.n, args.eta, d=args.d, rounds=args.rounds, outpath=args.out, seed=args.seed, input_file = args.input_file, folder=args.folder, decoder_workers=args.workers)


if __name__ == '__main__':
    # prefer CLI main for single-run usage
    main()
    
