import numpy as np
import math 
from scipy.linalg import expm
from scipy.optimize import differential_evolution, minimize
import matplotlib.pyplot as plt

from faster_gen_dual_rail import *
from group_rep_codes_classes import *
import os
import re

from qutip import Qobj, ket2dm
from qutip.wigner import spin_wigner

from sympy.physics import wigner as spw 

import pickle 


import math, numpy as np
from sympy import S, N




def encode(encoding_superoperator, density_matrix):
    dim_out = int(np.sqrt(encoding_superoperator.shape[0]))
    return evolve_under_superoperator(
        encoding_superoperator,
        density_matrix,
        dim_out,
    )

# Assuming the input is a density matrix of a pure state 
def purify_encoded_state(density_matrix):
    # Hermitian / density-matrix case (recommended)
    eigvals, eigvecs = np.linalg.eigh(density_matrix)   # eigvals sorted ascending
    idx = np.argmax(eigvals)                           # index of largest eigenvalue
    principal_vec = eigvecs[:, idx]                    # eigenvector (column)
    # normalize (should be normalized already from eigh, but safe):
    principal_vec = principal_vec / np.linalg.norm(principal_vec)

    # projector onto the principal eigenvector (pure state)
    proj = np.outer(principal_vec, principal_vec.conj())


    difference = np.linalg.norm(proj-density_matrix, ord='fro')
    if difference > 1e-6:
        print("Warning: significant difference between density matrix and purified state:", difference) 
    
    return principal_vec # return the projector onto the principal eigenvector (pure state)

def extract_amplitudes_and_phases(vector_state):
    magnitudes = np.abs(vector_state)
    phases = np.angle(vector_state)
    return magnitudes, phases


def make_orthogonal(vec_a, vec_b):
    # take two almost orthogonal vectors and return an orthogonal basis for their projector
    proj = projector(vec_a, vec_b)
    eigvals, eigvecs = np.linalg.eigh(proj)
    idx = np.argsort(eigvals)[::-1]  # descending
    basis_vec1 = eigvecs[:, idx[0]]
    basis_vec2 = eigvecs[:, idx[1]]
    return basis_vec1/ np.linalg.norm(basis_vec1), basis_vec2/ np.linalg.norm(basis_vec2)

def logical_states(encoding_superoperator):
    """
    Return the list of logical state vectors encoded by `encoding_superoperator`.
    The encoding superoperator maps d^2 -> (n+1)^2 so we infer logical dimension d
    from the input size and encode the computational basis density matrices |i><i|.

    Returns a list of length d where each element is a 1D numpy vector (encoded state).
    """
    # encoding_superoperator shape: (d_out^2, d_in^2)
    if encoding_superoperator.ndim != 2:
        raise ValueError('encoding_superoperator must be a 2D array')
    d_in_sq = encoding_superoperator.shape[1]
    d = int(round(np.sqrt(d_in_sq)))

    logicals = []
    for i in range(d):
        # basis density matrix |i><i|
        rho = np.zeros((d, d), dtype=complex)
        rho[i, i] = 1.0
        enc_dm = encode(encoding_superoperator, rho)
        vec = purify_encoded_state(enc_dm)
        logicals.append(vec)

    return logicals


def load_data(file):
    

    data = load_optimization_results(file)
    etas = []
    decs = []
    encs = []
    fes = []
    Ics = []    

    for i, item in enumerate(data):
            eta, dec, enc, fe, Ic = item
            etas.append(eta)
            decs.append(dec)
            encs.append(enc)
            fes.append(fe)
            Ics.append(Ic)

    return etas, decs, encs, fes, Ics


 


def spin_wigner_rotation_matrix(j_float, theta, phi):
    # make 2j an *integer* first; this handles half-integers exactly
    twoj = int(round(2*j_float))
    j = twoj / 2.0                  # numeric for sizes
    # Equivalent to SymPy's small-d matrix with both basis indices reversed.
    # The positive sign preserves the existing ascending-m viewing convention.
    d_asc = expm(1j * theta * SpinSector(twoj).Jy)

    # build D in your ordering (ascending m = -j..+j)
    dim = twoj + 1
    D = np.zeros((dim, dim), dtype=complex)
    m_vals = np.arange(-j, j+1)     # ascending

    for i, m in enumerate(m_vals):
        for k, mp in enumerate(m_vals):
            D[i, k] = d_asc[i, k] * np.exp(1j * phi * (mp - m))  # + sign

    return D


def spin_wigner_euler_rotation_matrix(j_float, alpha, beta, gamma):
    """
    Full z-y-z Euler spin-j rotation in the same ascending-m basis used by
    `schwinger_to_spin` and `spin_wigner_rotation_matrix`.

    The middle factor is exactly ``spin_wigner_rotation_matrix(j, beta, 0)``;
    this preserves the sign convention of the existing visual rotation. This
    is a separate helper, and the old two-angle rotation helper is unchanged.
    """
    twoj = int(round(2 * j_float))
    Ry = spin_wigner_rotation_matrix(j_float, beta, 0.0)
    Rz_left = SpinSector(twoj).Jz_rotation_matrix(alpha)
    Rz_right = SpinSector(twoj).Jz_rotation_matrix(gamma)
    return Rz_left @ Ry @ Rz_right

def extract_photon_number_from_filename(filename):
    m = re.search(r'n(\d+)', filename)
    return int(m.group(1)) if m else None




def schwinger_to_spin(state):
    # state with total photon number n has length n+1
    n_plus_1 = len(state)
    j = (n_plus_1 -1)/ 2

    # We know we map |j,m> to |n_a, n_b> with n_a = j+m, n_b = j-m
    new_state = np.zeros((n_plus_1,), dtype=complex)
    for i,m in enumerate(np.arange(-j,j+1 )):
        na = int(j+m)
        new_state[i] = state[na]

    return new_state




def commutator(A,B):
    return A @ B - B @ A

def commute(A,B, eps):
    C = commutator(A,B)
    norm = np.linalg.norm(C, ord='fro')
    return norm < eps

def rotate_state(state, theta, phi):
    n = len(state)
    j = (n - 1) / 2
    D = spin_wigner_rotation_matrix(j, theta, phi)
    return D @ state




def rotate_code(original_supoperator, theta, phi):
    logicals = logical_states(original_supoperator)
    spin_logicals = [schwinger_to_spin(l) for l in logicals]
    rotated = [rotate_state(s, theta, phi) for s in spin_logicals]
    return get_superoperator_from_logical_states(*rotated)


def rotate_code_euler(original_superoperator, alpha, beta, gamma):
    """Apply a full z-y-z Euler rotation to every logical codeword."""
    logicals = logical_states(original_superoperator)
    spin_logicals = [schwinger_to_spin(logical) for logical in logicals]
    j = (spin_logicals[0].shape[0] - 1) / 2
    rotation = spin_wigner_euler_rotation_matrix(j, alpha, beta, gamma)
    rotated = [rotation @ logical for logical in spin_logicals]
    return get_superoperator_from_logical_states(*rotated)


def reflect_code(original_superoperator):
    """Reflect a code by complex conjugation in the fixed Fock basis."""
    return original_superoperator.conj()


def rotate_projector(projector_matrix, theta, phi):
    """
    Rotate a code projector directly in the fixed-spin basis.
    """
    n = projector_matrix.shape[0]
    j = (n - 1) / 2
    D = spin_wigner_rotation_matrix(j, theta, phi)
    return D @ projector_matrix @ D.conj().T


def rotate_projector_euler(projector_matrix, alpha, beta, gamma):
    """
    Full-gauge version of `rotate_projector` using three Euler angles.

    This is a separate helper; the old two-angle projector rotation is unchanged.
    """
    n = projector_matrix.shape[0]
    j = (n - 1) / 2
    D = spin_wigner_euler_rotation_matrix(j, alpha, beta, gamma)
    return D @ projector_matrix @ D.conj().T




def compute_spin_wigner_for_rho(rho, theta, phi):
    """Return the spin Wigner function in the project's conventions.

    Project matrices use ascending magnetic number ``m=-j,...,j``, whereas
    QuTiP uses descending magnetic number ``m=j,...,-j``. QuTiP also returns
    its Wigner array with shape ``(len(phi), len(theta))``. Reverse the spin
    basis before calling QuTiP, then transpose the returned Wigner values so
    that all three returned arrays have shape ``(len(theta), len(phi))``.
    """
    rho_data = rho.full()
    if rho.isket:
        qutip_rho = Qobj(rho_data[::-1, :], dims=rho.dims)
    elif rho.isbra:
        qutip_rho = Qobj(rho_data[:, ::-1], dims=rho.dims)
    else:
        qutip_rho = Qobj(rho_data[::-1, ::-1], dims=rho.dims)

    result = spin_wigner(qutip_rho, theta, phi)
    W = np.asarray(result[0] if isinstance(result, (list, tuple)) else result)
    expected_qutip_shape = (len(phi), len(theta))
    if W.shape != expected_qutip_shape:
        raise ValueError(
            'Unexpected spin_wigner output shape '
            f'{W.shape}; expected {expected_qutip_shape}'
        )
    W = W.T
    TH, PH = np.meshgrid(theta, phi, indexing='ij')
    return TH, PH, W





def compute_wigners_for_encoding(encoding_superoperator, theta_res=300, phi_res=300, logicals=None, projector_only=False):
    """Compute the projector Wigner, optionally omitting logical-state Wigners."""
    if logicals is None:
        logicals = logical_states(encoding_superoperator)
    else:
        logicals = logicals

    spin_states = [schwinger_to_spin(s) for s in logicals]
    P_spin = sum(np.outer(s, s.conj()) for s in spin_states)
    d = len(spin_states)


    QP = Qobj(P_spin)
    Qs = [ket2dm(Qobj(s)) for s in spin_states]

    theta = np.linspace(0, np.pi, theta_res)
    phi = np.linspace(0, 2 * np.pi, phi_res)
    j = (len(spin_states[0]) - 1) / 2

    TH_proj, PH_proj, W_proj = compute_spin_wigner_for_rho(QP, theta, phi)

    results = {
        'projector': (TH_proj, PH_proj, W_proj),
        'theta': theta,
        'phi': phi,
        'j': j,
        'd': d
    }
    if projector_only:
        return results
    # add each logical state's Wigner under keys logical0, logical1, ... logical{d-1}
    for idx, Q in enumerate(Qs):
        TH_i, PH_i, W_i = compute_spin_wigner_for_rho(Q, theta, phi)
        results[f'logical{idx}'] = (TH_i, PH_i, W_i)

    return results


def _pi_fraction_label(val, denom=8):
    """Return a string label for val (radians) as a fraction of pi with denominator up to `denom`."""
    frac = val / np.pi
    nearest = int(np.round(frac * denom))
    num = nearest
    den = denom
    g = math.gcd(abs(num), den)
    if g != 0:
        num //= g
        den //= g
    if num == 0:
        return '0'
    sign = '-' if num < 0 else ''
    num = abs(num)
    if den == 1:
        return f'{sign}π' if num == 1 else f'{sign}{num}π'
    if num == 1:
        return f'{sign}π/{den}'
    return f'{sign}{num}π/{den}'




def spherical_surface_from_TH_PH(TH, PH, flip_phi=False):
    """
    Convert TH,PH (theta,phi grids) to Cartesian surface coordinates.
    flip_phi=False by default so sphere phi-orientation matches the planar heatmap.
    """
    if flip_phi:
        PHm = -PH
    else:
        PHm = PH
    X = np.sin(TH) * np.cos(PHm)
    Y = np.sin(TH) * np.sin(PHm)
    Z = np.cos(TH)
    return X, Y, Z

def plot_wigner_sphere(TH, PH, W, fname_base, title, j,d=None,
                       interactive=True, output_folder='.',
                       flip_phi=False, show_reference=False):
    """Static PNG + interactive HTML for the spherical Wigner plot.

    - Static Matplotlib sphere uses facecolors from W (no changes).
    - Interactive Plotly sphere:
        * Surface shows the colormap but does NOT capture hover events.
        * A transparent Scatter3d of slightly inflated vertices carries per-vertex
          customdata (x,y,z, φ, θ, W) so hover is stable and never occluded.
    """
    import os
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib import cm

    # ---------- Static Matplotlib sphere ----------
    X, Y, Z = spherical_surface_from_TH_PH(TH, PH, flip_phi=flip_phi)

    vmin = float(np.nanmin(W))
    vmax = float(np.nanmax(W))
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap('RdBu_r')

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.plot_surface(X, Y, Z, facecolors=cmap(norm(W)), rcount=100, ccount=100, antialiased=False, shade=False)
    dimension_label = f', K={d}' if d is not None and 'K=' not in title.replace(' ', '') else ''
    ax.set_title(title + (f', j={j}' if j is not None else '') + dimension_label)
    ax.set_axis_off()

    if show_reference:
        th0 = np.pi/2
        ph0 = 0.0
        ph_plot = -ph0 if flip_phi else ph0
        x0 = np.sin(th0) * np.cos(ph_plot)
        y0 = np.sin(th0) * np.sin(ph_plot)
        z0 = np.cos(th0)
        ax.scatter([x0], [y0], [z0], color='r', s=40)

    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(W.ravel())
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.6)
    cbar.set_label('Wigner value')

    os.makedirs(output_folder, exist_ok=True)
    fname_png = os.path.join(output_folder, fname_base + '_wigner_sphere.png')
    fig.savefig(fname_png, dpi=150)
    plt.close(fig)

    # ---------- Interactive Plotly sphere ----------
    html_name = None
    if interactive:
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots

            # Match plane orientation (φ along x, θ along y) by transposing arrays
            THs = TH                               # (nφ, nθ)
            PHs = (-PH) if flip_phi else PH     # (nφ, nθ)
            Zc  = W                                # (nφ, nθ) surface colors

            # Geometry
            Xs = np.sin(THs) * np.cos(PHs)
            Ys = np.sin(THs) * np.sin(PHs)
            Zs = np.cos(THs)

            vmin_i = float(np.nanmin(Zc))
            vmax_i = float(np.nanmax(Zc))

            figp = make_subplots(rows=1, cols=1, specs=[[{'type': 'surface'}]])

            # 1) Surface for the textured sphere (no hover so it won't steal events)
            surfly = go.Surface(
                x=Xs, y=Ys, z=Zs,
                surfacecolor=Zc,
                colorscale='RdBu',
                cmin=vmin_i, cmax=vmax_i,
                colorbar=dict(title='Wigner value'),
                showscale=True,
                hoverinfo='skip'  # <-- critical: do not capture hover on the surface
            )
            figp.add_trace(surfly, row=1, col=1)

            # 2) Vertex layer (slightly inflated) to ensure hover is never occluded
            eps = 5e-2  # tiny radial inflation; visually invisible but fixes occlusion
            Xi = ((1.0 + eps) * Xs).ravel()
            Yi = ((1.0 + eps) * Ys).ravel()
            Zi = ((1.0 + eps) * Zs).ravel()

            # Per-vertex angle + value (aligned with geometry indices)
            Ph_i = np.mod(PHs.ravel(), 2*np.pi)
            Th_i = THs.ravel()
            W_i  = Zc.ravel()

            # customdata: x,y,z, φ, θ, W  (vertex-aligned)
            cd_i = np.column_stack([Xi, Yi, Zi, Ph_i, Th_i, W_i])

            figp.add_trace(go.Scatter3d(
                x=Xi, y=Yi, z=Zi,
                mode='markers',
                marker=dict(
                    size=8,
                    color='rgba(0,0,0,0.0001)',  # nearly invisible but hoverable
                    line=dict(width=0)
                ),
                customdata=cd_i,
                hovertemplate=(
                    'φ(v): %{customdata[3]:.4f}<br>'
                    'θ(v): %{customdata[4]:.4f}<br>'
                    'W(v): %{customdata[5]:.4f}<extra></extra>'
                ),
                name='vertices',
                showlegend=False
            ))

            # True sphere aspect; axes off
            figp.update_scenes(
                xaxis_visible=False, yaxis_visible=False, zaxis_visible=False,
                aspectmode='data'
            )
            figp.update_layout(
                title_text=title + (f' (j={j})' if j is not None else ''),
                width=700, height=700
            )

            html_name = os.path.join(output_folder, fname_base + '_wigner_sphere_interactive.html')
            figp.write_html(html_name, include_plotlyjs='cdn')
            try:
                figp.show()
            except Exception:
                pass

        except Exception as e:
            print('Plotly not installed or failed; skipping interactive sphere.', e)

    return fname_png, html_name

def plot_wigner_plane(
    TH, PH, W, fname_base, title, j,d=None,
    interactive=True, output_folder='.',
    pi_step=np.pi/4, show_reference=False,
    # --- NEW: static (matplotlib) font-size knobs ---
    fig_size=(10, 4.5),
    title_fs=24,
    axis_label_fs=20,
    tick_fs=14,
    cbar_label_fs=16,
    cbar_tick_fs=14,
):
    """Save planar heatmap PNG and optional interactive HTML for the given Wigner (TH,PH,W)."""
    theta = TH[:, 0]
    phi = PH[0, :]

    vmin = float(np.nanmin(W))
    vmax = float(np.nanmax(W))

    # --- only change: allow tuning figure size + fonts ---
    fig, ax = plt.subplots(figsize=fig_size)
    im = ax.imshow(
        W, origin='lower', aspect='auto', cmap='RdBu_r',
        extent=[phi[0], phi[-1], theta[0], theta[-1]],
        vmin=vmin, vmax=vmax
    )

    # ticks spaced by pi_step
    phi_ticks = np.arange(0, 2 * np.pi + 1e-12, pi_step)
    theta_ticks = np.arange(0, np.pi + 1e-12, pi_step)
    phi_labels = [_pi_fraction_label(v, denom=8) for v in phi_ticks]
    theta_labels = [_pi_fraction_label(v, denom=8) for v in theta_ticks]

    ax.set_xticks(phi_ticks)
    ax.set_xticklabels(phi_labels, rotation=0, fontsize=tick_fs)
    ax.set_yticks(theta_ticks)
    ax.set_yticklabels(theta_labels, fontsize=tick_fs)

    if show_reference:
        ax.plot(0.0, np.pi/2.0, 'ro', markersize=6)

    ax.set_xlabel('φ (rad)', fontsize=axis_label_fs)
    ax.set_ylabel('θ (rad)', fontsize=axis_label_fs)
    dimension_label = f', K={d}' if d is not None and 'K=' not in title.replace(' ', '') else ''
    ax.set_title(title + (f', j={j}' if j is not None else '') + dimension_label, fontsize=title_fs, pad=12)
    ax.grid(True, linestyle=':', linewidth=0.6, axis='both')

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Wigner value', fontsize=cbar_label_fs)
    cbar.ax.tick_params(labelsize=cbar_tick_fs)

    os.makedirs(output_folder, exist_ok=True)
    plane_png = os.path.join(output_folder, fname_base + '_wigner_plane.png')
    fig.tight_layout()
    fig.savefig(plane_png, dpi=150)
    plt.close(fig)

    # ---- INTERACTIVE PART: EXACTLY YOUR ORIGINAL CODE BELOW ----
    html_name_plane = None
    if interactive:
        try:
            import plotly.graph_objects as go

            # Compute Cartesian coords for each (θ, φ)
            # Shapes: (len(theta), len(phi))
            TH_grid, PH_grid = np.meshgrid(theta, phi, indexing='ij')
            X = np.sin(TH_grid) * np.cos(PH_grid)
            Y = np.sin(TH_grid) * np.sin(PH_grid)
            Z = np.cos(TH_grid)

            custom = np.stack([X, Y, Z], axis=-1)

            heat = go.Figure(data=go.Heatmap(
                z=W,
                x=phi,
                y=theta,
                colorscale='RdBu',
                zmin=vmin,
                zmax=vmax,
                colorbar=dict(title='Wigner value'),
                customdata=custom,
                hovertemplate=(
                    'φ: %{x:.4f}<br>'
                    'θ: %{y:.4f}<br>'
                    'W: %{z:.4f}<extra></extra>'
                )
            ))

            heat.update_layout(
                title=title + (f' (j={j})' if j is not None else ''),
                xaxis_title='φ (rad)',
                yaxis_title='θ (rad)',
                height=500,
                width=900
            )
            heat.update_xaxes(tickmode='array', tickvals=phi_ticks.tolist(), ticktext=phi_labels)
            heat.update_yaxes(tickmode='array', tickvals=theta_ticks.tolist(), ticktext=theta_labels)

            if show_reference:
                heat.add_trace(go.Scatter(
                    x=[0.0], y=[np.pi/2.0],
                    mode='markers',
                    marker=dict(color='red', size=8),
                    name='ref'
                ))

            html_name_plane = os.path.join(output_folder, fname_base + '_wigner_plane_interactive.html')
            heat.write_html(html_name_plane, include_plotlyjs='cdn')
            try:
                heat.show()
            except Exception:
                pass
        except Exception:
            print('Plotly not installed; skipping interactive plane.')

    return plane_png, html_name_plane


def plot_logical_wigner_func(encoding_superoperator, out_prefix='logical_states', photon_number=None, interactive=True, output_folder=None, pi_step=np.pi/4, titles=None, projector_only=False):
    """
    High-level orchestrator: computes wigners for projector and logical states,
    then calls the small plotting helpers. Returns a dict of saved filenames.

    ``titles`` may be a dictionary overriding any of ``projector_sphere``,
    ``projector_plane``, ``logical{idx}_sphere``, or ``logical{idx}_plane``.
    Missing keys retain their existing default titles.

    Set ``projector_only=True`` to compute and plot only the code projector.

    You can call compute_wigners_for_encoding(), plot_wigner_sphere(), and
    plot_wigner_plane() individually if you prefer.
    """
    if titles is None:
        titles = {}
    elif not isinstance(titles, dict):
        raise TypeError('titles must be a dictionary or None')

    results = compute_wigners_for_encoding(encoding_superoperator, projector_only=projector_only)
    j = results['j']
    d = results['d']

    if output_folder:
        folder = output_folder
    else:
        folder = out_prefix + '_Wigner_Plots'
    os.makedirs(folder, exist_ok=True)
    base = os.path.basename(out_prefix)

    saved = {}

    # Projector
    THp, PHp, Wp = results['projector']
    fname_base_proj = base + '_projector'
    projector_sphere_title = titles.get('projector_sphere', 'Projector Spin Wigner')
    projector_plane_title = titles.get('projector_plane', 'Projector Spin Wigner Function')
    saved['projector_sphere'] = plot_wigner_sphere(THp, PHp, Wp, fname_base_proj, projector_sphere_title, j, d=d, interactive=interactive, output_folder=folder, flip_phi=False, show_reference=True)
    saved['projector_plane'] = plot_wigner_plane(THp, PHp, Wp, fname_base_proj, projector_plane_title, j, d=d, interactive=interactive, output_folder=folder, pi_step=pi_step, show_reference=True)

    # Logical states (support arbitrary d)
    logical_keys = sorted([k for k in results.keys() if k.startswith('logical')], key=lambda s: int(s.replace('logical','')))
    for key in logical_keys:
        THi, PHi, Wi = results[key]
        idx = int(key.replace('logical',''))
        fname_base_i = f"{base}_logical{idx}"
        sphere_title = titles.get(f'logical{idx}_sphere', f'Logical |{idx}> Spin Wigner')
        plane_title = titles.get(f'logical{idx}_plane', f'Logical |{idx}> Spin Wigner (plane)')
        saved[f'logical{idx}_sphere'] = plot_wigner_sphere(THi, PHi, Wi, fname_base_i, sphere_title, j, interactive=interactive, output_folder=folder, flip_phi=False, show_reference=True)
        saved[f'logical{idx}_plane'] = plot_wigner_plane(THi, PHi, Wi, fname_base_i, plane_title, j, interactive=interactive, output_folder=folder, pi_step=pi_step, show_reference=True)

    print('Saved Wigner outputs to', folder)
    return saved



    
 
def spin_projector_from_encoding(encoding_superoperator):
    """
    Build the optimized code projector in the spin basis used by the irrep projectors.
    Existing uses of `projector(logical_states(enc))` are unchanged elsewhere.
    """
    spin_logicals = [schwinger_to_spin(v) for v in logical_states(encoding_superoperator)]
    return projector(spin_logicals)







def plot_fidelities_vs_eta(etas, fidelities, source_file, out_suffix='_fidelity'):
    """Plot fidelity vs eta and save to <source_file_without_ext>_fidelity.png

    etas: list or array of eta values
    fidelities: list or array of fidelities corresponding to etas
    source_file: the path to the pickle/data file used (string)
    out_suffix: suffix to append to source filename (default '_fidelity')
    """
    try:
        plt.style.use('seaborn-v0_8-whitegrid')
    except Exception:
        pass

    dark_blue = '#0b3d91'
    purple = '#6a0dad'

    etas_arr = np.array(etas)
    fes_arr = np.array(fidelities)

    fig, ax = plt.subplots(figsize=(7,4))
    # solid line with filled markers, no contrasting outline
    ax.plot(etas_arr, fes_arr, marker='o', linestyle='-', color=dark_blue, linewidth=2,
        markersize=6, markerfacecolor=dark_blue, markeredgecolor='none')

    ax.set_xlabel(r'Loss transmissivity $\eta$ (eta)')
    ax.set_ylabel('Entanglement Fidelity')
    pn = extract_photon_number_from_filename(source_file)
    if pn is not None:
        ax.set_title(f'Entanglement Fidelity vs. eta (n={pn})')
    else:
        ax.set_title('Entanglement Fidelity vs. eta')
    ax.set_xlim(min(etas_arr) - 0.01, max(etas_arr) + 0.01)
    ax.set_ylim(0.0, 1.02)
    ax.grid(True, linestyle=':', linewidth=0.6)

    # annotate the numeric values for readability
    for x, y in zip(etas_arr, fes_arr):
        ax.text(x, y + 0.02, f'{y:.4f}', ha='center', va='bottom', fontsize=8, color='#222222')

    base = source_file[:-4] if source_file.endswith('.pkl') else source_file
    outname = base + out_suffix + '.png'
    fig.tight_layout()
    fig.savefig(outname, dpi=150)
    plt.close(fig)


def plot_coherent_information_vs_eta(etas, ics, source_file, out_suffix='_coherent_information'):
    """Plot coherent information vs eta and save to <source_file_without_ext>_coherent_information.png"""
    try:
        plt.style.use('seaborn-v0_8-whitegrid')
    except Exception:
        pass

    dark_blue = '#0b3d91'
    purple = '#6a0dad'

    etas_arr = np.array(etas)
    ics_arr = np.array(ics)

    fig, ax = plt.subplots(figsize=(7,4))
    # solid line with filled markers, no contrasting outline
    ax.plot(etas_arr, ics_arr, marker='s', linestyle='-', color=purple, linewidth=2,
        markersize=6, markerfacecolor=purple, markeredgecolor='none')

    ax.set_xlabel(r'Loss transmissivity $\eta$ (eta)')
    ax.set_ylabel('Coherent information (I_c)')
    pn = extract_photon_number_from_filename(source_file)
    if pn is not None:
        ax.set_title(f'Coherent information vs. eta (n={pn})')
    else:
        ax.set_title('Coherent information vs. eta')
    ax.grid(True, linestyle=':', linewidth=0.6)

    # horizontal zero line for reference
    ax.axhline(0.0, color='gray', linestyle='--', linewidth=1.0)

    # annotate points
    for x, y in zip(etas_arr, ics_arr):
        ax.text(x, y + 0.02 if y >= 0 else y - 0.03, f'{y:.4f}', ha='center', va='bottom' if y >= 0 else 'top', fontsize=8, color='#222222')

    base = source_file[:-4] if source_file.endswith('.pkl') else source_file
    outname = base + out_suffix + '.png'
    fig.tight_layout()
    fig.savefig(outname, dpi=150)
    plt.close(fig)





def load_single_file_data(folder,file):
    full_path = os.path.join(folder, file)
    file = open(full_path, 'rb')
    data = pickle.load(file)
    file.close()
    return data


# Test performance of a code for different noise strength 
# after optimizing decoding
def test_code(file, new_eta, eps=1e-8, *, return_optimizer=False):
    """Evaluate a saved encoder with a freshly optimized decoder; do not save.

    The encoder remains fixed. Return (Fe, I_c), or the evaluated optimizer
    when return_optimizer=True so callers can inspect sector diagnostics.
    eps sets the SCS tolerance and strict validation targets; a rejected solve
    raises instead of returning a fallback decoder's performance.
    """
    if not np.isfinite(new_eta) or not 0 <= new_eta <= 1:
        raise ValueError('new_eta must be between zero and one')
    if not np.isfinite(eps) or eps <= 0:
        raise ValueError('eps must be positive and finite')
    with open(file, 'rb') as fh:
            data = pickle.load(fh)
            n, eta, d, rounds, seed, dec, enc, fe, I_c = data
    
    loss_op = TruncatedLoss(n, new_eta)
    optimizer = SDPCodeOptimizer(loss_op, enc, n, d, eps=eps)
    optimizer.strict_solver_mode = True
    optimizer.target_gap = eps
    optimizer.target_residual = 10 * eps
    optimizer.target_physical = 10 * eps
    for info in optimizer.decoder_solver_info.values():
        info['max_iters'] = 100000
    if not optimizer.solve_decoder_problem():
        raise RuntimeError('Decoder optimization did not pass validation')
    optimizer.calculate_performance()

    return optimizer if return_optimizer else (optimizer.Fe, optimizer.I_c)




def plot_file_simple(
    pkl_path,
    out_folder=None,
    interactive=True,
    rotate=False,
    theta=0,
    phi=0,
    titles=None,
    reflect=False,
    alignment_euler=None,
    projector_only=False,
):
    """
    Plot Wigner (projector) and projector heatmaps for a single .pkl result file.
    Expects the simple tuple format: (n, eta, d, rounds, seed, dec, enc, fe, I_c).
    If out_folder is None, creates <pkl_basename>_Wigner_Plots next to the pickle.
    If ``alignment_euler=(alpha, beta, gamma)`` is supplied, first applies that
    fixed z-y-z Euler pre-alignment. If ``rotate=True``, the existing visual
    rotation specified by the sphere coordinates ``theta`` and ``phi`` is then
    applied unchanged. If ``reflect=True``, complex-conjugates the code before
    both rotations. Thus the transformation order is reflection, fixed Euler
    pre-alignment, and visual ``theta``, ``phi`` rotation.
    ``titles`` is an optional dictionary whose keys may be ``projector_sphere``,
    ``projector_plane``, ``logical{idx}_sphere``, or ``logical{idx}_plane``.
    Set ``projector_only=True`` to omit logical-state Wigner computation and plots.
    The projector sphere and plane are still saved as PNGs, plus interactive HTML
    when ``interactive=True``. The default retains the full standard collection.
    Returns the dictionary of paths produced by the Wigner plot functions.
    """
    import pickle
    if not os.path.isfile(pkl_path):
        print(f'File not found: {pkl_path}')
        return None

    try:
        with open(pkl_path, 'rb') as fh:
            data = pickle.load(fh)
    except Exception as e:
        print(f'Failed to load {pkl_path}: {e}')
        return None

    try:
        n, eta, d, rounds, seed, dec, enc, fe, I_c = data
        print(n,eta, fe, I_c, rounds)
    except Exception:
        print(f'{pkl_path} does not match expected tuple format, skipping.')
        return None

    if reflect:
        enc = reflect_code(enc)

    if alignment_euler is not None:
        try:
            alpha, beta, gamma = alignment_euler
        except (TypeError, ValueError) as exc:
            raise ValueError(
                'alignment_euler must contain exactly three angles '
                '(alpha, beta, gamma)'
            ) from exc
        enc = rotate_code_euler(
            enc,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
        )

    if rotate:
        try:
            enc = rotate_code(enc, theta=theta, phi=phi)
        except Exception as e:
            print(f'Rotation failed for {pkl_path}: {e}')
            # Continue with the pre-aligned, reflected, or original encoder.

    pn = int(n) if n is not None else None
    base = os.path.splitext(os.path.basename(pkl_path))[0]
    if out_folder is None:
        out_folder = os.path.join(os.path.dirname(pkl_path) or '.', base + '_Wigner_Plots')
    os.makedirs(out_folder, exist_ok=True)

    wigner_prefix = os.path.join(out_folder, base + '_logical_states')
    
    wigner_saved = plot_logical_wigner_func(enc, out_prefix=wigner_prefix, photon_number=pn, interactive=interactive, output_folder=out_folder, titles=titles, projector_only=projector_only)
    print(f'Wigner saved for {base}: {wigner_saved}')
   



    return wigner_saved


def batch_plot_folder_simple(results_folder='0.95 results', interactive=True, rotate=False, theta=0, phi=0):
    """
    Simple batch runner that calls plot_file_simple on every .pkl in results_folder.
    Keeps behaviour minimal: looks for the expected tuple format and skips invalid files.
    """
    import pickle

    if not os.path.isdir(results_folder):
        print(f'Folder not found: {results_folder}')
        return

    files = sorted([f for f in os.listdir(results_folder) if f.endswith('.pkl')])
    if not files:
        print(f'No .pkl files in {results_folder}')
        return

    for fname in files:
        full = os.path.join(results_folder, fname)
        print(f'Processing {fname} ...')
        plot_file_simple(full, out_folder=os.path.join(results_folder, os.path.splitext(fname)[0] + '_Wigner_Plots'),
                         interactive=interactive, rotate=rotate, theta=theta, phi=phi)

    print('Done.')



def gather_performance_from_folders(folder_list):
    # Use n as keys.
    # For each n, keep one entry per eta: the d=2 result with the highest fidelity.
    performance_dict = {}

    for folder in folder_list:
        files = sorted([f for f in os.listdir(folder) if f.endswith('.pkl')])
        for file in files:
            full_path = os.path.join(folder, file)
            with open(full_path, 'rb') as fh:
                data = pickle.load(fh)
                n, eta, d, rounds, seed, dec, enc, fe, I_c = data

                # Only include qubit codes; skip higher-dimensional logical codes.
                if d != 2:
                    continue

                if n not in performance_dict:
                    performance_dict[n] = {}

                # For each (n, eta), keep the run with the highest fidelity.
                if eta not in performance_dict[n] or fe > performance_dict[n][eta][0]:
                    performance_dict[n][eta] = (fe, I_c)

    return {
        n: [(eta, fe, I_c) for eta, (fe, I_c) in eta_dict.items()]
        for n, eta_dict in performance_dict.items()
    }



def plot_all_fidelities(folder_list, out_folder='Joint_Performance_Plots'):
    performance_dict = gather_performance_from_folders(folder_list)
    if not performance_dict:
        print('No performance data found for given folders.')
        return

    try:
        plt.style.use('seaborn-v0_8-whitegrid')
    except Exception:
        pass

    os.makedirs(out_folder, exist_ok=True)

    # build a reproducible ordered list of n values
    n_keys = sorted(performance_dict.keys())

    # get a pretty palette (seaborn if available, fallback to matplotlib)
    try:
        import seaborn as sns
        palette = sns.color_palette("husl", len(n_keys))  # pleasant distinct colors
    except Exception:
        cmap = plt.get_cmap('tab10')
        palette = [cmap(i % 10) for i in range(len(n_keys))]

    markers = ['o', 's', 'D', '^', 'v', 'P', 'X', '<', '>', '*']  # cycle if many n

    for quantity, ylabel, out_suffix in [
        ('fe', 'Entanglement Fidelity', '_fidelity'),
        ('I_c', 'Coherent Information', '_coherent_information')
    ]:
        fig, ax = plt.subplots(figsize=(12, 8))

        all_etas = []
        all_vals = []

        # diagnostic summary to help find suspicious data
        print(f'--- Summary for {quantity} ---')

        for idx, n in enumerate(n_keys):
            data = performance_dict[n]
            data = sorted(data, key=lambda t: t[0])
            etas = np.array([item[0] for item in data])
            if etas.size == 0:
                continue
            if quantity == 'fe':
                vals = np.array([item[1] for item in data])
            else:
                vals = np.array([item[2] for item in data])

            # basic diagnostics
            try:
                vmin_n = float(np.nanmin(vals))
                vmax_n = float(np.nanmax(vals))
                med_n = float(np.nanmedian(vals))
            except Exception:
                vmin_n = vmax_n = med_n = float('nan')
            print(f'n={n}: η range [{etas.min():.3f},{etas.max():.3f}]  {quantity} min/med/max = {vmin_n:.6g} / {med_n:.6g} / {vmax_n:.6g}')

            all_etas.append(etas)
            all_vals.append(vals)

            color = palette[idx]
            marker = markers[idx % len(markers)]
            ax.plot(etas, vals,
                    marker=marker, linestyle='-',
                    color=color, linewidth=2,
                    markersize=6, markerfacecolor=color, markeredgecolor='none',
                    label=f'n={n}')

        if not all_etas:
            print(f'No data to plot for {quantity}')
            plt.close(fig)
            continue

        # compute global data ranges across all n
        all_etas_flat = np.concatenate(all_etas)
        all_vals_flat = np.concatenate(all_vals)

        x_min, x_max = float(all_etas_flat.min()), float(all_etas_flat.max())
        x_pad = max(1e-6, (x_max - x_min) * 0.02)
        ax.set_xlim(x_min - x_pad, x_max + x_pad)

        # -- tighten y-limits to the actual data with a small, controlled padding --
        y_min_data, y_max_data = float(all_vals_flat.min()), float(all_vals_flat.max())
        data_range = max(y_max_data - y_min_data, 0.0)
        pad = max(1e-4, data_range * 0.03)

        ymin = y_min_data - pad
        ymax = y_max_data + pad

        if quantity == 'fe':
            ymin = max(0.0, ymin)
            ymax = min(1.0, ymax)
            if ymax - ymin < 1e-3:
                center = (ymax + ymin) / 2.0
                ymin = max(0.0, center - 5e-4)
                ymax = min(1.0, center + 5e-4)
        else:
            if ymax - ymin < 1e-6:
                ymin -= 1e-3
                ymax += 1e-3

        ax.set_ylim(ymin, ymax)

        ax.set_xlabel(r'Loss transmissivity $\eta$', fontsize=24)
        ax.set_ylabel(ylabel, fontsize=24, labelpad=18)
        ax.set_title(f'{ylabel} vs. η for various n', fontsize=30, pad=20)
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.grid(True, linestyle=':', linewidth=0.6)
        ax.legend(title='n', loc='best', frameon=True, fontsize=18, title_fontsize=20)

        # inset zoom: larger inset placed a bit higher; plot only for eta in [0.85,1.0]
        try:
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            # bigger inset (~1.5-2x previous small versions), and shifted up
            axins = inset_axes(ax, width="92%", height="92%",
                               bbox_to_anchor=(0.52, 0.12, 0.48, 0.48),
                               bbox_transform=ax.transAxes, loc='lower left')

            eps = 1e-16
            inset_eta_min = 0.85
            inset_eta_max = 1.0

            any_plotted = False
            for idx, n in enumerate(n_keys):
                if idx >= len(all_etas):
                    continue
                etas = all_etas[idx]
                vals = all_vals[idx]
                mask = (etas >= inset_eta_min) & (etas <= inset_eta_max)
                if not np.any(mask):
                    continue
                etas_sub = etas[mask]
                vals_sub = vals[mask]
                color = palette[idx]

                if quantity == 'fe':
                    y = 1.0 - vals_sub
                else:
                    y = np.abs(1.0 - vals_sub)
                # clip to positive for log scale (use very small eps)
                y = np.clip(y, eps, None)
                axins.plot(etas_sub, y, marker='.', linestyle='-', color=color, markersize=5)
                any_plotted = True

            if any_plotted:
                axins.set_yscale('log')
                axins.set_xlim(inset_eta_min, inset_eta_max)

                # build tight inset y-limits from actual inset data and add modest multiplicative padding
                all_y_for_inset = []
                for idx in range(len(all_etas)):
                    etas = all_etas[idx]
                    vals = all_vals[idx]
                    mask = (etas >= inset_eta_min) & (etas <= inset_eta_max)
                    if not np.any(mask):
                        continue
                    vals_sub = vals[mask]
                    if quantity == 'fe':
                        all_y_for_inset.append(np.clip(1.0 - vals_sub, eps, None))
                    else:
                        all_y_for_inset.append(np.clip(np.abs(1.0 - vals_sub), eps, None))

                if all_y_for_inset:
                    all_y_flat = np.concatenate(all_y_for_inset)
                    ymin_raw = float(all_y_flat.min())
                    ymax_raw = float(all_y_flat.max())
                    # multiplicative padding that fits data tightly but leaves breathing room
                    pad_factor = 2.0
                    ymin_ins = max(eps, ymin_raw / pad_factor)
                    ymax_ins = max(ymax_raw * pad_factor, ymin_ins * 10.0)
                    # ensure valid order
                    if ymax_ins <= ymin_ins:
                        ymax_ins = ymin_ins * 10.0 + eps
                    axins.set_ylim(ymin_ins, ymax_ins)

                axins.grid(True, linestyle=':', linewidth=0.4)
                if quantity == 'fe':
                    axins.set_ylabel('1 - fe (log scale)', fontsize=18)
                else:
                    axins.set_ylabel('|1 - I_c| (log scale)', fontsize=18)
                axins.tick_params(axis='both', which='major', labelsize=16)
            else:
                axins.remove()
        except Exception as e:
            print('Could not create inset axes:', e)

        base = 'all_n_results'
        outname = os.path.join(out_folder, base + out_suffix + '.png')
        fig.tight_layout()
        fig.savefig(outname, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved {outname}')



def update_file(path, data):
    with open(path, "wb") as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

# Rewrite this for arbitrary d, not just qubits
def get_superoperator_from_logical_states(*logical_states):
    logical_states = [np.asarray(state) for state in logical_states]
    new_isometry = np.column_stack(logical_states)
    new_superoperator = np.kron(new_isometry, new_isometry.conj())

    return new_superoperator


def plot_literature_code(
    n,
    author,
    group=None,
    instance=None,
    out_folder=None,
    interactive=True,
    titles=None,
    eta=None,
    eps=1e-9,
):
    """Construct and plot a selected two-mode code from the literature.

    Catalog keys have the form ``(n, author, group)``.  When a source gives
    more than one inequivalent code with the same first three fields, later
    instances use ``(n, author, group, instance)``; the unnumbered key is the
    default first instance.  Author and group labels are case-insensitive,
    and underscores in group labels are ignored.

    A vector entry with index ``n_a`` is the coefficient of the Fock state
    ``|n_a, n-n_a>``.  The AAB label is shared here by the
    Aydin--Albert--Barg entries retained from the previous helper and the
    Aydin--Alekseyev--Barg ``n=19`` and ``n=22`` permutation-invariant
    codes; their group field is ``None`` because these entries do not carry
    a finite-subgroup attribution.  The KT entries are numerical outputs of the
    Kubischta--Teixeira Mathematica notebook, converted from its
    descending-``m`` spin basis to this file's ascending-``n_a`` basis.  The
    OG entry is the approximate ``n=25`` binary-octahedral code of
    Omanakuttan--Gross. The CLY entry is the ``n=4`` Chuang--Leung--Yamamoto
    code, with its binary-octahedral ``2O`` symmetry.

    Examples
    --------
    ``plot_literature_code(4, 'CLY', '2O', eta=0.9)``
    ``plot_literature_code(9, 'AAB')``
    ``plot_literature_code(19, 'AAB')``
    ``plot_literature_code(22, 'AAB')``
    ``plot_literature_code(9, 'KT', '2D4')``
    ``plot_literature_code(21, 'KT', '2D_4', instance=2)``
    ``plot_literature_code(25, 'OG', '2O')``
    ``plot_literature_code(25, 'OG', '2O', eta=0.9)``

    Returns the selected catalog key, codewords, code projector, encoding
    superoperator, dictionary of saved Wigner-plot paths, and optionally the
    fixed-encoder performance after one decoder optimization at ``eta``.
    Supplying ``eta`` also saves the standard nine-field result pickle in
    ``out_folder``, named ``<base>_eta<eta>.pkl``. Its path is returned as
    ``code_file`` (None when ``eta`` is omitted).
    ``eps`` defaults to 1e-9 and sets the strict decoder-validation targets,
    as in test_code; rejected solves raise instead of saving a result.
    """

    def _state(total_photons, amplitudes):
        state = np.zeros(total_photons + 1, dtype=complex)
        for n_a, amplitude in amplitudes.items():
            state[n_a] = amplitude
        return state

    # Aydin--Alekseyev--Barg, Eq. (31) and the approximate q_{2i}
    # coefficients on p. 25.  Dicke weight w maps to n_a=w, and Eq. (31)
    # supplies the additional sqrt(binomial(19, w)) factor.  The source states
    # are unnormalized and the printed q_{2i} are rounded, so normalize after
    # assembling the even-weight word and obtain the odd-weight word by the
    # reflection prescribed in Eq. (31).
    aab_pi_n19_q = np.array(
        [
            1.0,
            0.0477572,
            -0.0267249,
            -0.00506367,
            0.00332914,
            0.00527235,
            -0.000947223,
            0.0152707,
            0.00888631,
            0.32678,
        ],
        dtype=float,
    )
    aab_pi_n19_logical_0 = _state(
        19,
        {
            weight: coefficient * np.sqrt(math.comb(19, weight))
            for weight, coefficient in zip(range(0, 20, 2), aab_pi_n19_q)
        },
    )
    aab_pi_n19_logical_0 /= np.linalg.norm(aab_pi_n19_logical_0)
    aab_pi_n19_logical_1 = aab_pi_n19_logical_0[::-1].copy()

    # Omanakuttan--Gross, Eqs. (75)--(76).  Equation (76) prints the last two
    # magnetic labels as -13/2 and -21/2.  Those labels break the required
    # spacing-by-four support and do not define a 2O-invariant code.  The
    # symmetry-consistent labels are -15/2 and -23/2, corresponding here to
    # n_a=5 and n_a=1.  The paper gives approximate coefficients, so normalize
    # the assembled state before constructing its J_z-reflected logical mate.
    og_n25_multiplicity_basis = np.array(
        [
            [
                -np.sqrt(1377 / 4132),
                -np.sqrt(1 / 674),
                -np.sqrt(109 / 1169),
                -np.sqrt(803 / 1918),
                -np.sqrt(103 / 690),
                -np.sqrt(1 / 263),
                -np.sqrt(1 / 3608),
            ],
            [
                np.sqrt(1 / 4402),
                -np.sqrt(2 / 839),
                -np.sqrt(293 / 983),
                -np.sqrt(11 / 1264),
                np.sqrt(913 / 2925),
                np.sqrt(21 / 412),
                -np.sqrt(1069 / 3264),
            ],
            [
                -np.sqrt(1 / 61408),
                np.sqrt(1750 / 2781),
                -np.sqrt(325 / 3548),
                np.sqrt(43 / 763),
                -np.sqrt(47 / 551),
                np.sqrt(183 / 1349),
                np.sqrt(2 / 1011),
            ],
        ],
        dtype=float,
    )
    og_n25_multiplicity_coefficients = np.array(
        [
            -np.sqrt(267 / 1213),
            np.sqrt(701 / 1457),
            np.sqrt(337 / 1128),
        ],
        dtype=float,
    )
    og_n25_amplitudes = (
        og_n25_multiplicity_coefficients @ og_n25_multiplicity_basis
    )
    og_n25_logical_0 = _state(
        25,
        dict(zip((25, 21, 17, 13, 9, 5, 1), og_n25_amplitudes)),
    )
    og_n25_logical_0 /= np.linalg.norm(og_n25_logical_0)
    og_n25_logical_1 = og_n25_logical_0[::-1].copy()

    literature_codes = {
        # Chuang--Leung--Yamamoto, Phys. Rev. A 56, 1114 (1997), Eq. (3.1).
        (4, 'CLY', '2O'): (
            _state(4, {0: 1 / np.sqrt(2), 4: 1 / np.sqrt(2)}),
            _state(4, {2: 1.0}),
        ),
        # Aydin--Albert--Barg codes retained from the previous helper.
        (9, 'AAB', None): (
            _state(9, {9: 1 / 2, 3: np.sqrt(3) / 2}),
            _state(9, {6: np.sqrt(3) / 2, 0: 1 / 2}),
        ),
        (11, 'AAB', None): (
            _state(11, {0: np.sqrt(5) / 4, 8: np.sqrt(11) / 4}),
            _state(11, {3: np.sqrt(11) / 4, 11: np.sqrt(5) / 4}),
        ),
        # Aydin--Alekseyev--Barg ((19, 2, 5)) PI code from p. 25.
        (19, 'AAB', None): (
            aab_pi_n19_logical_0,
            aab_pi_n19_logical_1,
        ),
        (21, 'AAB', None): (
            _state(
                21,
                {
                    0: np.sqrt(5 / 68),
                    8: np.sqrt(7 / 12),
                    17: np.sqrt(35 / 102),
                },
            ),
            _state(
                21,
                {
                    4: np.sqrt(35 / 102),
                    13: -np.sqrt(7 / 12),
                    21: -np.sqrt(5 / 68),
                },
            ),
        ),
        # Aydin--Alekseyev--Barg Q_{4,2,5,-}, obtained by substituting
        # (g, m, delta, epsilon) = (4, 2, 5, -1) in Construction 5.1.
        # Proposition 5.5 shows that it corrects four deletions.
        (22, 'AAB', None): (
            _state(
                22,
                {
                    0: np.sqrt(1 / 12),
                    8: np.sqrt(11 / 20),
                    18: np.sqrt(11 / 30),
                },
            ),
            _state(
                22,
                {
                    4: np.sqrt(11 / 30),
                    14: -np.sqrt(11 / 20),
                    22: -np.sqrt(1 / 12),
                },
            ),
        ),
        # KT BD_4, tau_1 branch matching the supplied Mathematica output.
        # This branch is reproduced by SeedRandom[3].
        (9, 'KT', '2D4'): (
            _state(
                9,
                {
                    0: 0.5023441438459086,
                    4: -0.6085233950206410,
                    8: -0.6142879120219393,
                },
            ),
            _state(
                9,
                {
                    1: -0.6142879120219394,
                    5: -0.6085233950206411,
                    9: 0.5023441438459088,
                },
            ),
        ),
        # Kubischta--Teixeira BD_4, tau_2 code from SeedRandom[1].
        (21, 'KT', '2D4'): (
            _state(
                21,
                {
                    1: 0.4409063836264944,
                    5: -0.4426807032398304,
                    9: 0.3483203759716437,
                    13: 0.3590980152295407,
                    17: 0.4972514075055183,
                    21: 0.3348102816991354,
                },
            ),
            _state(
                21,
                {
                    0: 0.3348102816991353,
                    4: 0.4972514075055182,
                    8: 0.3590980152295406,
                    12: 0.3483203759716436,
                    16: -0.4426807032398303,
                    20: 0.4409063836264943,
                },
            ),
        ),
        # Inequivalent KT tau_2 notebook solution from SeedRandom[12].
        (21, 'KT', '2D4', 2): (
            _state(
                21,
                {
                    1: 0.4311755035059772,
                    5: -0.3633073493147888,
                    9: -0.5147975457890357,
                    13: 0.2128749665776041,
                    17: -0.5308375112757139,
                    21: -0.2999578756875792,
                },
            ),
            _state(
                21,
                {
                    0: -0.2999578756875792,
                    4: -0.5308375112757139,
                    8: 0.2128749665776041,
                    12: -0.5147975457890357,
                    16: -0.3633073493147888,
                    20: 0.4311755035059772,
                },
            ),
        ),
        # Omanakuttan--Gross [[25, 1, 5]] binary-octahedral code.
        (25, 'OG', '2O'): (
            og_n25_logical_0,
            og_n25_logical_1,
        ),
    }

    author_key = str(author).upper()
    if group is None:
        group_key = None
    else:
        group_key = str(group).upper().replace('_', '').replace(' ', '')

    base_key = (n, author_key, group_key)
    if instance in (None, 1):
        code_key = base_key
    else:
        code_key = base_key + (instance,)

    if code_key not in literature_codes:
        available = ', '.join(repr(key) for key in literature_codes)
        raise ValueError(
            f'Unknown literature code {code_key!r}. '
            f'Available keys are: {available}'
        )

    logical_0, logical_1 = literature_codes[code_key]

    code_projector = projector(logical_0, logical_1)
    encoding_superoperator = get_superoperator_from_logical_states(
        logical_0,
        logical_1,
    )

    base_parts = [f'n{n}', author_key.lower()]
    if group_key is not None:
        base_parts.append(group_key.lower())
    if len(code_key) == 4:
        base_parts.append(str(code_key[3]))
    base = 'literature_' + '_'.join(base_parts)
    if out_folder is None:
        out_folder = base + '_Wigner_Plots'
    os.makedirs(out_folder, exist_ok=True)

    performance = None
    code_file = None
    if eta is not None:
        code_file = os.path.join(out_folder, f'{base}_eta{float(eta):.6f}.pkl')
        _, fe, coherent_information = get_performance(
            n,
            2,
            encoding_superoperator,
            float(eta),
            reencode=False,
            file_output=code_file,
            eps=eps,
        )
        performance = {
            'eta': float(eta),
            'fe': float(fe),
            'I_c': float(coherent_information),
        }
        print(
            f'Performance for {code_key} at eta={float(eta):.6f}: '
            f'Fe={float(fe):.12f}, I_c={float(coherent_information):.12f}'
        )
        print(f'Code saved to {code_file}')

    wigner_prefix = os.path.join(out_folder, base + '_logical_states')
    wigner_saved = plot_logical_wigner_func(
        encoding_superoperator,
        out_prefix=wigner_prefix,
        photon_number=n,
        interactive=interactive,
        output_folder=out_folder,
        titles=titles,
    )
    print(f'Wigner saved for {base}: {wigner_saved}')

    return {
        'n': n,
        'code_key': code_key,
        'logical_states': (logical_0, logical_1),
        'projector': code_projector,
        'encoding_superoperator': encoding_superoperator,
        'performance': performance,
        'code_file': code_file,
        'plots': wigner_saved,
    }



def get_performance(n,d, encoding, eta, reencode=False, file_output=None, *, eps=1e-9):
    """Evaluate with strict solver validation; optionally save a result pickle.

    eps sets the SCS tolerance, gap target, and residual/physicality targets
    as in test_code. A rejected solve raises before evaluation or saving.
    The encoder stays fixed unless reencode=True requests one encoder update.
    """
    if not np.isfinite(eta) or not 0 <= eta <= 1:
        raise ValueError('eta must be between zero and one')
    if not np.isfinite(eps) or eps <= 0:
        raise ValueError('eps must be positive and finite')
    loss_op = TruncatedLoss(n, eta)
    optimizer = SDPCodeOptimizer(loss_op, encoding, n, d, eps=eps, file_output=file_output)
    optimizer.strict_solver_mode = True
    optimizer.target_gap = eps
    optimizer.target_residual = 10 * eps
    optimizer.target_physical = 10 * eps
    for info in optimizer.decoder_solver_info.values():
        info['max_iters'] = 100000
    if not optimizer.solve_decoder_problem():
        raise RuntimeError('Decoder optimization did not pass validation')

    if reencode:
        optimizer.encoder_solver_info['max_iters'] = 100000
        if not optimizer.solve_encoder_problem():
            raise RuntimeError('Encoder optimization did not pass validation')
    
    optimizer.calculate_performance()
    if file_output is not None:
        optimizer.save()

    return optimizer.T_E, optimizer.Fe, optimizer.I_c







def projector_basis(projector_matrix, rank=None, tol=1e-8):
    eigvals, eigvecs = np.linalg.eigh(projector_matrix)
    if rank is None:
        keep = eigvals > tol
    else:
        idx = np.argsort(eigvals)[-rank:]
        keep = np.zeros_like(eigvals, dtype=bool)
        keep[idx] = True
    basis = eigvecs[:, keep]
    return basis




def spin_projector_from_result_file(codefile):
    folder = os.path.dirname(codefile) or '.'
    filename = os.path.basename(codefile)
    n, eta, d, rounds, seed, dec, enc, fe, I_c = load_single_file_data(folder, filename)
    return spin_projector_from_encoding(enc), {
        'n': n,
        'eta': eta,
        'd': d,
        'rounds': rounds,
        'seed': seed,
        'fe': fe,
        'I_c': I_c,
    }



def _support_projectors_from_files(file_a, file_b):
    """Load compatible saved codes as exact support projectors and their rank."""
    projectors = []
    parameters = []
    for path in (file_a, file_b):
        raw, metadata = spin_projector_from_result_file(os.fspath(path))
        n, rank = metadata['n'], metadata['d']
        if raw.shape != (n + 1, n + 1) or not 1 <= rank <= n + 1:
            raise ValueError('Saved photon number, rank and support dimensions disagree')
        basis = projector_basis((raw + raw.conj().T) / 2, rank=rank)
        projectors.append(basis @ basis.conj().T)
        parameters.append((n, rank))
    if parameters[0] != parameters[1]:
        raise ValueError('The two files must have the same photon number and code rank')
    return projectors, parameters[0][1]


def projector_overlap_from_files(
    file_a, file_b, *,
    alignment_euler_a=None, alignment_euler_b=None,
    reflect_a=False, reflect_b=False,
):
    """Return normalized support-projector overlap for two saved codes.

    Both files must have the same photon number and code rank. As in the group
    analysis, each saved support is replaced by its exact rank-d projector.
    Each file is independently conjugated first, then rotated by its supplied
    (alpha, beta, gamma) angles in the plot_file_simple convention (radians).
    None means no rotation. Display-level theta/phi rotations are not applied.
    No file is modified and no alignment or decoder optimization is performed.
    Use find_projector_alignment_from_files to search for an alignment.
    """
    supports, rank = _support_projectors_from_files(file_a, file_b)
    projectors = []
    for support, angles, reflected in zip(
        supports, (alignment_euler_a, alignment_euler_b), (reflect_a, reflect_b),
    ):
        if reflected:
            support = support.conj()
        if angles is not None:
            angles = np.asarray(angles, dtype=float)
            if angles.shape != (3,) or not np.all(np.isfinite(angles)):
                raise ValueError('Each Euler rotation must contain three finite angles')
            support = rotate_projector_euler(support, *angles)
        projectors.append(support)
    return float(np.trace(projectors[0] @ projectors[1]).real / rank)


def find_projector_alignment(
    projector_a, projector_b, *, allow_reflection=True,
    seeds=(0, 1, 2), maxiter=100, popsize=12, initial_alignments=(),
):
    """Find a high-overlap rotation of projector_a onto fixed projector_b.

    Inputs are equal-rank support projectors in the occupation/spin basis.
    The returned alignment_euler and reflect apply to A, in plot_file_simple's
    convention: conjugation followed by the fixed Euler rotation. B stays fixed;
    pretransform B before calling if a particular canonical orientation is wanted.

    initial_alignments is a sequence of dictionaries with alignment_euler and
    reflect keys, for example a previous return value. Each candidate is evaluated
    before polishing, so a poorer search cannot replace it. seeds=() skips the
    global searches. This is a best-found alignment, not a certified maximum.
    The final overlap is recomputed with rotate_projector_euler.
    """
    projectors = [np.asarray(p, dtype=complex) for p in (projector_a, projector_b)]
    bases = []
    for p in projectors:
        if p.ndim != 2 or p.shape[0] != p.shape[1] or not np.all(np.isfinite(p)):
            raise ValueError('Projectors must be finite square matrices')
        rank = int(round(float(np.trace(p).real)))
        if not 1 <= rank <= len(p) or not np.isclose(np.trace(p), rank, atol=1e-6, rtol=0):
            raise ValueError('Projectors must have positive integer rank')
        if not np.allclose(p, p.conj().T, atol=1e-6, rtol=0) or not np.allclose(p @ p, p, atol=1e-6, rtol=0):
            raise ValueError('Inputs must be Hermitian support projectors')
        bases.append(projector_basis((p + p.conj().T) / 2, rank=rank))
    if bases[0].shape != bases[1].shape:
        raise ValueError('Projectors must have the same dimension and rank')
    rank = bases[0].shape[1]
    j = (bases[0].shape[0] - 1) / 2
    exact_a, exact_b = [basis @ basis.conj().T for basis in bases]
    best = None

    def retain(angles, reflected):
        nonlocal best
        angles = np.asarray(angles, dtype=float)
        if angles.shape != (3,) or not np.all(np.isfinite(angles)):
            raise ValueError('Each initial alignment must contain three finite angles')
        aligned = rotate_projector_euler(exact_a.conj() if reflected else exact_a, *angles)
        overlap = float(np.trace(aligned @ exact_b).real / rank)
        if best is None or overlap > best['overlap']:
            best = dict(overlap=overlap, alignment_euler=angles.tolist(), reflect=bool(reflected))

    initial_alignments = list(initial_alignments)
    seeds = tuple(seeds)
    if not allow_reflection and any(hint.get('reflect', False) for hint in initial_alignments):
        raise ValueError('Reflected initial alignments require allow_reflection=True')
    for reflected in ((False, True) if allow_reflection else (False,)):
        source = bases[0].conj() if reflected else bases[0]

        def objective(angles):
            rotation = spin_wigner_euler_rotation_matrix(j, *angles)
            cross = bases[1].conj().T @ rotation @ source
            return 1.0 - float(np.vdot(cross, cross).real / rank)

        starts = [[0.0, 0.0, 0.0]] + [hint['alignment_euler'] for hint in initial_alignments
                                     if bool(hint.get('reflect', False)) == reflected]
        for angles in starts:
            retain(angles, reflected)
            result = minimize(objective, angles, method='BFGS', options=dict(gtol=1e-10, maxiter=250))
            retain(result.x, reflected)
        for seed in seeds:
            result = differential_evolution(
                objective, [(0.0, 2*np.pi), (0.0, np.pi), (0.0, 2*np.pi)],
                seed=seed, maxiter=maxiter, popsize=popsize, polish=True,
                tol=1e-9, atol=1e-12, workers=1,
            )
            retain(result.x, reflected)
    return best


def find_projector_alignment_from_files(
    file_a, file_b, *, allow_reflection=True,
    seeds=(0, 1, 2), maxiter=100, popsize=12, initial_alignments=(),
):
    """Find a best-found alignment of the first saved code onto the second.

    Accepts standard nine-field pickle paths with equal photon number and code
    rank. Numerical supports are replaced by exact rank-d projectors, as in
    projector_overlap_from_files, then passed to find_projector_alignment.
    All search options have the same meaning as in that function.

    Returns its dictionary with overlap, alignment_euler (radians), and reflect.
    Apply reflect and then alignment_euler to file_a; file_b stays fixed. These
    arguments use the plot_file_simple convention. Replay the overlap with
    projector_overlap_from_files(..., alignment_euler_a=result['alignment_euler'],
    reflect_a=result['reflect']). No input files are changed, no decoder is
    optimized, and the search does not certify a global maximum.
    """
    projectors, _ = _support_projectors_from_files(file_a, file_b)
    return find_projector_alignment(
        *projectors, allow_reflection=allow_reflection, seeds=seeds,
        maxiter=maxiter, popsize=popsize, initial_alignments=initial_alignments,
    )


# === Group-commutator alignment and representation-block diagnostics ===

def euler_rotation(model, alpha, beta, gamma):
    """Return the spin-n/2 ZYZ Euler rotation for a ``GroupQECModel``."""
    spin = model.spin_sector
    return (
        expm(-1j * alpha * spin.Jz)
        @ expm(-1j * beta * spin.Jy)
        @ expm(-1j * gamma * spin.Jz)
    )


def group_commutator_loss(projector_matrix, representations):
    """Return the normalized average squared group-commutator norm.

    For a rank-d projector P and a unitary representation rho, this is

        1 / (2 d |G|) sum_g ||[P, rho(g)]||_F^2.

    It vanishes exactly when the code projector is group invariant.
    """
    projector_matrix = (
        projector_matrix + projector_matrix.conj().T
    ) / 2
    rank = float(np.real(np.trace(projector_matrix)))
    if rank <= 0:
        raise ValueError('projector_matrix must have positive trace')
    if len(representations) == 0:
        raise ValueError('representations must be nonempty')

    total = 0.0
    for representation in representations:
        commutator_matrix = (
            projector_matrix @ representation
            - representation @ projector_matrix
        )
        total += np.linalg.norm(
            commutator_matrix,
            ord='fro',
        ) ** 2

    return float(total / (2.0 * rank * len(representations)))


def align_projector_by_group_commutator(
    projector_matrix,
    model,
    cyclic=False,
    allow_reflection=True,
    seeds=(0, 1, 2),
    maxiter=100,
    popsize=12,
):
    """Align a code projector by minimizing its group commutators.

    If ``cyclic`` is true, only the two angles specifying the cyclic axis are
    optimized. For nonabelian groups all three ZYZ Euler angles are used.
    If ``allow_reflection`` is true, the optimization is repeated for P*.
    """
    projector_matrix = (
        projector_matrix + projector_matrix.conj().T
    ) / 2
    representations = model.physical_representations

    if cyclic:
        bounds = [
            (0.0, 2.0 * np.pi),
            (0.0, np.pi),
        ]

        def unpack_angles(x):
            return float(x[0]), float(x[1]), 0.0
    else:
        bounds = [
            (0.0, 2.0 * np.pi),
            (0.0, np.pi),
            (0.0, 2.0 * np.pi),
        ]

        def unpack_angles(x):
            return float(x[0]), float(x[1]), float(x[2])

    best = None
    reflection_choices = (
        (False, True) if allow_reflection else (False,)
    )

    for reflected in reflection_choices:
        working_projector = (
            projector_matrix.conj()
            if reflected
            else projector_matrix
        )

        def evaluate(x):
            alpha, beta, gamma = unpack_angles(x)
            rotation = euler_rotation(
                model,
                alpha,
                beta,
                gamma,
            )
            aligned_projector = (
                rotation.conj().T
                @ working_projector
                @ rotation
            )
            loss = group_commutator_loss(
                aligned_projector,
                representations,
            )
            return loss, aligned_projector, rotation

        def objective(x):
            loss, _, _ = evaluate(x)
            return loss

        for seed in seeds:
            result = differential_evolution(
                objective,
                bounds=bounds,
                seed=seed,
                maxiter=maxiter,
                popsize=popsize,
                polish=True,
                updating='immediate',
                workers=1,
            )
            loss, aligned_projector, rotation = evaluate(result.x)

            if best is None or loss < best['commutator_loss']:
                individual_commutators = np.array([
                    np.linalg.norm(
                        aligned_projector @ representation
                        - representation @ aligned_projector,
                        ord='fro',
                    )
                    for representation in representations
                ])
                best = {
                    'commutator_loss': float(loss),
                    'commutator_rms': float(np.sqrt(np.mean(
                        individual_commutators ** 2
                    ))),
                    'commutator_max': float(np.max(
                        individual_commutators
                    )),
                    'angles': unpack_angles(result.x),
                    'reflected': reflected,
                    'aligned_projector': aligned_projector,
                    'rotation': rotation,
                    'optimizer_result': result,
                }

    return best


def representation_block_report(aligned_projector, model):
    """Report isotypic weights and multiplicity-space block structure.

    For each occurring irrep C, the compressed block is compared with
    M_C tensor I_C, where

        M_C = Tr_C(S_C^dagger P S_C) / dim(C).

    For an exactly invariant projector, M_C is an orthogonal projector and
    the reported Schur error vanishes.
    """
    aligned_projector = (
        aligned_projector + aligned_projector.conj().T
    ) / 2
    report = {}

    for irrep_name, multiplicity in model.decomposition.items():
        component = model.isotypic_component(irrep_name)
        embedding = component.embedding_operator
        multiplicity = component.multiplicity
        irrep_dimension = component.dimension

        compressed_matrix = (
            embedding.conj().T
            @ aligned_projector
            @ embedding
        )
        compressed_tensor = compressed_matrix.reshape(
            multiplicity,
            irrep_dimension,
            multiplicity,
            irrep_dimension,
        )

        multiplicity_trace = np.einsum(
            'aibi->ab',
            compressed_tensor,
        )
        multiplicity_trace = (
            multiplicity_trace + multiplicity_trace.conj().T
        ) / 2
        multiplicity_operator = (
            multiplicity_trace / irrep_dimension
        )

        eigenvalues = np.linalg.eigvalsh(multiplicity_operator)
        eigenvalues = np.sort(np.real(eigenvalues))[::-1]

        schur_reconstruction = np.kron(
            multiplicity_operator,
            np.eye(irrep_dimension),
        )
        schur_error = np.linalg.norm(
            compressed_matrix - schur_reconstruction,
            ord='fro',
        )

        report[irrep_name] = {
            'irrep_dimension': irrep_dimension,
            'multiplicity': multiplicity,
            'weight': float(np.real(np.trace(multiplicity_trace))),
            'multiplicity_eigenvalues': eigenvalues,
            'multiplicity_operator': multiplicity_operator,
            'schur_error': float(schur_error),
        }

    return report


def reconstruct_invariant_projector(
    aligned_projector,
    model,
    allocation,
):
    """Reconstruct the nearest invariant projector for a fixed allocation.

    ``allocation`` maps each selected irrep filename to the desired rank of
    its multiplicity-space projector P_C.
    """
    aligned_projector = (
        aligned_projector + aligned_projector.conj().T
    ) / 2
    invariant_projector = np.zeros_like(
        aligned_projector,
        dtype=complex,
    )

    for irrep_name, multiplicity_rank in allocation.items():
        component = model.isotypic_component(irrep_name)
        embedding = component.embedding_operator
        multiplicity = component.multiplicity
        irrep_dimension = component.dimension

        if multiplicity_rank < 1 or multiplicity_rank > multiplicity:
            raise ValueError(
                f'{irrep_name} has multiplicity {multiplicity}, but '
                f'rank {multiplicity_rank} was requested'
            )

        compressed_matrix = (
            embedding.conj().T
            @ aligned_projector
            @ embedding
        )
        compressed_tensor = compressed_matrix.reshape(
            multiplicity,
            irrep_dimension,
            multiplicity,
            irrep_dimension,
        )
        multiplicity_trace = np.einsum(
            'aibi->ab',
            compressed_tensor,
        )
        multiplicity_trace = (
            multiplicity_trace + multiplicity_trace.conj().T
        ) / 2

        eigenvalues, eigenvectors = np.linalg.eigh(
            multiplicity_trace
        )
        selected = np.argsort(eigenvalues)[::-1][
            :multiplicity_rank
        ]
        selected_vectors = eigenvectors[:, selected]
        multiplicity_projector = (
            selected_vectors @ selected_vectors.conj().T
        )

        invariant_projector += (
            embedding
            @ np.kron(
                multiplicity_projector,
                np.eye(irrep_dimension),
            )
            @ embedding.conj().T
        )

    return (
        invariant_projector + invariant_projector.conj().T
    ) / 2


def analyze_group_representation_from_file(
    result_file,
    group,
    seeds=(0, 1, 2),
    maxiter=100,
    popsize=12,
    allow_reflection=True,
    verbose=True,
):
    """Run the complete group-representation analysis for one result file.

    The function loads the numerical code projector, aligns it by minimizing
    its group-commutator loss, reports every isotypic/multiplicity block,
    selects the rank-compatible representation content with the largest
    multiplicity-eigenvalue weight, and reconstructs the corresponding exact
    invariant projector.

    Only ``result_file`` and ``group`` are required. A group folder whose
    basename starts with ``C`` is treated as cyclic, so only its symmetry-axis
    orientation is optimized.
    """
    numerical_projector, metadata = spin_projector_from_result_file(
        result_file
    )

    # Replace the slightly imperfect support matrix by an exact rank-d
    # projector before evaluating group commutators.
    numerical_basis = projector_basis(
        (numerical_projector + numerical_projector.conj().T) / 2,
        rank=metadata['d'],
    )
    numerical_projector = (
        numerical_basis @ numerical_basis.conj().T
    )

    group_model = GroupQECModel(group, n=metadata['n'])
    group_name = os.path.basename(os.path.normpath(group))
    cyclic = group_name.startswith('C')

    alignment = align_projector_by_group_commutator(
        numerical_projector,
        group_model,
        cyclic=cyclic,
        allow_reflection=allow_reflection,
        seeds=seeds,
        maxiter=maxiter,
        popsize=popsize,
    )

    block_report = representation_block_report(
        alignment['aligned_projector'],
        group_model,
    )

    # Dynamic programming over the allowed ranks rank(P_C). Choosing one
    # multiplicity direction in C costs dim(C) physical code dimensions and
    # contributes dim(C) times its normalized multiplicity eigenvalue to
    # Tr(P Q). The selected allocation must have total rank d.
    target_rank = metadata['d']
    states = {0: (0.0, {})}

    for irrep_name, block in block_report.items():
        irrep_dimension = block['irrep_dimension']
        multiplicity = block['multiplicity']
        eigenvalues = block['multiplicity_eigenvalues']
        new_states = dict(states)

        for used_rank, (score, allocation) in states.items():
            max_multiplicity_rank = min(
                multiplicity,
                (target_rank - used_rank) // irrep_dimension,
            )
            for multiplicity_rank in range(
                1,
                max_multiplicity_rank + 1,
            ):
                physical_rank = (
                    used_rank
                    + multiplicity_rank * irrep_dimension
                )
                contribution = (
                    irrep_dimension
                    * np.sum(eigenvalues[:multiplicity_rank])
                )
                candidate_score = score + contribution
                previous = new_states.get(physical_rank)

                if previous is None or candidate_score > previous[0]:
                    candidate_allocation = dict(allocation)
                    candidate_allocation[irrep_name] = multiplicity_rank
                    new_states[physical_rank] = (
                        float(candidate_score),
                        candidate_allocation,
                    )

        states = new_states

    if target_rank not in states:
        raise RuntimeError(
            f'No rank-{target_rank} invariant representation projector '
            f'can be assembled for {group_name} at n={metadata["n"]}'
        )

    _, allocation = states[target_rank]
    invariant_projector = reconstruct_invariant_projector(
        alignment['aligned_projector'],
        group_model,
        allocation,
    )
    invariant_overlap = float(np.real(
        np.trace(
            alignment['aligned_projector']
            @ invariant_projector
        ) / target_rank
    ))
    invariant_commutator_loss = group_commutator_loss(
        invariant_projector,
        group_model.physical_representations,
    )
    invariant_idempotence_error = float(np.linalg.norm(
        invariant_projector @ invariant_projector
        - invariant_projector,
        ord='fro',
    ))

    if verbose:
        print('Group-representation analysis')
        print('  file:', result_file)
        print('  group:', group_name)
        print('  normalized commutator loss:', alignment['commutator_loss'])
        print('  maximum commutator:', alignment['commutator_max'])
        print('  Euler angles:', alignment['angles'])
        print('  reflected:', alignment['reflected'])
        print('Representation blocks')
        for irrep_name, block in block_report.items():
            print(' ', irrep_name)
            print('    weight:', block['weight'])
            print(
                '    multiplicity eigenvalues:',
                block['multiplicity_eigenvalues'],
            )
            print('    Schur error:', block['schur_error'])
        print('Selected representation allocation:', allocation)
        print('Reconstructed invariant-code overlap:', invariant_overlap)
        print(
            'Reconstructed commutator loss:',
            invariant_commutator_loss,
        )
        print(
            'Reconstructed idempotence error:',
            invariant_idempotence_error,
        )

    return {
        'result_file': result_file,
        'group': group_name,
        'metadata': metadata,
        'model': group_model,
        'numerical_projector': numerical_projector,
        'alignment': alignment,
        'block_report': block_report,
        'allocation': allocation,
        'invariant_projector': invariant_projector,
        'invariant_overlap': invariant_overlap,
        'invariant_commutator_loss': invariant_commutator_loss,
        'invariant_idempotence_error': invariant_idempotence_error,
    }


if __name__ == '__main__':
    # Example: perform the complete analysis using only a result file and a
    # group folder. The rank-two representation allocation is inferred from
    # the aligned multiplicity blocks.
    
    '''analysis = analyze_group_representation_from_file(
        'repro_sweep/n4/n4_eta0.90_d2_r21000_s1.pkl',
        '2O',
    )'''

    plot_file_simple('n7/n7_eta0.95_d2_r21000_s0.pkl', reflect=False, rotate= False)
    #plot_literature_code(25, 'OG', '2O', eta=0.9)
    #plot_literature_code(22, 'AAB', eta=0.9)
      
#alpha=4.084131192108156, beta=0.8240296164093108, gamma=2.562622339313475 for n = 13
