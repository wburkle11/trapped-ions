
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.optimize import root
from scipy import constants as const


# ----------------------------
# Physics: equilibrium + modes
# ----------------------------

def equilibrium_positions_1d(N, max_iter=2000):
    x0 = np.linspace(-(N - 1) / 2, (N - 1) / 2, N)

    def F(u):
        f = np.zeros_like(u)
        for i in range(N):
            s = 0.0
            for j in range(N):
                if i == j:
                    continue
                rij = u[i] - u[j]
                s += rij / (abs(rij) ** 3)
            f[i] = u[i] - s
        return f

    sol = root(F, x0, method="hybr", options={"maxfev": max_iter})
    if not sol.success:
        raise RuntimeError("Equilibrium solver failed: {}".format(sol.message))
    u = np.array(sol.x, dtype=float)
    u.sort()
    return u


def transverse_modes(N, wz, wx):
    """
    Returns omega_m (rad/s) and eigenvectors b (N,N) with columns b[:,m].
    """
    u = equilibrium_positions_1d(N)

    C = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            C[i, j] = 1.0 / (abs(u[i] - u[j]) ** 3)
    D = np.diag(C.sum(axis=1))

    r = (wx / wz) ** 2
    A = r * np.eye(N) - D + C

    lam, b = np.linalg.eigh(A)
    omega_m = wz * np.sqrt(lam)
    return omega_m, b


# ----------------------------
# Lamb–Dicke + MS couplings
# ----------------------------

def delta_k_from_raman_geometry(lambda_m, crossing_angle_deg):
    k = 2.0 * np.pi / lambda_m
    theta = np.deg2rad(crossing_angle_deg)
    return 2.0 * k * np.sin(theta / 2.0)


def compute_eta_im(omega_m, b, delta_k, mass_kg):
    x0_m = np.sqrt(const.hbar / (2.0 * mass_kg * omega_m))  # (M,) meters
    eta_im = delta_k * (b * x0_m[None, :])                  # (N,M)
    return eta_im


def build_g_im_from_eta(eta_im, Omega_hz):
    Omega = 2.0 * np.pi * Omega_hz
    return 0.5 * eta_im * Omega


def magnus_chi_and_alpha(omega_m, g_im, omega_com, detuning_hz, t_gate):
    mu = omega_com + 2.0 * np.pi * detuning_hz
    delta = mu - omega_m

    eps_hz = 1.0
    eps = 2.0 * np.pi * eps_hz

    # If delta is too small, replace it by ±eps with the same sign as delta.
    # (If delta is exactly 0, choose +eps by convention.)
    sgn = np.where(delta >= 0.0, 1.0, -1.0)
    delta_safe = np.where(np.abs(delta) < eps, sgn * eps, delta)

    t = float(t_gate)
    F = (t / delta_safe) - (np.sin(delta_safe * t) / (delta_safe ** 2))

    chi = (2.0 * (g_im * F[None, :])) @ g_im.T
    chi = 0.5 * (chi + chi.T)
    np.fill_diagonal(chi, 0.0)

    alpha = (g_im / delta_safe[None, :]) * (1.0 - np.exp(1j * delta_safe[None, :] * t))
    return chi, alpha, delta


def effective_J_from_chi(chi, t_gate):
    return chi / float(t_gate)


# ----------------------------
# Plotting helpers (improved)
# ----------------------------

def plot_mode_spectrum_auto(omega_m, detunings_kHz, annotate=True):
    """
    Auto-scaled mode spectrum plot:
      - plots (ω_m - ω_COM)/2π in kHz
      - overlays detunings (kHz) as red sticks
      - auto x-limits/ticks based on data range
    """
    omega_com = omega_m[np.argmax(omega_m)]
    f_rel_kHz = (omega_m - omega_com) / (2.0 * np.pi * 1e3)

    # gather all x-values we want in view
    det_vals = np.array(list(detunings_kHz.values()), dtype=float) if detunings_kHz else np.array([])
    all_x = np.concatenate([f_rel_kHz, det_vals]) if det_vals.size else f_rel_kHz.copy()

    xmin, xmax = float(np.min(all_x)), float(np.max(all_x))
    span = max(1.0, xmax - xmin)
    pad = 0.08 * span
    xmin -= pad
    xmax += pad

    # choose tick step nicely
    target_ticks = 8
    raw_step = span / target_ticks
    nice_steps = np.array([10, 20, 50, 100, 200, 500, 1000, 2000, 5000], dtype=float)
    step = nice_steps[np.argmin(np.abs(nice_steps - raw_step))]
    tick_min = step * np.floor(xmin / step)
    tick_max = step * np.ceil(xmax / step)
    ticks = np.arange(tick_min, tick_max + 0.5 * step, step)

    fig, ax = plt.subplots(figsize=(9.0, 2.2))
    ax.axhline(0, color="0.6", linewidth=0.8)

    # mode spectrum (black)
    for x in f_rel_kHz:
        ax.vlines(x, 0, 1.0, colors="k", linewidth=2.0)

    # detunings (red)
    if detunings_kHz:
        for label, dk in detunings_kHz.items():
            ax.vlines(dk, 0, 0.75, colors="red", linewidth=2.0)
            if annotate:
                ax.text(dk + 0.02 * span, 0.78, str(label), fontsize=10)

    ax.set_ylim(-0.05, 1.05)
    ax.set_yticks([])
    ax.set_xlabel("Mode frequency - COM frequency (kHz)")
    ax.set_xlim(xmin, xmax)

    ax.set_xticks(ticks)
    # label 0 tick as "COM" if it's very close to 0
    ticklabels = []
    for t in ticks:
        if abs(t) < 1e-9:
            ticklabels.append("COM (0)")
        else:
            ticklabels.append("{:.0f}".format(t))
    ax.set_xticklabels(ticklabels)

    for spine in ["left", "right", "top"]:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    return fig


def bar3d_J(ax, J_kHz_over_2pi, title, zlim=None):
    N = J_kHz_over_2pi.shape[0]
    xs, ys = np.meshgrid(np.arange(1, N + 1), np.arange(1, N + 1), indexing="ij")
    x = xs.ravel()
    y = ys.ravel()
    z0 = np.zeros_like(x, dtype=float)
    dz = J_kHz_over_2pi.ravel()

    dx = 0.7 * np.ones_like(dz)
    dy = 0.7 * np.ones_like(dz)

    vmax = np.max(np.abs(dz)) + 1e-12
    colors = plt.cm.Spectral((dz + vmax) / (2 * vmax))

    ax.bar3d(x - 0.35, y - 0.35, z0, dx, dy, dz,
             shade=True, color=colors, linewidth=0.0)

    ax.set_title(title, pad=4, fontsize=10)
    ax.set_xlabel(r"$j$", labelpad=-6)
    ax.set_ylabel(r"$i$", labelpad=-6)
    ax.set_zlabel(r"$J_{i,j}\;(\mathrm{kHz})/2\pi$", labelpad=6)

    ax.set_xlim(0.5, N + 0.5)
    ax.set_ylim(0.5, N + 0.5)
    # Force symmetric z-axis so negative values are visible
    if zlim is None:
        zmax = np.max(np.abs(dz)) + 1e-12
        ax.set_zlim(-zmax, zmax)
    else:
        ax.set_zlim(zlim[0], zlim[1])

    ax.plot_surface(
        xs, ys,
        np.zeros_like(xs),
        color='k',
        alpha=0.05
        )
    ax.view_init(elev=25, azim=-55)


def plot_single_coupling_panel(omega_m, g_im, omega_com, detuning_kHz, t_gate,
                               title_prefix="COM", zlim=None):
    """
    Makes a single 3D coupling plot for one detuning choice.
    Returns (fig, residual_metric, J_matrix)
    """
    det_hz = float(detuning_kHz) * 1e3
    chi, alpha, delta = magnus_chi_and_alpha(
        omega_m=omega_m,
        g_im=g_im,
        omega_com=omega_com,
        detuning_hz=det_hz,
        t_gate=t_gate
    )
    J = effective_J_from_chi(chi, t_gate)
    J_kHz_over_2pi = (J / (2.0 * np.pi)) / 1e3

    fig = plt.figure(figsize=(5.2, 4.6))
    ax = fig.add_subplot(111, projection="3d")
    title = "{} {:+.0f} kHz".format(title_prefix, detuning_kHz)
    bar3d_J(ax, J_kHz_over_2pi, title=title, zlim=zlim)
    fig.tight_layout()

    residual = float(np.max(np.abs(alpha)))
    return fig, residual, J

def plot_normalized_J_heatmap(J, title="Normalized J_ijs", cmap="viridis"):
    """
    Plot a heatmap of normalized |J_ij| in [0,1].
    - Uses off-diagonal entries to set the normalization scale.
    - Diagonal is set to 0 for display.
    """
    J = np.array(J, dtype=float)
    N = J.shape[0]

    # Use magnitude (often what you want for a "strength" heatmap)
    A = np.abs(J).copy()
    np.fill_diagonal(A, 0.0)

    # Normalize by max off-diagonal
    off = A[~np.eye(N, dtype=bool)]
    max_off = float(np.max(off)) if off.size else 0.0

    if max_off > 0:
        A_norm = A / max_off
    else:
        A_norm = A  # all zeros

    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    im = ax.imshow(A_norm, origin="upper", interpolation="nearest", vmin=0.0, vmax=1.0, cmap=cmap)
    ax.set_title(title)
    ax.set_xlabel("j")
    ax.set_ylabel("i")

    # match your example: ticks 0..N-1
    ax.set_xticks(np.arange(N))
    ax.set_yticks(np.arange(N))

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Normalized |J_ij|")

    fig.tight_layout()
    return fig

def distance_profile_J(J, use_abs=True, exclude_edges=0):
    """
    Build an averaged coupling vs distance r = |i-j|.

    exclude_edges: if >0, only use ions i,j in [exclude_edges, N-1-exclude_edges]
                   to reduce edge effects.
    Returns:
      r_vals: array [1..N-1]
      J_mean: mean coupling magnitude at each r
      J_std:  std at each r (spread due to non-translation invariance)
    """
    J = np.array(J, dtype=float)
    N = J.shape[0]
    if use_abs:
        M = np.abs(J)
    else:
        M = J.copy()

    i0 = int(exclude_edges)
    i1 = N - int(exclude_edges)

    rs = np.arange(1, N)  # distances
    J_mean = np.zeros_like(rs, dtype=float)
    J_std = np.zeros_like(rs, dtype=float)

    for k, r in enumerate(rs):
        vals = []
        for i in range(i0, i1):
            j = i + r
            if j < i1:
                vals.append(M[i, j])
        vals = np.array(vals, dtype=float)
        if vals.size > 0:
            J_mean[k] = np.mean(vals)
            J_std[k] = np.std(vals)
        else:
            J_mean[k] = np.nan
            J_std[k] = np.nan

    return rs, J_mean, J_std


def fit_power_law(r_vals, J_vals, r_min=1, r_max=None):
    """
    Fit J(r) ≈ J0 / r^alpha using a log-log linear regression.
    Returns (J0, alpha).
    """
    r_vals = np.array(r_vals, dtype=float)
    J_vals = np.array(J_vals, dtype=float)

    if r_max is None:
        r_max = np.nanmax(r_vals)

    mask = (r_vals >= r_min) & (r_vals <= r_max) & np.isfinite(J_vals) & (J_vals > 0)
    x = np.log(r_vals[mask])
    y = np.log(J_vals[mask])

    if x.size < 2:
        raise ValueError("Not enough points to fit power law. Adjust r_min/r_max.")

    # y = log(J0) - alpha * log(r)
    slope, intercept = np.polyfit(x, y, 1)
    alpha = -slope
    J0 = np.exp(intercept)
    return J0, alpha

def per_mode_alpha_closure(omega_m, g_im, omega_com, detuning_hz, t_gate, eps_hz=1.0):
    """
    Compute per-mode closure metrics at time t_gate.

    Returns:
      delta_rad_s : (M,) detunings in rad/s
      alpha_im    : (N,M) complex displacements at t_gate
      max_abs     : (M,) max_i |alpha_im|
      rms_abs     : (M,) sqrt(mean_i |alpha_im|^2)
    """
    mu = omega_com + 2.0 * np.pi * detuning_hz
    delta = mu - omega_m  # (M,) rad/s

    # Safe floor to avoid divide-by-zero; if you're truly at resonance,
    # loop closure isn't well-defined in the effective-Ising picture anyway.
    eps = 2.0 * np.pi * eps_hz
    delta_safe = np.where(np.abs(delta) < eps, np.sign(delta) * eps + eps, delta)

    t = float(t_gate)
    alpha_im = (g_im / delta_safe[None, :]) * (1.0 - np.exp(1j * delta_safe[None, :] * t))

    abs_alpha = np.abs(alpha_im)
    max_abs = np.max(abs_alpha, axis=0)
    rms_abs = np.sqrt(np.mean(abs_alpha**2, axis=0))
    return delta, alpha_im, max_abs, rms_abs


def print_loop_closure_report(omega_m, g_im, omega_com, detuning_kHz, t_gate,
                              threshold=0.10, top_k=10):
    """
    Print whether each mode's phase-space loop is 'adequately closed' at t_gate.

    threshold: closure criterion on max_i |alpha_im(T)| (dimensionless).
               Typical: 0.05 (strict), 0.1 (reasonable), 0.2 (loose)
    top_k: print the worst top_k modes sorted by max |alpha|.
    """
    detuning_hz = float(detuning_kHz) * 1e3
    delta, alpha_im, max_abs, rms_abs = per_mode_alpha_closure(
        omega_m=omega_m,
        g_im=g_im,
        omega_com=omega_com,
        detuning_hz=detuning_hz,
        t_gate=t_gate,
        eps_hz=1.0
    )

    # Convert to more readable units
    delta_hz = delta / (2.0 * np.pi)

    # Sort modes by worst closure
    order = np.argsort(-max_abs)

    print("\n--- Loop-closure report ---")
    print(f"Detuning relative to COM: {detuning_kHz:+.1f} kHz")
    print(f"Gate time T = {t_gate*1e6:.3f} µs")
    print(f"Closure threshold on max_i |alpha_im(T)|: {threshold:.3f}")
    print("")
    print("Mode  (ω_m/2π MHz)   δ_m/2π (kHz)   max_i|α_im|    rms_i|α_im|   CLOSED?")

    for rank, m in enumerate(order[:min(top_k, len(order))]):
        om_mhz = omega_m[m] / (2*np.pi) / 1e6
        dm_khz = delta_hz[m] / 1e3
        closed = (max_abs[m] < threshold)
        flag = "YES" if closed else "NO"
        print(f"{m:>3d}   {om_mhz:>10.6f}   {dm_khz:>12.3f}   {max_abs[m]:>10.3e}   {rms_abs[m]:>10.3e}   {flag}")

    # Also print a quick global summary
    n_closed = int(np.sum(max_abs < threshold))
    print(f"\nClosed modes: {n_closed}/{len(omega_m)} (by max_i |alpha| < {threshold})")

    # Return arrays in case you want to plot them
    return delta, max_abs, rms_abs



# ----------------------------
# Main (example)
# ----------------------------

def main():
    # ---- system parameters ----
    N = 7
    wz = 2.0 * np.pi * 0.2e6
    wx = 2.0 * np.pi * 1.2e6

    # ---- physical constants ----
    mass_yb171 = 171.0 * const.atomic_mass

    # ---- Raman Δk ----
    lambda_raman = 355e-9
    crossing_angle_deg = 90.0
    delta_k = delta_k_from_raman_geometry(lambda_raman, crossing_angle_deg)

    # ---- drive ----
    Omega_hz = 50e3

    # ---- gate time ----
    t_gate = 200e-6

    # Detunings relative to COM (kHz) (still used for red markers in spectrum plot)
    detunings_kHz = {
        "MS-detuning": 10,
    }

    # choose ONE detuning to plot J_ij for (pick key from detunings_kHz)
    detuning_key = "MS-detuning"     # e.g. "b"
    detuning_kHz = detunings_kHz[detuning_key]

    # ---- compute modes ----
    omega_m, b = transverse_modes(N, wz=wz, wx=wx)
    omega_com = omega_m[np.argmax(omega_m)]

    # ---- eta and g ----
    eta_im = compute_eta_im(omega_m, b, delta_k=delta_k, mass_kg=mass_yb171)
    g_im = build_g_im_from_eta(eta_im, Omega_hz=Omega_hz)
    
    # Print loop-closure diagnostics per mode
    print_loop_closure_report(
        omega_m=omega_m,
        g_im=g_im,
        omega_com=omega_com,
        detuning_kHz=detuning_kHz,
        t_gate=t_gate,
        threshold=0.1,   # try 0.05 for strict, 0.2 for loose
        top_k=N           # print all modes if you want
        )

    # ----- improved mode spectrum plot (auto-scaled) -----
    plot_mode_spectrum_auto(omega_m, detunings_kHz, annotate=True)

    # ----- single coupling plot -----
    figJ, residual, J = plot_single_coupling_panel(
        omega_m=omega_m,
        g_im=g_im,
        omega_com=omega_com,
        detuning_kHz=detuning_kHz,
        t_gate=t_gate,
        title_prefix="COM",
        zlim=None  # or set e.g. (-2,2) to fix scale
    )
    
    plot_normalized_J_heatmap(J, title="Normalized J_ijs")
    
    # --- power-law extraction ---
    r, Jmean, Jstd = distance_profile_J(J, use_abs=True, exclude_edges=0)

    # Choose fit window (often ignore r=1 if very non-power-law; up to you)
    J0, alpha = fit_power_law(r, Jmean, r_min=1, r_max=None)

    print(f"\nPower-law fit |J(r)| ≈ J0 / r^alpha")
    print(f"  J0 = {J0:.3e} (same units as J)")
    print(f"  alpha = {alpha:.3f}")

    # Plot J(r)
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.errorbar(r, Jmean, yerr=Jstd, fmt="o", capsize=3)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("r = |i-j|")
    ax.set_ylabel("mean |J_ij| at r")
    ax.set_title(f"Power-law fit: alpha={alpha:.2f}")
    fig.tight_layout()

    #print("Residual displacement metric: max |alpha_im(t_gate)| = {:.3e}".format(residual))
    plt.show()


if __name__ == "__main__":
    main()


