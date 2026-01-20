import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.constants as scipy

# ---------- constants & params ----------
q_e = scipy.elementary_charge
eps_o = scipy.epsilon_0
M_Yb = 171*1.67262192e-27
delta_k = np.sqrt(2)*2*np.pi/(355e-9)
hbar = scipy.hbar

fz = .7   # MHz   axial
fx = 2.0  # MHz   radial
omega_z = 2*np.pi*fz*1e6
omega_x = 2*np.pi*fx*1e6

N = 5

# ---------- helpers from your code ----------
def equil_positions(N):
    def coupled_equations(u):
        eqs = []
        for m in range(N):
            sum1 = sum(1/(u[m]-u[n])**2 for n in range(m))
            sum2 = sum(1/(u[m]-u[n])**2 for n in range(m+1, N))
            eqs.append(u[m] - sum1 + sum2)
        return eqs
    initial_guess = [m/10 for m in range(1, N+1)]
    sol = fsolve(coupled_equations, initial_guess)
    return np.round(sol, 4)

def eigensystem(frad, fax, N):
    u = equil_positions(N)
    Anm = np.empty((N, N), dtype=float)
    for n in range(N):
        for m in range(N):
            if n == m:
                Anm[n, n] = (frad/fax)**2 - np.sum(
                    [1/abs(u[m]-u[p])**3 for p in range(N) if p != m]
                )
            else:
                Anm[n, m] = 1/abs(u[m]-u[n])**3
    vals, vecs = np.linalg.eig(Anm)
    vals = np.sqrt(vals)*fz  # MHz
    vecs = vecs.T            # modes × ions
    return vals, vecs

# ---------- compute spectrum & mode shapes ----------
evals_MHz, evects = eigensystem(fx, fz, N)          # evals in MHz
idx = np.argsort(evals_MHz)                          # sort by freq
freqs_MHz = np.asarray(evals_MHz)[idx]
modes = np.asarray(evects)[idx]                      # shape: (N modes, N ions)

# normalize each mode so max |amplitude| = 1 (just for plotting)
modes_norm = modes / np.max(np.abs(modes), axis=1, keepdims=True)

# identify COM (highest frequency)
COM_col = np.argmax(freqs_MHz)

# ---------- plotting ----------
fig = plt.figure(figsize=(10, 6))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 2], hspace=0.35)

# (1) stick plot of frequencies
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(freqs_MHz, np.zeros_like(freqs_MHz), 'ko', ms=5)
ax1.axvline(freqs_MHz[-1], color='tab:blue', ls='--', lw=1,
            label=f"COM ≈ {freqs_MHz[-1]:.2f} MHz")
ax1.set_yticks([])
ax1.set_xlabel("Frequency (MHz)")
ax1.set_title("1D Chain Normal-Mode Frequencies")
ax1.legend(loc='upper left', fontsize=9)

# bandwidth annotation
bandwidth = freqs_MHz[-1] - freqs_MHz[0]
ax1.annotate(f'Bandwidth: {bandwidth:.2f} MHz',
             xy=(0.5*(freqs_MHz[-1]+freqs_MHz[0]), 0.02),
             xycoords=('data','axes fraction'),
             ha='center', va='bottom', fontsize=9, color='tab:blue')

# (2) heatmap of signed mode amplitudes (eigenvectors)
#ax2 = fig.add_subplot(gs[1, 0])
#im = ax2.imshow(modes_norm.T, aspect='auto', cmap='coolwarm',
                #vmin=-1, vmax=1, origin='lower')
# axes: rows = ion index, columns = mode (sorted by freq)
#ax2.set_ylabel("Ion index")
#ax2.set_xlabel("Mode (sorted by frequency)")

# put frequency labels under each column
#ax2.set_xticks(np.arange(N))
#ax2.set_xticklabels([f"{f:.2f}" for f in freqs_MHz], rotation=45, ha='right')
#ax2.set_title("Normal-Mode Amplitudes (signed eigenvectors)")

# colorbar
#cbar = fig.colorbar(im, ax=ax2, pad=0.02)
#cbar.set_label("Amplitude (relative)")

plt.tight_layout()
plt.show()