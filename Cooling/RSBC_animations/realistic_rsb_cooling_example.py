import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.special import eval_genlaguerre, gammaln

# ------------------- knobs -------------------
nbar0   = 4.0
n_max   = 100
eta     = 0.08
Omega0  = 2*np.pi*250e3         # carrier-scale Rabi [rad/s]
n_steps = 230                   # more cycles with adaptive τ
alpha   = 0.00                  # recoil (0 => no recoil)
use_D2  = True                  # also apply 2nd RSB occasionally
D2_every= 0                     # do a Δn=2 pulse every M cycles

# ------------------- initial thermal -------------------
n  = np.arange(n_max+1)
p0 = (nbar0**n) / ((nbar0+1.0)**(n+1.0))
p0 /= p0.sum()

# ----- exact Laguerre Rabi rate Ω_{n,n-Δ} -----
def rabi_rate_delta(n, d, Omega0, eta):
    n = np.asarray(n, dtype=int)
    np_ = n - d
    valid = (np_ >= 0)
    out = np.zeros_like(n, dtype=float)
    if not np.any(valid): 
        return out
    n_less  = np.minimum(n[valid], np_[valid])
    n_great = np.maximum(n[valid], np_[valid])
    dn      = np.abs(d)
    logfac  = 0.5*(gammaln(n_great+1) - gammaln(n_less+1))
    sqrt_ratio = np.exp(logfac)
    L = eval_genlaguerre(n_less, dn, eta**2)
    out[valid] = Omega0 * np.exp(-eta**2/2) * (eta**dn) * sqrt_ratio * np.abs(L)
    return out

def rsb_step_general(p, d=1, tau=0.0):
    """Apply a red-sideband pulse n->n-d with duration tau (incoherent sin^2 map)."""
    out = p.copy()
    idx = np.arange(out.size)
    valid = idx >= d
    Om = rabi_rate_delta(idx, d, Omega0, eta)
    f  = np.sin(0.5*Om*tau)**2
    f[~valid] = 0.0
    moved = f * out
    out -= moved
    # add to (n-d); shift by d
    out[:-d] += moved[d:]
    return out

# def recoil_step(p, alpha):
#     if alpha <= 0: 
#         return p
#     out = p.copy()
#     up  = alpha * out
#     out -= up
#     out[1:] += up[:-1]
#     out[-1] += up[-1]
#     return out

def adaptive_tau(p, d=1):
    """Choose τ for a π-pulse at the modal n (max-probability bin) for Δn=d."""
    k = np.argmax(p)
    k = max(k, d)              # ensure allowed transition
    Om = rabi_rate_delta(np.array([k]), d, Omega0, eta)[0]
    if Om == 0:                # guard against a Laguerre node
        # nudge to nearest neighbor with nonzero coupling
        for dk in (1, -1, 2, -2):
            jj = k+dk
            if d <= jj < len(p):
                Om = rabi_rate_delta(np.array([jj]), d, Omega0, eta)[0]
                if Om > 0: break
    return np.pi/Om if Om>0 else 0.0

# ------------- simulate with adaptive pulses -------------
states = [p0]
for step in range(1, n_steps+1):
    p = states[-1]

    # 1st RSB with π on current mode
    tau1 = adaptive_tau(p, d=1)
    p = rsb_step_general(p, d=1, tau=tau1)

    # occasional 2nd RSB to clear even/odd pockets (optional)
    #if use_D2 and (step % D2_every == 0):
        #tau2 = adaptive_tau(p, d=2)
        #p = rsb_step_general(p, d=2, tau=tau2)

    # optical pumping recoil (set alpha=0 for none)
    #p = recoil_step(p, alpha)

    # normalize & store
    s = p.sum()
    states.append(p/s if s>0 else p)

# ------------------- animation -------------------
ymax = min(1.05, 1.05*max(frame.max() for frame in states))
fig, ax = plt.subplots(figsize=(10,4))
bars = ax.bar(n, states[0], width=0.9, edgecolor='black', linewidth=0.3)
txt  = ax.text(0.01, 0.97, "", transform=ax.transAxes, va='top', ha='left', family='monospace')

ax.set_xlim(-0.5, min(n_max, int(max(8*nbar0+40, 2*nbar0+60))))
ax.set_ylim(0, ymax)
ax.set_xlabel("Motional state $n$")
ax.set_ylabel("Probability")
ax.grid(True, linestyle='--', alpha=0.35)
ax.set_title("Resolved sideband cooling with Laguerre rates")

def update(k):
    pk = states[k]
    for rect, h in zip(bars, pk):
        rect.set_height(h)
    txt.set_text(f"cycle: {k:3d}  |  ⟨n⟩ = {np.dot(n, pk):6.3f}")
    return bars, txt

anim = FuncAnimation(fig, update, frames=len(states), interval=140, blit=False)
plt.tight_layout()
anim.save("Z:\\Users\\Wes\\Code\\Cooling\\rsc_adaptive_no_recoil.gif", writer="pillow", dpi=120)  # or show:

