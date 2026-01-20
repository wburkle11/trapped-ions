import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ------------------- initial thermal -------------------
nbar0 = 10.0
n_max = 90
n = np.arange(n_max+1)
p0 = (nbar0**n) / ((nbar0+1.0)**(n+1.0))
p0 /= p0.sum()

# -------- perfect first-RSB π pulse (no recoil) ----------
def perfect_rsb_step(p):
    """
    Perfect Δn=1 π-pulse:
      all population in n>=1 transfers to n-1; n=0 stays.
      Equivalent to: p'(0)=p(0)+p(1); p'(i)=p(i+1) for i>=1; last bin -> 0.
    """
    p = np.asarray(p, float)
    out = np.zeros_like(p)
    if p.size == 1:
        out[0] = 1.0
        return out
    out[0]      = p[0] + p[1]
    out[1:-1]   = p[2:]
    out[-1]     = 0.0
    s = out.sum()
    return out/s if s>0 else out

# (Optional) perfect Δn=2 pulse you can sprinkle in if you want:
def perfect_rsb_step_d(p, d=1):
    p = np.asarray(p, float)
    out = np.zeros_like(p)
    out[0] = p[:d+1].sum()                # n=0 keeps its pop + inflow from 1..d
    if d+1 < len(p):
        out[1:len(p)-d] = p[d+1:]         # shift left by d
    s = out.sum()
    return out/s if s>0 else out

# ------------------- simulate -------------------
n_steps = 60
states = [p0]
for k in range(n_steps):
    # only perfect Δn=1 pulses:
    states.append(perfect_rsb_step(states[-1]))
    # If you want a perfect Δn=2 every so often, replace with:
    # if (k+1) % 5 == 0:
    #     states.append(perfect_rsb_step_d(states[-1], d=2))
    # else:
    #     states.append(perfect_rsb_step(states[-1]))

# ------------------- animation -------------------
ymax = min(1.05, 1.05*max(f.max() for f in states))
fig, ax = plt.subplots(figsize=(10,4))
bars = ax.bar(n, states[0], width=0.9, edgecolor='black', linewidth=0.3)
txt  = ax.text(0.01, 0.97, "", transform=ax.transAxes, va='top', ha='left', family='monospace')

ax.set_xlim(-0.5, min(60, int(max(8*nbar0+40, 2*nbar0+60))))
#ax.set_ylim(0, ymax)
ax.set_xlabel("Motional state $n$")
ax.set_ylabel("Probability")
ax.grid(True, linestyle='--', alpha=0.35)
ax.set_title("Perfect RSB π-pulses")

def update(k):
    pk = states[k]
    nbar = np.dot(n, pk)

    for rect, h in zip(bars, pk):
        rect.set_height(h)

    txt.set_text(
    fr"pulse: {k:3d}  |  $\bar{{n}}_{{\mathrm{{dop}}}}$ = {nbar0:.2f}  |  $\langle n\rangle$ = {nbar:6.3f}  |  $P_0$ = {pk[0]:6.3f}"
    )

    return bars, txt
anim = FuncAnimation(fig, update, frames=len(states), interval=140, blit=False)
plt.tight_layout()


anim.save("Z:\\Users\\Wes\\Code\\Cooling\\rsc_adaptive_no_recoil.gif", writer="pillow", dpi=120)  # or show:

