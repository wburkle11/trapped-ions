import numpy as np, qutip as qt, matplotlib.pyplot as plt
from matplotlib.patches import Patch

MHz = 2*np.pi*1e6
Ωσm, Ωπ = 25*MHz, 1*MHz
Δ       = 50*MHz

g_m, g_p = Ωσm/2, Ωπ/2
H = qt.Qobj(np.array([[0.0,  g_m,  -g_p],
                      [g_m,  Δ,     0.0],
                      [-g_p, 0.0,   Δ + 5*MHz ]], dtype=complex))

E, V = H.eigenstates()
E = np.array([float(ev) for ev in E])
bare = [qt.basis(3,0), qt.basis(3,1), qt.basis(3,2)]  # |e>,|+>,|0>
W = np.array([[abs((b.dag()*v)[0,0])**2 for v in V] for b in bare])

o = np.argsort(E); E, W = E[o], W[:,o]
x = E/(2*np.pi*1e6)

plt.rcParams.update({"font.size": 11})
fig, ax = plt.subplots(figsize=(8.6, 2.6), dpi=150)

gaps = np.diff(np.sort(x)); gaps = gaps[gaps>0]
cluster = 0.35*(gaps.min() if gaps.size else (x.ptp() if x.ptp()>0 else 1.0))
bar_w, step = 0.35*cluster, 0.40*cluster
offsets = np.array([-step, 0.0, step])

colors = {'|0⟩':'#f49c20','|+⟩':'#2ca02c','|e⟩':'#666666'}
rows, cols = [2,1,0], [colors['|0⟩'], colors['|+⟩'], colors['|e⟩']]

for j, xj in enumerate(x):
    for r, col, off in zip(rows, cols, offsets):
        ax.bar(xj+off, W[r,j], width=bar_w, color=col, edgecolor='k', linewidth=0.4, zorder=3)

ax.set_yscale('log'); ax.set_ylim(1e-4, 1.1)
ax.set_xlim(min(-10, x.min()-2), x.max()+2)
ax.set_xlabel("Energy (MHz)"); ax.set_ylabel("Composition")
ax.set_title(rf"Dressed-state Compositions")
ax.grid(True, axis='y', linestyle='--', alpha=0.35, zorder=0)
ax.legend(handles=[Patch(color=colors['|0⟩'], label='|0⟩'),
                   Patch(color=colors['|+⟩'], label='|+⟩'),
                   Patch(color=colors['|e⟩'], label='|e⟩')],
          loc='center', frameon=False, ncol=3)
plt.tight_layout(); plt.show()




