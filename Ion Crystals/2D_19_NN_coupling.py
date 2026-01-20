import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from math import sqrt, exp

# ---------- Hex (triangular) lattice: axial coords ----------
R = 2                    # hex "radius": 1 + 6 + 12 = 19 sites
a = 1.0                  # lattice spacing

# axial (q,r) -> 2D coordinates (pointy-top hex basis)
def axial_to_xy(q, r, a=1.0):
    x = a * (sqrt(3) * (q + r/2))
    y = a * (3/2) * r
    return np.array([x, y])

# cube distance for hex grid
def hex_distance(q1, r1, q2, r2):
    s1 = -q1 - r1
    s2 = -q2 - r2
    return (abs(q1-q2) + abs(r1-r2) + abs(s1-s2)) // 2

# generate axial coordinates in a hex of radius R
axial = []
for q in range(-R, R+1):
    rmin = max(-R, -q-R)
    rmax = min(R, -q+R)
    for r in range(rmin, rmax+1):
        axial.append((q, r))

# map to xy
pts = np.array([axial_to_xy(q, r, a) for q, r in axial])

# ---------- Build edge lists ----------
# nearest neighbors (hex distance = 1)
nn_edges = []
all_edges = []
for (i, (q1, r1)), (j, (q2, r2)) in combinations(list(enumerate(axial)), 2):
    d_hex = hex_distance(q1, r1, q2, r2)
    all_edges.append((i, j))
    if d_hex == 1:
        nn_edges.append((i, j))

# Euclidean distances (for fading all-to-all lines)
dists = {(i, j): np.linalg.norm(pts[i]-pts[j]) for (i, j) in all_edges}
maxd = max(dists.values())

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(5, 5))

# brighter long-range edges (alpha decays with distance, but higher baseline)
for (i, j) in all_edges:
    x = [pts[i,0], pts[j,0]]
    y = [pts[i,1], pts[j,1]]
    # fade approximately exponentially with distance
    alpha = 0.12 + 0.35 * exp(-dists[(i, j)] / (a*1.3))
    ax.plot(x, y, lw=0.9, color=(0.4, 0.4, 0.4, alpha))  # light grey

# dark nearest-neighbor edges
for (i, j) in nn_edges:
    x = [pts[i,0], pts[j,0]]
    y = [pts[i,1], pts[j,1]]
    ax.plot(x, y, lw=2.3, color="#7b1c1c")  # dark red

# nodes
ax.scatter(pts[:,0], pts[:,1], s=120, color="black", zorder=3)

ax.set_aspect("equal")
ax.axis("off")
plt.tight_layout()
plt.show()


