# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 14:05:38 2025

@author: iontrap
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch  # <-- for custom legend

# Mean occupation number
nbar = 10.0

# Occupation number states
n = np.arange(0, 60)

# Choose a motional state n_target >> nbar
#n_target = 35  # Example: 3×n̄

# Thermal distribution
pth = (nbar**n) / ((nbar + 1)**(n + 1))
pth /= np.sum(pth)  # Normalize

# Color list: red if n == n_target, else steelblue
#colors = ['red' if ni == n_target else 'steelblue' for ni in n]

# Plot
plt.figure(figsize=(10, 4))
plt.bar(n, pth, width=0.8, edgecolor='black')

# Title and axes
plt.title(r"Probability of Occupation vs. Motional State", fontsize=14)
plt.xlabel("Motional State $n$", fontsize=12)
plt.ylabel("Probability $p_{\mathrm{th}}(n)$", fontsize=12)

# Custom legend handles
legend_handles = [
    Patch(color='steelblue', label=fr"$\bar{{n}}$ = {nbar}"),
    #Patch(color='red', label=f"Number of RSB Pulses = {n_target}")
]
plt.legend(handles=legend_handles, loc='upper right')

# Aesthetics
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()