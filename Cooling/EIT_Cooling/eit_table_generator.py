# laser_params_plot_table_simple.py
import importlib.util
import math
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# 1) Full path to your parameter file (.py)
# ------------------------------------------------------------
file_path = r"Z:\Users\Wes\Code\Cooling\EIT_Cooling\final_eit_master_theory.py"

# ------------------------------------------------------------
# 2) Load the module from file
# ------------------------------------------------------------
spec = importlib.util.spec_from_file_location("eit_params", file_path)
src = importlib.util.module_from_spec(spec)
spec.loader.exec_module(src)

# ------------------------------------------------------------
# 3) Define the parameters to display
# ------------------------------------------------------------
TWOPI = 2 * math.pi
def to_MHz(omega_rad_s): return omega_rad_s / (TWOPI * 1e6)

params = [
    ("Omega_s_minus"),
    ("Omega_s_plus"),
    ("Omega_pi"),
    ("Delta_d"),
    ("delta_p"),       
]

rows = []
for key, sym in params:
    if hasattr(src, key):
        val = getattr(src, key)
        if np.isscalar(val):
            val_MHz = to_MHz(float(val))
            # Format adaptively: 1 decimal place if between 0.1–10, else 3
            if abs(val_MHz) < 10:
                val_str = f"{val_MHz:.1f} MHz"
            else:
                val_str = f"{val_MHz:.1f} MHz"
            rows.append([sym, val_str])

# ------------------------------------------------------------
# 4) Display table in Matplotlib
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5, len(rows)*0.6 + 1))
ax.axis("off")

table = ax.table(
    cellText=rows,
    cellLoc="center",
    loc="center"
)

table.auto_set_font_size(False)
table.set_fontsize(13)
table.scale(0.8, 1.4)

ax.set_title("Laser Parameters", fontsize=14, weight="bold", pad=20)
plt.tight_layout()
plt.show()

