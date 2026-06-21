"""
AnalyzeSweep.py

Performs a sweep of evolutions
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.colors import LogNorm
import random

basefilename = "data/SweepHiddenVector1"
with open(basefilename + ".dat") as f:
    data = f.readlines()
results = np.array([list(map(float, line.split(", "))) for line in data]).transpose()

mx    = results[0]
beta  = results[1]
r     = results[2]
gX    = results[3]
Ydm   = results[4]
Omegah2 = results[5]
Switch  = results[6]

data = np.array(list(zip(
    results[0], results[1], results[2],
    results[3], results[4], results[5], results[6]
)), dtype=[
    ('mx', 'f8'),
    ('beta', 'f8'),
    ('r', 'f8'),
    ('gX', 'f8'),
    ('Ydm', 'f8'),
    ('Omegah2', 'f8'),
    ('Switch', 'i4')
])

print("beta min:", beta.min())
print("beta max:", beta.max())
print("gX min:", gX.min())
print("gX max:", gX.max())

#Data without errors (Switch=0)

good_data = data[(data['Switch'] == 0) & (data['Omegah2'] > 0) & np.isfinite(data['Omegah2'])]

print("Total r in file:", len(np.unique(data['r'])))
print("Total r after filter:", len(np.unique(good_data['r'])))

unique_mx   = np.unique(good_data['mx'])
print("unique_mx:", unique_mx)
unique_beta = np.unique(good_data['beta'])
print("unique_beta:", unique_beta)
unique_r = np.unique(good_data['r'])
print("unique_r:", unique_r)

relic_rows = []
basefilename = "data/SweepHiddenVector_Relic"
separator = ', '
datafile = open(basefilename + ".csv", "w")
datafile.write("mx,beta,r,gX,Ydm,Omegah2,Switch\n")
target = 0.12


MAX_DEBUG_PLOTS = 10
debug_count = 0

valid_ternas = []

for r_u in unique_r:
    for mx_u in unique_mx:
        for beta_u in unique_beta:

            sub_data = good_data[
                (good_data['mx'] == mx_u) &
                (good_data['beta'] == beta_u) &
                (good_data['r'] == r_u)
            ]

            if len(sub_data) < 2:
                continue

            valid_ternas.append((r_u, mx_u, beta_u))

            # order gx before to look for points
            sub_data = np.sort(sub_data, order='gX')

            Omegah2_sub = sub_data['Omegah2']
            gX_sub = sub_data['gX']

            OmegaDiff = Omegah2_sub - target
            cross = np.where((OmegaDiff[:-1] * OmegaDiff[1:]) < 0)[0]

            if len(cross) > 1:
                print("multiple crossings for (r,mx,beta)=", r_u, mx_u, beta_u, "cross=", cross)

            for i in cross:
                # ilinear interpolation
                t = (target - Omegah2_sub[i])/(Omegah2_sub[i+1] - Omegah2_sub[i])
                gX_target = gX_sub[i] + t * (gX_sub[i+1] - gX_sub[i])

                #Ydm formulam
                Ydm_target = target/(2.8282*10**8 * mx_u)

                relic_rows.append((mx_u, beta_u, r_u, gX_target, Ydm_target, target, 0))

                text = (
                    str(mx_u) + separator
                    + str(beta_u) + separator
                    + str(r_u) + separator
                    + str(gX_target) + separator
                    + str(Ydm_target) + separator
                    + str(target) + separator
                    + str(0) + "\n"
                )
                datafile.write(text)

datafile.close()

relic_data = np.array(relic_rows, dtype=good_data.dtype)
print("relic_data:", relic_data)

# --- Debug plots: elegir 10 ternas al azar y plotear Omega vs gX ---
if len(valid_ternas) == 0:
    print("No valid ternas to debug-plot.")
else:
    N_DEBUG = min(MAX_DEBUG_PLOTS, len(valid_ternas))
    debug_selection = random.sample(valid_ternas, N_DEBUG)

    for j, (r_u, mx_u, beta_u) in enumerate(debug_selection, start=1):
        sub_data = good_data[
            (good_data['mx'] == mx_u) &
            (good_data['beta'] == beta_u) &
            (good_data['r'] == r_u)
        ]
        if len(sub_data) < 2:
            continue

        sub_data = np.sort(sub_data, order='gX')
        Omegah2_sub = sub_data['Omegah2']
        gX_sub = sub_data['gX']

        gx_fine = np.linspace(gX_sub.min(), gX_sub.max(), 200)
        om_fine = np.interp(gx_fine, gX_sub, Omegah2_sub)

        plt.figure()
        plt.plot(gX_sub, Omegah2_sub, "o", label="data")
        plt.plot(gx_fine, om_fine, "-", label="interp (linear)")
        plt.axhline(target, ls="--", lw=1.0, label=f"target={target}")
        plt.xlabel("gX")
        plt.ylabel("Omega h^2")
        plt.title(f"Random Debug {j}/{N_DEBUG}: r={r_u:.3g}, mx={mx_u:.3g}, sin(beta)={np.sin(beta_u):.3e}")
        plt.tight_layout()

        fname = f"data/debug_Omega_vs_gX_{j:02d}_random.png"
        plt.savefig(fname, dpi=200)
        plt.close()
        print("Saved:", fname)

# --- PLOTS: para cada r fijo -> puntos (mx, sin(beta)) coloreados por gX ---

if len(relic_data) == 0:
    print("No relic points found.")
    exit()



for r_u in np.unique(relic_data['r']):
    sel = (relic_data['r'] == r_u)
    if np.sum(sel) == 0:
        continue

    mx_vals_all = relic_data['mx'][sel]
    sinb_vals_all = np.sin(relic_data['beta'][sel])
    gx_vals_all = relic_data['gX'][sel]

    # filtro para poder usar escala log en y
    mask = (sinb_vals_all > 0.0) & (sinb_vals_all < 0.9)

    mx_vals = mx_vals_all[mask]
    sinb_vals = sinb_vals_all[mask]
    gx_vals = gx_vals_all[mask]

    if len(mx_vals) == 0:
        print(f"No valid points with 0 < sin(beta) < 1 for r = {r_u:.3g}")
        continue

    plt.figure()

    sc = plt.scatter(
        mx_vals,
        sinb_vals,
        c=gx_vals,
        norm=LogNorm(),
        cmap="viridis",
        s=60
    )


    plt.colorbar(sc, label="gX (Omega h^2 = 0.12)")

    plt.xlabel("mx [GeV]")
    plt.ylabel("sin(beta)")
    plt.title(f"r = {r_u:.3g}")

    plt.tight_layout()
    filename = f"data/gX_points_r{r_u:.3g}.pdf"
    plt.savefig(filename, dpi=200)
    plt.close()
    print("Saved:", filename)

print("Done.")