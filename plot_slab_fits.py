"""Figure and table for slab_drude_fit.py results.

(a) measured |t| for every scan, with the fit for the fitted ones
(b) power absorption coefficient: point by point from the envelope (eps_L fixed) and from the fit
(c) extinction coefficient kappa, same way

usage: python plot_slab_fits.py out_prefix "label=fit.npz" ... [--limit "label=fit.npz=sigma"]
  --limit scans are drawn in (a) only; sigma (S/m) is their lower limit, shown in (b) and (c).
Writes out_prefix.png/.pdf and out_prefix_kappa_alpha.csv (model n, kappa, alpha per fitted scan).
"""
import argparse, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import brentq
from slab_drude_fit import t_rel, eps_complex, C

COLORS = ['#0072B2', '#56B4E9', '#D55E00', '#E69F00', '#009E73', '#CC79A7']

ap = argparse.ArgumentParser()
ap.add_argument('out'); ap.add_argument('fits', nargs='+')
ap.add_argument('--limit', action='append', default=[])
a = ap.parse_args()

plt.rcParams.update({'font.family': 'Arial', 'font.size': 11, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': 0.8, 'legend.frameon': False})
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.4))
rows, ci = [], 0


def nk(z, f, sigma):
    N = np.sqrt(eps_complex(f, z['p'][0], sigma, float(z['tand']), float(z['tau'])))
    return N.real, N.imag, 4 * np.pi * f * 1e12 * N.imag / C / 100


for item in a.fits:
    lab, path = item.split('=')
    z = np.load(path); col = COLORS[ci]; ci += 1
    f, T, keep, d = z['f'], z['T'], z['keep'], float(z['d'])
    a1.semilogy(f * 1000, T, '-', color=col, lw=0.6, alpha=0.8, label=f'{lab}, measured')
    a1.semilogy(z['fk'] * 1000, z['model_T'], '--', color='#222222', lw=1.2)
    # point by point: the sigma that reproduces each measured |t| with eps_L fixed
    idx = np.where(keep)[0][::50]
    sp = [brentq(lambda s: abs(t_rel(np.array([f[i]]), eps_complex(np.array([f[i]]), z['p'][0], s,
                 float(z['tand']), float(z['tau'])), d))[0] - T[i], 1e-3, 1e5) for i in idx]
    _, kp, ap_ = nk(z, f[idx], np.array(sp))
    fg = z['fg']; nr, kg, ag = nk(z, fg, z['p'][1])
    a2.plot(f[idx] * 1000, ap_, 'o', ms=3, mfc='none', mec=col)
    a2.plot(fg * 1000, ag, '-', color=col, lw=2, label=f'{lab}: σ = {z["p"][1]:.0f} S/m')
    a3.plot(f[idx] * 1000, kp, 'o', ms=3, mfc='none', mec=col)
    a3.plot(fg * 1000, kg, '-', color=col, lw=2, label=f'{lab}, ε = {z["p"][0]:.2f}')
    rows.append((lab, fg, nr, kg, ag))

for item in a.limit:
    lab, path, sig = item.split('=')
    z = np.load(path); col = COLORS[ci]; ci += 1
    a1.semilogy(z['f'] * 1000, z['T'], '-', color=col, lw=0.6, alpha=0.8, label=f'{lab}, measured')
    if not sig:                        # several scans can share one lower-limit line
        continue
    fg = z['fg']; _, kg, ag = nk(z, fg, float(sig))
    a2.plot(fg * 1000, ag, ':', color=col, lw=2, label=f'{lab.split(",")[0]}: σ ≥ {float(sig):.0f} S/m (lower limit)')
    a3.plot(fg * 1000, kg, ':', color=col, lw=2, label=f'{lab.split(",")[0]}, lower limit')

a1.plot([], [], '--', color='#222222', label='fit')
a1.set_xlim(200, 1000); a1.set_ylim(1e-3, 1)
a1.set_xlabel('Frequency (GHz)'); a1.set_ylabel('Transmission |t| (field)')
a1.set_title('(a) Measured transmission', loc='left', fontsize=11); a1.legend(fontsize=8, loc='upper right')
a2.set_xlim(200, 1000); a2.set_ylim(bottom=0)
a2.set_xlabel('Frequency (GHz)'); a2.set_ylabel('Absorption coefficient α (cm$^{-1}$, power)')
a2.set_title('(b) Absorption: points = each frequency, line = fit', loc='left', fontsize=11)
a2.legend(fontsize=8, loc='center right')
a3.set_xlim(200, 1000); a3.set_ylim(bottom=0)
a3.set_xlabel('Frequency (GHz)'); a3.set_ylabel('Extinction coefficient κ')
a3.set_title('(c) κ(ω)', loc='left', fontsize=11); a3.legend(fontsize=8, loc='upper right')
fig.tight_layout()
fig.savefig(a.out + '.png', dpi=200); fig.savefig(a.out + '.pdf')

hdr = 'f_GHz,' + ','.join(f'{r[0]} n,{r[0]} kappa,{r[0]} alpha_cm-1' for r in rows)
tab = np.column_stack([rows[0][1] * 1000] + [c for r in rows for c in r[2:]])
np.savetxt(a.out + '_kappa_alpha.csv', tab, delimiter=',', header=hdr, comments='', fmt='%.5g')
print('saved', a.out + '.png')
