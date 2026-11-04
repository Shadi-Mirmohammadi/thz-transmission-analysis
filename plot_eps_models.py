"""Figure for eps_models.py: (a) eps' along c from the literature, the phonon model and the fit;
(b) measured bare transmission with the fitted linear and phonon-shape models.

usage: python plot_eps_models.py out_prefix_fits.npz tonouchi_points.csv figure_prefix
tonouchi_points.csv: f_THz,eps1_c,eps1_a read from Fig. 2b of Jiang et al., Photonics 9, 233 (2022).
"""
import sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from eps_models import lorentz_eps, slab_t, eps_model, NPAR

fits, tono, out = sys.argv[1:4]
z = np.load(fits)
f, T, keep, d_mm = z['f'], z['T'], z['keep'], float(z['d_mm'])
names = ['const', 'linear', 'quadratic', 'lorentz']
rms = dict(zip(names, z['rms']))
tp = np.genfromtxt(tono, delimiter=',', skip_header=1)

plt.rcParams.update({'font.family': 'Arial', 'font.size': 11, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': 0.8, 'legend.frameon': False})
VERM, BLUE, INK, GREY = '#D55E00', '#0072B2', '#222222', '#8C8C8C'
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.3), gridspec_kw={'width_ratios': [1.15, 1]})

# (a) permittivity along c
fg = np.linspace(0, 3, 400)
band = (f >= 0.3) & (f <= 0.65)
a1.axvspan(f[band][0], f[band][-1], color='#EEEEEE', lw=0)
a1.text(0.475, 15.15, 'your scan', ha='center', va='top', fontsize=9, color='#555555')
a1.plot(tp[:, 0], tp[:, 1], 'o', ms=4.5, mfc='white', mec=INK, mew=1, label='Tonouchi 2022, THz-TDS (dark)')
a1.plot(fg, lorentz_eps(fg, 'c').real, '-', color=BLUE, lw=2, label='phonon model, Schubert 2016')
p = z['lorentz']
a1.plot(fg, eps_model('lorentz', p, fg), '--', color=BLUE, lw=2,
        label=f'phonon shape fitted to your level ({p[0]:.2f} at 0 THz)')
fb = np.linspace(0.3, 0.65, 50)
a1.plot(fb, 11.86 + 0.367 * fb, '-', color=VERM, lw=3, label='your HFSS equation, 11.86 + 0.367 f')
a1.errorbar([0.03], [12.4], yerr=[0.4], fmt='s', ms=5, color=INK, capsize=3, lw=1,
            label='Fiedler 2019, static 12.4 ± 0.4')
a1.plot([0.2, 1.0], [12.4, 12.4], ':', color=INK, lw=1.5, label='Gopalan 2020, THz (no dispersion)')
a1.set_xlim(0, 3); a1.set_ylim(11, 15.2)
a1.set_xlabel('Frequency (THz)'); a1.set_ylabel("ε' along c")
a1.set_title('(a) Permittivity along c', loc='left', fontsize=11)
a1.legend(loc='upper left', fontsize=8.5, bbox_to_anchor=(0.0, 0.93))

# (b) transmission: data and two models
fm = f[(f >= 0.3) & (f <= 0.65)]
Tm = T[(f >= 0.3) & (f <= 0.65)]
a2.plot(fm * 1000, Tm, '-', color=GREY, lw=1.2, label='measured, bare 686 µm')
pl = z['linear']; k = NPAR['linear']
a2.plot(fm * 1000, slab_t(fm, eps_model('linear', pl, fm), np.clip(pl[k] + pl[k + 1] * fm, 0, None), d_mm),
        '-', color=VERM, lw=1.6, label=f'straight-line ε, rms {rms["linear"]:.3f}')
k = NPAR['lorentz']
a2.plot(fm * 1000, slab_t(fm, eps_model('lorentz', p, fm), np.clip(p[k] + p[k + 1] * fm, 0, None), d_mm),
        '--', color=BLUE, lw=1.6, label=f'phonon-shape ε, rms {rms["lorentz"]:.3f}')
for lo, hi in ((448, 495), (545, 572)):
    a2.axvspan(lo, hi, color='#F2F2F2', lw=0)
a2.text(471, 1.07, 'not fitted', ha='center', va='top', fontsize=8, color='#777777')
a2.text(558, 1.07, 'water', ha='center', va='top', fontsize=8, color='#777777')
a2.set_xlim(300, 650); a2.set_ylim(0.15, 1.08)
a2.set_xlabel('Frequency (GHz)'); a2.set_ylabel('Transmission (field)')
a2.set_title('(b) Both models fit the bare scan equally', loc='left', fontsize=11)
a2.legend(loc='lower left', fontsize=8.5, ncol=1)

fig.tight_layout()
fig.savefig(out + '.png', dpi=200); fig.savefig(out + '.pdf')
print('saved', out + '.png')
