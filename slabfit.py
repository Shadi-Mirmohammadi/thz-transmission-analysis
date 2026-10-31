"""CW-THz (photomixing) transmission of a bare slab: envelope ratio, fringe analysis, and
least-squares fit of n(f), kappa(f), thickness to the coherent Fabry-Perot slab formula."""
import sys, numpy as np
from scipy.signal import hilbert
from scipy.optimize import least_squares
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ref_f, samp_f = sys.argv[1], sys.argv[2]; d0_um = float(sys.argv[3]); out = sys.argv[4]
def load(p):
    a = np.loadtxt(p, skiprows=1); a = a[np.argsort(a[:, 0])]; return a[:, 0] / 1000.0, a[:, 1]   # THz, nA
fr, yr = load(ref_f); fs, ys = load(samp_f)
f = np.linspace(max(fr[0], fs[0]), min(fr[-1], fs[-1]), 8000)
yr = np.interp(f, fr, yr); ys = np.interp(f, fs, ys)
def env(y, win=41):
    a = np.abs(hilbert(y - np.convolve(y, np.ones(win) / win, 'same')))
    return np.convolve(a, np.ones(win) / win, 'same')
T = env(ys) / env(yr)
m = (f > f[0] + 0.01) & (f < f[-1] - 0.01)
f, T = f[m], T[m]
# fringe period -> n*d
from scipy.signal import find_peaks
pk, _ = find_peaks(T, distance=int(0.02 / (f[1] - f[0])), prominence=0.03)
fsr = np.median(np.diff(f[pk])); nd_mm = 299.792458 / (2 * fsr * 1e3)   # c/(2*FSR) in mm
print(f'fringe spacing {fsr*1e3:.1f} GHz -> n*d = {nd_mm:.3f} mm; with d = {d0_um} um -> n = {nd_mm/(d0_um/1000):.3f}')
def slab_T(f, n, k, d_mm):
    N = n + 1j * k; k0 = 2 * np.pi * f * 1e12 / 299792458.0 * 1e-3 * 1e0   # rad/mm: 2*pi*f/c with c in mm/s
    k0 = 2 * np.pi * f * 1e12 / 2.99792458e11
    t12 = 2 / (1 + N); t21 = 2 * N / (1 + N); r = (N - 1) / (N + 1)
    ph = np.exp(1j * N * k0 * d_mm)
    return np.abs(t12 * t21 * ph / (1 - r**2 * ph**2) * np.exp(-1j * k0 * d_mm))
def model(p, f):
    n = p[0] + p[1] * (f - 0.6) + p[2] * (f - 0.6)**2
    k = np.abs(p[3] + p[4] * (f - 0.6))
    return slab_T(f, n, k, p[5])
p0 = [nd_mm / (d0_um / 1000), 0.0, 0.0, 0.01, 0.0, d0_um / 1000]
res = least_squares(lambda p: model(p, f) - T, p0, bounds=([2, -2, -5, 0, -0.5, 0.3], [5, 2, 5, 0.2, 0.5, 1.5]))
p = res.x; n = p[0] + p[1] * (f - 0.6) + p[2] * (f - 0.6)**2; k = np.abs(p[3] + p[4] * (f - 0.6))
eps1 = n**2 - k**2; tand = 2 * n * k / eps1
print('fit: n(f) = %.4f %+.4f*(f-0.6) %+.4f*(f-0.6)^2 ; kappa = %.4f %+.4f*(f-0.6) ; d = %.1f um ; rms resid %.4f' % (p[0], p[1], p[2], p[3], p[4], p[5] * 1000, np.sqrt(np.mean(res.fun**2))))
for ff in (0.3, 0.4, 0.5, 0.6, 0.8, 1.0):
    if f[0] <= ff <= f[-1]:
        i = np.argmin(abs(f - ff)); print(f'  {ff:.1f} THz: n={n[i]:.3f} eps={eps1[i]:.2f} tand={tand[i]:.4f} alpha={4*np.pi*k[i]*ff*1e12/2.99792458e10:.2f} /cm')
np.savetxt(out + '_dataset.csv', np.c_[f, eps1, tand], delimiter=',', header='f_THz,eps_r,tan_delta', comments='')
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
ax[0].plot(f, T, 'b', lw=1, label='experiment (envelope ratio)'); ax[0].plot(f, model(p, f), 'r', lw=1, label='slab fit'); ax[0].set_ylabel('amplitude transmission'); ax[0].legend()
ax[1].plot(f, eps1, 'k'); ax[1].set_ylabel("eps'"); ax2 = ax[1].twinx(); ax2.plot(f, tand, 'g'); ax2.set_ylabel('tan delta', color='g'); ax[1].set_xlabel('Frequency (THz)')
fig.tight_layout(); fig.savefig(out + '.png', dpi=110)
