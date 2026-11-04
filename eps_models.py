"""Bare beta-Ga2O3 slab: fit eps'(f) models to a CW-THz transmission scan, and compare with the
phonon (Lorentz) model of Schubert et al., Phys. Rev. B 93, 125209 (2016).

Models for eps'(f), each with tan_delta = t0 + t1*f (f in THz), thickness fixed:
  const      eps0
  linear     eps0 + s*f
  quadratic  eps0 + K*f^2           (low-frequency limit of the Lorentz model)
  lorentz    eps0 + [eps_c(f) - eps_c(0)] from Schubert's phonon parameters (only eps0 is fitted)

usage: python eps_models.py air.txt bare.txt thickness_um out_prefix [--tonouchi points.csv]
       [--exclude 448-495,545-572]   (GHz ranges left out of the fit)
"""
import sys, argparse, numpy as np
from scipy.signal import hilbert
from scipy.optimize import least_squares

C_MM = 2.99792458e11          # speed of light, mm/s
THZ_TO_CM = 33.35641          # 1 THz in cm^-1

# Schubert et al. 2016, Table II (TO modes) and Table IV (eps_inf).
# Axes: x || c, y || a-hat (perpendicular to c in the a-c plane), z || b. alpha = angle from c.
BU = dict(A=[256.4, 426.8, 820.3, 792.8, 358.3, 161.7, 485.6, 520.7],
          w=[743.5, 692.4, 572.5, 432.5, 356.8, 279.1, 262.3, 213.7],
          g=[10.4, 6.4, 12.3, 10.0, 3.79, 1.85, 1.75, 1.98],
          alpha=[48.7, 5.4, 106.0, 21.0, 144.0, 0.0, 158.5, 80.9])
AU = dict(A=[542.5, 718.8, 579.8, 72.0], w=[663.2, 448.6, 296.6, 154.8], g=[3.23, 10.28, 14.31, 2.1])
EPS_INF = dict(c=3.89, a=2.90, b=3.87)


def lorentz_eps(f_thz, axis):
    """Complex permittivity along c, a-hat or b from the phonon model (no free carriers)."""
    w = np.asarray(f_thz) * THZ_TO_CM
    if axis == 'b':
        modes, wts = AU, np.ones(4)
    else:
        al = np.radians(BU['alpha'])
        modes, wts = BU, (np.cos(al)**2 if axis == 'c' else np.sin(al)**2)
    eps = np.full(w.shape, EPS_INF[axis], dtype=complex)
    for A, w0, g, wt in zip(modes['A'], modes['w'], modes['g'], wts):
        eps += wt * A**2 / (w0**2 - w**2 - 1j * w * g)
    return eps


def load(path):
    a = np.loadtxt(path, skiprows=1)
    a = a[np.argsort(a[:, 0])]
    return a[:, 0] / 1000.0, a[:, 1]      # THz, nA


def envelope(y, win=41):
    a = np.abs(hilbert(y - np.convolve(y, np.ones(win) / win, 'same')))
    return np.convolve(a, np.ones(win) / win, 'same')


def slab_t(f, eps1, tand, d_mm):
    """|t| of a slab in air at normal incidence, all echoes included."""
    N = np.sqrt(eps1 * (1 + 1j * tand))
    k0 = 2 * np.pi * f * 1e12 / C_MM
    r = (N - 1) / (N + 1)
    ph = np.exp(1j * N * k0 * d_mm)
    return np.abs((1 - r**2) * ph / (1 - r**2 * ph**2))


def eps_model(name, p, f):
    if name == 'const':
        return p[0] + 0 * f
    if name == 'linear':
        return p[0] + p[1] * f
    if name == 'quadratic':
        return p[0] + p[1] * f**2
    if name == 'lorentz':
        shape = lorentz_eps(f, 'c').real - lorentz_eps(0.0, 'c').real
        return p[0] + shape
    raise ValueError(name)


NPAR = dict(const=1, linear=2, quadratic=2, lorentz=1)


def fit(name, f, T, d_mm):
    k = NPAR[name]

    def resid(p):
        tand = np.clip(p[k] + p[k + 1] * f, 0, None)
        return slab_t(f, eps_model(name, p, f), tand, d_mm) - T

    p0 = [12.0] + [0.1] * (k - 1) + [0.002, 0.0]
    lo = [8.0] + [-5.0] * (k - 1) + [0.0, -0.05]
    hi = [16.0] + [5.0] * (k - 1) + [0.05, 0.05]
    best = None
    for e0 in np.arange(10.5, 13.6, 0.1):      # the fringe order makes the problem multi-modal
        p0[0] = e0
        r = least_squares(resid, p0, bounds=(lo, hi))
        if best is None or r.cost < best.cost:
            best = r
    return best.x, np.sqrt(np.mean(best.fun**2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('air'); ap.add_argument('bare'); ap.add_argument('thickness_um', type=float)
    ap.add_argument('out'); ap.add_argument('--tonouchi'); ap.add_argument('--exclude', default='448-495,545-572')
    a = ap.parse_args()
    d_mm = a.thickness_um / 1000

    fr, yr = load(a.air); fs, ys = load(a.bare)
    f = np.linspace(max(fr[0], fs[0]), min(fr[-1], fs[-1]), 8000)
    T = envelope(np.interp(f, fs, ys)) / envelope(np.interp(f, fr, yr))
    keep = (f > f[0] + 0.01) & (f < f[-1] - 0.01)
    for rng in a.exclude.split(','):
        if rng:
            lo, hi = (float(v) / 1000 for v in rng.split('-'))
            keep &= ~((f >= lo) & (f <= hi))
    fit_f, fit_T = f[keep], T[keep]

    rows = []
    for name in ('const', 'linear', 'quadratic', 'lorentz'):
        p, rms = fit(name, fit_f, fit_T, d_mm)
        k = NPAR[name]
        e = [eps_model(name, p, np.array(x)) for x in (0.3, 0.5, 0.65)]
        rows.append((name, p, rms))
        print(f'{name:10s} params {np.round(p[:k], 4)}  tand = {p[k]:.5f} + {p[k+1]:.5f}*f  rms {rms:.4f}  '
              f'eps(0.3/0.5/0.65 THz) = {e[0]:.2f}/{e[1]:.2f}/{e[2]:.2f}')

    print('\nSchubert phonon model, real eps and tan_delta:')
    for ax in ('c', 'a', 'b'):
        vals = [lorentz_eps(x, ax) for x in (0.0, 0.3, 0.65, 1.0, 3.0)]
        print(f'  {ax}: ' + '  '.join(f'{v.real:.3f}' for v in vals) + '   (0/0.3/0.65/1/3 THz)'
              f'   tand at 0.5 THz {lorentz_eps(0.5, ax).imag / lorentz_eps(0.5, ax).real:.1e}')

    np.savez(a.out + '_fits.npz', f=f, T=T, keep=keep, d_mm=d_mm,
             **{n: p for n, p, _ in rows}, rms=np.array([r for _, _, r in rows]))


if __name__ == '__main__':
    main()
