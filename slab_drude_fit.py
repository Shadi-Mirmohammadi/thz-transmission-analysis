"""Conductive slab: fit lattice eps' and DC conductivity to a CW-THz scan pair (sample, air),
using both the envelope (amplitude) and the fringe phase.

Model: eps(f) = eps_L*(1 + i*tand_L) + i*sigma(f)/(w*eps0), sigma(f) = sigma0/(1 - i*w*tau)
  tau follows from the assumed mobility (tau = mu*m*/q); below 1 THz it hardly matters.
Outputs sigma0, eps_L, carrier density n = sigma0/(q*mu), kappa(f) and the power absorption
coefficient alpha(f) = 4*pi*f*kappa/c.

usage: python slab_drude_fit.py air.txt sample.txt thickness_um out_prefix
       [--band 0.22-0.95] [--exclude 440-456,545-570,742-762] [--mu 100] [--tand 0.009] [--snr 5]
"""
import argparse, numpy as np
from scipy.signal import hilbert
from scipy.optimize import least_squares
from eps_models import load

C = 2.99792458e8; EPS0 = 8.8541878e-12; Q = 1.602177e-19; ME = 9.10938e-31
MSTAR = 0.28                   # electron effective mass of beta-Ga2O3, units of m0


def analytic(path, f, win=41):
    fr, y = load(path)
    y = np.interp(f, fr, y)
    y = y - np.convolve(y, np.ones(win) / win, 'same')
    z = hilbert(y)
    env = np.convolve(np.abs(z), np.ones(win) / win, 'same')
    return z, env


def eps_complex(f_thz, eps_l, sigma0, tand, tau):
    w = 2 * np.pi * f_thz * 1e12
    sig = sigma0 / (1 - 1j * w * tau)
    return eps_l * (1 + 1j * tand) + 1j * sig / (w * EPS0)


def t_rel(f_thz, eps, d_m):
    """Slab transmission relative to the same length of air, all echoes included."""
    N = np.sqrt(eps)
    k0 = 2 * np.pi * f_thz * 1e12 / C
    r = (N - 1) / (N + 1)
    ph = np.exp(1j * N * k0 * d_m)
    return (1 - r**2) * ph / (1 - r**2 * ph**2) * np.exp(-1j * k0 * d_m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('air'); ap.add_argument('sample'); ap.add_argument('thickness_um', type=float)
    ap.add_argument('out')
    ap.add_argument('--band', default='0.22-0.95')
    ap.add_argument('--exclude', default='440-456,545-570,742-762')   # water lines, GHz
    ap.add_argument('--mu', type=float, default=100.0)                # cm^2/Vs, only for n and tau
    ap.add_argument('--tand', type=float, default=0.009)              # lattice loss, fixed
    ap.add_argument('--snr', type=float, default=5.0)
    ap.add_argument('--eps', type=float, help='fix eps_L (e.g. when the phase drifted)')
    ap.add_argument('--amp-only', action='store_true', help='fit the envelope only')
    a = ap.parse_args()
    d = a.thickness_um * 1e-6
    tau = a.mu * 1e-4 * MSTAR * ME / Q

    f = np.arange(0.2, 0.99996, 0.00004)
    za, ea = analytic(a.air, f)
    zs, es = analytic(a.sample, f)
    floor = np.median(es[f > 0.96]) if np.median(es[f > 0.96]) < 0.02 * es.max() else 0.0
    lo, hi = (float(v) for v in a.band.split('-'))
    keep = (f > lo) & (f < hi) & (es > a.snr * max(floor, 1e-12))
    for rng in a.exclude.split(','):
        if rng:
            x0, x1 = (float(v) / 1000 for v in rng.split('-'))
            keep &= ~((f >= x0) & (f <= x1))

    T = es / ea
    ratio = zs * np.conj(za)
    ratio = np.convolve(ratio / np.abs(ratio), np.ones(501) / 501, 'same')
    psi = np.unwrap(np.angle(ratio[keep]))
    fk, Tk = f[keep], T[keep]

    def full(p):                       # p = [eps_L, sigma0, phase offset]
        return [a.eps, p[0], p[1]] if a.eps else p

    def model(p):
        p = full(p)
        tr = t_rel(fk, eps_complex(fk, p[0], p[1], a.tand, tau), d)
        return np.log(np.abs(tr)), np.unwrap(np.angle(tr))

    def resid(p):
        lt, ph = model(p)
        r = (np.log(Tk) - lt) / 0.05
        return r if a.amp_only else np.concatenate([r, (psi - ph - full(p)[2]) / 0.1])

    best = None
    for e0 in ([None] if a.eps else (9.5, 10.5, 11.5, 12.5)):
        for s0 in (1, 10, 100, 1000):
            for c0 in np.arange(-np.pi, np.pi, np.pi / 2) + psi[0]:
                x0 = ([] if a.eps else [e0]) + [s0, c0]
                lo_b = ([] if a.eps else [6]) + [0, -1e3]
                hi_b = ([] if a.eps else [16]) + [1e5, 1e3]
                r = least_squares(resid, x0, bounds=(lo_b, hi_b))
                if best is None or r.cost < best.cost:
                    best = r
    J = best.jac; cov = np.linalg.pinv(J.T @ J) * (2 * best.cost / max(1, len(best.fun) - len(best.x)))
    err = np.sqrt(np.diag(cov))
    p = np.array(full(best.x), float)
    if a.eps:
        err = np.concatenate([[0.0], err])
    n_cm3 = p[1] / (Q * a.mu * 1e-4) / 1e6
    lt, ph = model(best.x)
    rms_T = np.sqrt(np.mean((Tk - np.exp(lt))**2)); rms_ph = np.sqrt(np.mean((psi - ph - p[2])**2))

    fg = np.linspace(0.2, 1.0, 401)
    N = np.sqrt(eps_complex(fg, p[0], p[1], a.tand, tau))
    alpha = 4 * np.pi * fg * 1e12 * N.imag / C / 100            # 1/cm, power
    print(f"eps_L {p[0]:.2f} +- {err[0]:.2f}   sigma0 {p[1]:.3g} +- {err[1]:.2g} S/m   "
          f"n {n_cm3:.2e} cm^-3 (mu {a.mu:g})   band {fk[0]:.3f}-{fk[-1]:.3f} THz   "
          f"rms T {rms_T:.4f}  rms phase {rms_ph:.2f} rad   floor {floor:.3f} nA")
    for x in (0.3, 0.5, 0.7, 0.9):
        i = np.argmin(abs(fg - x))
        print(f"   {x:.1f} THz: n {N.real[i]:.3f}  kappa {N.imag[i]:.4f}  alpha {alpha[i]:.1f} 1/cm")
    np.savez(a.out + '_fit.npz', f=f, T=T, keep=keep, psi=psi, fk=fk, p=p, err=err, n_cm3=n_cm3,
             mu=a.mu, tand=a.tand, tau=tau, d=d, fg=fg, nr=N.real, kappa=N.imag, alpha=alpha,
             model_T=np.exp(lt), model_ph=ph + p[2], floor=floor)


if __name__ == '__main__':
    main()
