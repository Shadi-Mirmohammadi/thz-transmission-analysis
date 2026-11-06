# CW-THz Transmission Analysis

MATLAB and Python scripts I use to analyze continuous-wave (photomixing) THz transmission of crystals and metasurfaces from 0.2 to 1 THz. I used them for Fe-doped β-Ga₂O₃ crystals and for the ITO resonator arrays I fabricated on them. Project page: [optically tunable THz metasurface](https://shadi-mirmohammadi.github.io/projects/thz-metasurface/).

## Files

| File | What it does |
|---|---|
| `Plot_THz_Transmission.m` | Transmission of a sample from its scan and an air reference, by dividing the two envelopes. |
| `slabfit.py` | Bare slab: the fringe spacing gives n·d, then a least-squares fit of n(f), κ(f) and the thickness to the coherent Fabry–Pérot slab formula. |
| `eps_models.py` | Bare β-Ga₂O₃ slab: fits four forms of ε′(f) (constant, straight line, ε₀ + K·f², and the phonon-model shape from Schubert et al., Phys. Rev. B 93, 125209 (2016)) with the thickness fixed, and prints the phonon-model ε along c, a and b. |
| `plot_eps_models.py` | Figure for `eps_models.py`: ε′ along c from the literature and the fits, and the measured transmission with two of the models. |
| `slab_drude_fit.py` | Conductive slab: fits the lattice ε and the DC conductivity to a sample/air scan pair using both the envelope and the fringe phase, then gives the carrier density for an assumed mobility, κ(f) and the absorption coefficient α(f). Options fix ε or fit the envelope only when the phase drifted. |
| `plot_slab_fits.py` | Figure and κ/α table for `slab_drude_fit.py`: measured transmission with fits, α and κ point by point and from the fit, and lower limits for samples that are too opaque to fit. |
| `compare_device_hfss.m` | Measured transmission of a patterned sample, as deposited and annealed, against an HFSS conductivity sweep. |
| `photocurrent_vs_power.m` | Photocurrent against pump power from a folder of CSV files: baseline from the 0 W files, steady-state statistics, log-log fits and the slope of each regime. |

## Input

THz scans are text files with the header line `Frequency (GHz)	THz Photocurrent (nA)` followed by two columns. HFSS results are CSV exports.

## Running

- MATLAB: put the scans in one folder, edit the file names at the top of the script, and run it.
- Python (numpy, scipy, matplotlib): `python slabfit.py air.txt sample.txt 705 fit.png` (reference scan, sample scan, thickness guess in µm, output figure).
- `python eps_models.py air.txt bare.txt 686 out` (thickness fixed, in µm), then `python plot_eps_models.py out_fits.npz points.csv figure`, where `points.csv` holds literature ε′ values read from a published figure (`f_THz,eps1_c,eps1_a`).
- `python slab_drude_fit.py air.txt sample.txt 650 out --band 0.35-0.95` (add `--eps 10.7 --amp-only` if the phase drifted), then `python plot_slab_fits.py figure "label=out_fit.npz" --limit "label=other_fit.npz=140"`.

The measurement data is not included.
