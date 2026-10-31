# CW-THz Transmission Analysis

MATLAB and Python scripts I use to analyze continuous-wave (photomixing) THz transmission of crystals and metasurfaces from 0.2 to 1 THz. I used them for Fe-doped β-Ga₂O₃ crystals and for the ITO resonator arrays I fabricated on them. Project page: [optically tunable THz metasurface](https://shadi-mirmohammadi.github.io/projects/thz-metasurface/).

## Files

| File | What it does |
|---|---|
| `Plot_THz_Transmission.m` | Transmission of a sample from its scan and an air reference, by dividing the two envelopes. |
| `slabfit.py` | Bare slab: the fringe spacing gives n·d, then a least-squares fit of n(f), κ(f) and the thickness to the coherent Fabry–Pérot slab formula. |
| `compare_device_hfss.m` | Measured transmission of a patterned sample, as deposited and annealed, against an HFSS conductivity sweep. |
| `photocurrent_vs_power.m` | Photocurrent against pump power from a folder of CSV files: baseline from the 0 W files, steady-state statistics, log-log fits and the slope of each regime. |

## Input

THz scans are text files with the header line `Frequency (GHz)	THz Photocurrent (nA)` followed by two columns. HFSS results are CSV exports.

## Running

- MATLAB: put the scans in one folder, edit the file names at the top of the script, and run it.
- Python (numpy, scipy, matplotlib): `python slabfit.py air.txt sample.txt 705 fit.png` (reference scan, sample scan, thickness guess in µm, output figure).

The measurement data is not included.
