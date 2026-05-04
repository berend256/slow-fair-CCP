# Slow Coupon Collecting under Fair Mechanisms

This repository contains all code used to generate the exact and numerical
results in the paper

D. Barak‑Pelleg and D. Berend,  
*Fano Geometry and Slow Coupon Collecting*.

The repository is intended as a reproducibility and data supplement to the
paper.

## Contents

- Exact recovery‑set enumeration for projective planes PG(2,q)
- Exact recovery‑set enumeration for affine planes AG(2,q)
- Exact expected coverage times via the Grunbaum–Yaakobi identity
- Monte Carlo simulations for larger parameters
- Exact and simulated full‑model baselines
- Exact computations for the star mechanism

## Generating exact projective‑plane results

Run:
```bash
python3 compute_recovery_sets_PG.py
```

and enter the following values when prompted:

- (p, k) = (2, 1) to generate results for q = 2
- (p, k) = (3, 1) to generate results for q = 3
- (p, k) = (2, 2) to generate results for q = 4

The script produces the files:

- PG_q2.csv
- PG_q3.csv
- PG_q4.csv

Each CSV file contains the recovery‑set counts by size together with the
exact expected coverage time (both as a rational number and as a decimal
approximation).

## Generating exact affine‑plane results

Analogously, run:
```bash
python3 compute_recovery_sets_AG.py
```

and enter the corresponding values of (p, k) to generate the files:

- AG_q2.csv
- AG_q3.csv
- AG_q4.csv

## Notes

Exact recovery‑set enumeration is feasible only for small values of q
(up to q = 4). Monte Carlo simulations for larger parameters may require
up to 10^8 iterations and can take several hours.
