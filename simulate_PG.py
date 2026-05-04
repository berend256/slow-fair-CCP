
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
simulate_PG.py
==============

Memory-safe Monte‑Carlo simulation of the coverage‑time process in PG(2,q),
with optional verbosity and return values for batch scripts.
"""

import numpy as np
import galois


# ==================================================
# Construct PG(2,q)
# ==================================================
def generate_pg2q(p, k):
    """
    Construct the incidence matrix A for PG(2,q).
    Returns:
        A : boolean numpy array of shape (n_lines, n_points)
            with A[i,j] = True iff point j lies on line i.
    """
    q = p**k
    GF = galois.GF(q)
    elems = GF.elements
    zero, one = GF(0), GF(1)

    # ---------- Points ----------
    pts = []
    # [1 : y : z]
    for y in elems:
        for z in elems:
            pts.append([one, y, z])
    # [0 : 1 : z]
    for z in elems:
        pts.append([zero, one, z])
    # [0 : 0 : 1]
    pts.append([zero, zero, one])

    P = GF(pts)
    n = P.shape[0]

    # ---------- Lines ----------
    lns = []
    # [1 : b : c]
    for b in elems:
        for c in elems:
            lns.append([one, b, c])
    # [0 : 1 : c]
    for c in elems:
        lns.append([zero, one, c])
    # [0 : 0 : 1]
    lns.append([zero, zero, one])

    L = GF(lns)
    assert L.shape[0] == n

    # ---------- Incidence ----------
    M = L @ P.T
    A = (M == 0)
    return A


# ==================================================
# Monte‑Carlo Simulation (memory‑minimal)
# ==================================================
def simulate_coupon(p, k, num_sims=10000, verbose=True):
    """
    Perform num_sims Monte‑Carlo runs of the coverage process on PG(2,q).

    Args:
        p, k     : integers defining q = p^k
        num_sims : number of independent trials
        verbose  : if True, print results; if False, suppress printing

    Returns:
        {
            "mean": float,
            "std": float,
            "trunc": int,
            "num_sims": num_sims
        }
    """

    q = p**k
    A = generate_pg2q(p, k)
    n = A.shape[0]

    max_steps = 20 * n  # safety cutoff

    draws = []
    rng = np.random.default_rng()  # fast modern RNG

    for sim in range(num_sims):
        covered = np.zeros(n, dtype=bool)

        for t in range(max_steps):
            line = rng.integers(0, n)
            covered |= A[line]

            if covered.all():
                draws.append(t + 1)
                break
        else:
            # truncated
            draws.append(max_steps)

    draws = np.array(draws, dtype=float)
    num_trunc = int(np.sum(draws == max_steps))

    mean_val = float(draws.mean())
    std_val  = float(draws.std())

    if verbose:
        print(f"\nPG(2,{q}) with n={n}")
        print(f"Average draws: {mean_val:.6f}")
        print(f"Std deviation: {std_val:.6f}")
        print(f"Truncated runs: {num_trunc}/{num_sims}")

    return {
        "mean": mean_val,
        "std": std_val,
        "trunc": num_trunc,
        "num_sims": num_sims
    }


# ==================================================
# Standalone mode
# ==================================================
if __name__ == "__main__":
    p = int(input("Enter prime p: "))
    k = int(input("Enter power k: "))
    simulate_coupon(p, k, num_sims=10000, verbose=True)

