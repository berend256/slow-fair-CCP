
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
simulate_AG.py
==============

Monte‑Carlo simulation of the coverage process on AG(2,q),
with optional verbosity and return values for batch scripts.
"""

import numpy as np
import galois


# ======================================================
# Construct AG(2,q)
# ======================================================

def generate_AG2q(p, k):
    q = p**k
    GF = galois.GF(q)
    elems = GF.elements

    # Points: all (x, y)
    points = [(x, y) for x in elems for y in elems]
    n_points = len(points)  # q^2

    # Lines: non‑vertical + vertical
    lines = []

    # non‑vertical: y = a x + b
    for a in elems:
        for b in elems:
            lines.append(("nv", a, b))

    # vertical: x = c
    for c in elems:
        lines.append(("v", c))

    n_lines = len(lines)  # q^2 + q

    A = np.zeros((n_lines, n_points), dtype=bool)

    for i, line in enumerate(lines):
        if line[0] == "nv":
            _, a, b = line
            for j, (x, y) in enumerate(points):
                if y == a*x + b:
                    A[i, j] = True
        else:
            _, c = line
            for j, (x, y) in enumerate(points):
                if x == c:
                    A[i, j] = True

    return A


# ======================================================
# Monte‑Carlo simulation (memory‑minimal)
# ======================================================

def simulate_coupon(p, k, num_sims=10000, verbose=True):
    """
    Run num_sims Monte‑Carlo trials for AG(2,q).

    Returns a dictionary:
        {"mean": float, "std": float, "trunc": int, "num_sims": num_sims}
    """

    q = p**k
    A = generate_AG2q(p, k)
    n_points = A.shape[1]
    n_lines = A.shape[0]

    max_steps = 20 * n_points  # safety cutoff

    draws = []
    rng = np.random.default_rng()

    for _ in range(num_sims):
        covered = np.zeros(n_points, dtype=bool)

        for t in range(max_steps):
            line = rng.integers(0, n_lines)
            covered |= A[line]

            if covered.all():
                draws.append(t + 1)
                break
        else:
            draws.append(max_steps)

    draws = np.array(draws, dtype=float)
    num_trunc = int(np.sum(draws == max_steps))

    mean_val = float(draws.mean())
    std_val  = float(draws.std())

    if verbose:
        print(f"\nAG(2,{q}) with n={n_points}")
        print(f"Average draws: {mean_val:.6f}")
        print(f"Std deviation: {std_val:.6f}")
        print(f"Truncated: {num_trunc}/{num_sims}")

    return {
        "mean": mean_val,
        "std": std_val,
        "trunc": num_trunc,
        "num_sims": num_sims
    }


# ======================================================
# Standalone interface
# ======================================================

if __name__ == "__main__":
    p = int(input("Enter prime p: "))
    k = int(input("Enter power k: "))
    simulate_coupon(p, k, num_sims=10000, verbose=True)

