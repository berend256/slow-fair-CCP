
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
compute_recovery_sets_PG.py
===========================

Exact recovery-set enumeration for PG(2,q).

Given q = p^k, the program:
1. Generates PG(2,q) using the galois library.
2. Computes the incidence matrix A (lines x points).
3. Enumerates ALL subsets of lines.
4. Computes:
   - alpha[s] = number of recovery sets of size s
   - M        = largest size of a non-recovery set
5. Applies the Grunbaum–Yaakobi identity to compute E[T].

# q is specified at runtime via input p,k with q = p**k
Runs realistically for q = 2, 3, 4.
"""

import numpy as np
import galois
import itertools
import math
from math import comb

# -------------------------------------------------------------
# Generate PG(2,q): returns incidence matrix A (shape n_lines x n_points)
# -------------------------------------------------------------
def generate_PG2q(p, k):
    q = p**k
    GF = galois.GF(q)
    elems = GF.elements
    zero, one = GF(0), GF(1)

    # ---------- Points ----------
    pts = []
    # normalize representation: first non-zero coordinate is scaled to 1
    for y in elems:
        for z in elems:
            pts.append([one, y, z])       # (1, y, z)
    for z in elems:
        pts.append([zero, one, z])        # (0, 1, z)
    pts.append([zero, zero, one])         # (0, 0, 1)
    P = GF(pts)
    n = P.shape[0]

    # ---------- Lines ----------
    lns = []
    for b in elems:
        for c in elems:
            lns.append([one, b, c])       # (1, b, c)
    for c in elems:
        lns.append([zero, one, c])        # (0, 1, c)
    lns.append([zero, zero, one])         # (0, 0, 1)
    L = GF(lns)
    assert L.shape[0] == n  # n lines = n points

    # ---------- Incidence matrix A[i,j] = point j lies on line i ----------
    M = L @ P.T
    A = (M == 0)  # Boolean matrix
    return A


# -------------------------------------------------------------
# Enumerate recovery sets
# -------------------------------------------------------------
def compute_recovery_sets(A):
    """
    A: Boolean incidence matrix (N lines × n points)
    Returns:
        alpha: list where alpha[s] = number of recovery sets of size s
        M: largest size of non-recovery set
    """
    N = A.shape[0]
    n_points = A.shape[1]

    alpha = [0] * (N + 1)
    max_non_recovery = -1

    # Represent each line as a bitmask of covered points
    masks = []
    for i in range(N):
        mask = 0
        for j, val in enumerate(A[i]):
            if val:
                mask |= (1 << j)
        masks.append(mask)

    full_mask = (1 << n_points) - 1

    # Enumerate all subsets of the N lines
    for s in range(N + 1):
        for combo in itertools.combinations(range(N), s):
            union_mask = 0
            for idx in combo:
                union_mask |= masks[idx]

            if union_mask == full_mask:
                alpha[s] += 1
            else:
                if s > max_non_recovery:
                    max_non_recovery = s

    return alpha, max_non_recovery


# -------------------------------------------------------------
# GY identity
# -------------------------------------------------------------
def expected_coverage_time(alpha, M):
    """
    Apply the GY formula:
        E[T] = N (H_N - H_{N-M-1}) - sum_{s=0}^M [ alpha[s] / C(N-1, s) ].
    """
    N = len(alpha) - 1

    # Harmonic numbers
    def H(k):
        return sum(1.0 / j for j in range(1, k + 1))

    val = N * (H(N) - H(N - M - 1))
    for s in range(M + 1):
        if comb(N - 1, s) > 0:
            val -= alpha[s] / comb(N - 1, s)
    return val


# -------------------------------------------------------------
# Main driver
# -------------------------------------------------------------
if __name__ == "__main__":
    print("Exact recovery-set computation for PG(2,q)")
    p = int(input("Enter prime p: "))
    k = int(input("Enter power k: "))

    A = generate_PG2q(p, k)
    N = A.shape[0]
    print(f"\nGenerated PG(2,{p**k}) with N = {N} lines and points.")

    alpha, M = compute_recovery_sets(A)

    print("\nRecovery-set statistics:")
    print(f"M (largest non-recovery size) = {M}")
    print("alpha(s) values:")
    for s, val in enumerate(alpha):
        if val > 0:
            print(f"  s = {s}: {val}")


# -------------------------------------------------------------
# Print expected coverage time in rational form a/b ≈ decimal
# -------------------------------------------------------------
ET = expected_coverage_time(alpha, M)

# Convert to rational a/b using Python's Fraction
from fractions import Fraction
ET_frac = Fraction(ET).limit_denominator()

a, b = ET_frac.numerator, ET_frac.denominator
ET_float = float(ET)


print("\nExpected coverage time:")
print(f"  E[T] = {a}/{b} ≈ {float(ET):.6f}")

# -------------------------------------------------------------
# Optional CSV Saving
# -------------------------------------------------------------
save = input("\nSave results to CSV file? (y/n): ").strip().lower()

if save == "y":
    q = p**k
    filename = f"PG_q{q}.csv"
    with open(filename, "w") as f:
        # Header for the alpha table
        f.write("s,alpha_s\n")
        for s, val in enumerate(alpha):
            f.write(f"{s},{val}\n")
        # Expected value
        f.write(f"E_T_fraction,{a}/{b}\n")
        f.write(f"E_T_approx,{float(ET):.12f}\n")

    print(f"\nCSV results saved to {filename}\n")
else:
    print()
