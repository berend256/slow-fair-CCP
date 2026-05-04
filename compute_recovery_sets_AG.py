
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
compute_recovery_sets_AG.py
===========================

Exact recovery-set enumeration for AG(2,q).

Given q = p^k, this program:
1. Generates AG(2,q) using the galois library.
2. Builds the incidence matrix A (lines x points).
3. Enumerates all subsets of lines.
4. Computes:
     - alpha[s] = number of recovery sets of size s
     - M        = largest size of a non-recovery set
5. Applies the Grunbaum–Yaakobi identity to compute E[T].

# q is specified at runtime via input p,k with q = p**k
Runs realistically for q = 2 and q = 3.
q = 4 gives 2^(16+4)=1,048,576 subsets: borderline but possible.
"""

import numpy as np
import galois
import itertools
from math import comb
from fractions import Fraction


# -------------------------------------------------------------
# Generate AG(2,q): return incidence matrix A
# -------------------------------------------------------------
def generate_AG2q(p, k):
    q = p**k
    GF = galois.GF(q)
    elems = GF.elements
    zero = GF(0)

    # ---------------- Points (x,y) ----------------
    points = [(x, y) for x in elems for y in elems]
    n_points = len(points)

    # ---------------- Lines ----------------
    lines = []

    # Non-vertical: y = a x + b
    for a in elems:
        for b in elems:
            lines.append(("non-vertical", a, b))

    # Vertical: x = c
    for c in elems:
        lines.append(("vertical", c))

    n_lines = len(lines)

    # ---------------- Incidence matrix ----------------
    A = np.zeros((n_lines, n_points), dtype=bool)

    for i, line in enumerate(lines):
        if line[0] == "non-vertical":
            _, a, b = line
            for j, (x, y) in enumerate(points):
                if y == a*x + b:
                    A[i, j] = True
        else:  # vertical
            _, c = line
            for j, (x, y) in enumerate(points):
                if x == c:
                    A[i, j] = True

    return A


# -------------------------------------------------------------
# Compute recovery-set statistics
# -------------------------------------------------------------
def compute_recovery_sets(A):
    N = A.shape[0]       # number of lines
    n_points = A.shape[1]
    full_mask = (1 << n_points) - 1

    # Convert each row of A to a bitmask
    masks = []
    for i in range(N):
        mask = 0
        for j, val in enumerate(A[i]):
            if val:
                mask |= (1 << j)
        masks.append(mask)

    alpha = [0] * (N + 1)
    max_non_recovery = -1

    # Enumerate all subsets
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
    N = len(alpha) - 1

    def H(k):
        return sum(1.0 / j for j in range(1, k + 1))

    val = N * (H(N) - H(N - M - 1))
    for s in range(M + 1):
        if comb(N - 1, s) > 0:
            val -= alpha[s] / comb(N - 1, s)
    return val


# -------------------------------------------------------------
# Main
# -------------------------------------------------------------
if __name__ == "__main__":
    print("Exact recovery-set computation for AG(2,q)")
    p = int(input("Enter prime p: "))
    k = int(input("Enter power k: "))

    A = generate_AG2q(p, k)
    q = p**k
    print(f"\nGenerated AG(2,{q}) with {A.shape[1]} points and {A.shape[0]} lines.")

    alpha, M = compute_recovery_sets(A)

    print("\nRecovery-set statistics:")
    print(f"M (largest non-recovery size) = {M}")
    print("alpha(s) values:")
    for s,val in enumerate(alpha):
        if val > 0:
            print(f"  s = {s}: {val}")

    ET = expected_coverage_time(alpha, M)
    ET_frac = Fraction(ET).limit_denominator()
    a, b = ET_frac.numerator, ET_frac.denominator


print("\nExpected coverage time:")
print(f"  E[T] = {a}/{b} ≈ {float(ET):.6f}")

# -------------------------------------------------------------
# Optional CSV Saving
# -------------------------------------------------------------
save = input("\nSave results to CSV file? (y/n): ").strip().lower()

if save == "y":
    q = p**k
    filename = f"AG_q{q}.csv"

    with open(filename, "w") as f:
        # Header
        f.write("s,alpha_s\n")
        # α(s) table
        for s, val in enumerate(alpha):
            f.write(f"{s},{val}\n")
        # Expected time (exact + approx)
        f.write(f"E_T_fraction,{a}/{b}\n")
        f.write(f"E_T_approx,{float(ET):.12f}\n")

    print(f"\nCSV results saved to {filename}\n")
else:
    print()
