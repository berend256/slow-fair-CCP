#!/usr/bin/env python3

# =====================================================================
# FULL vs STAR comparison on Ω = binom([m], k)
# =====================================================================
#
# PURPOSE
# -------
#
# This program compares the FULL and STAR coupon-collection mechanisms
# on the coupon space Ω = binom([m], k).
#
# The FULL mechanism draws ℓ = binom(m-1, k-1) coupons uniformly from Ω
# at each step. The STAR mechanism has a closed-form expectation.
#
#
# NUMERICAL DIFFICULTY
# --------------------
#
# For large m and k, the parameters
#
#     n = binom(m, k),   ℓ = binom(m-1, k-1)
#
# can be astronomically large. Exact hypergeometric sampling is only
# possible when the RNG parameters fit into machine integers.
#
# Therefore, the program automatically chooses one of the following:
#
#   (1) EXACT hypergeometric transition
#   (2) NORMAL approximation with finite-population correction
#   (3) POISSON approximation
#   (4) REFUSAL, with explanation
#
# depending on the asymptotic regime.
#
#
# ITERATIONS
# ----------
#
# ITER is fixed at 1000 by default.
# Increase it for higher precision (slower runtime).
#
# =====================================================================

import numpy as np
from math import comb, sqrt
import statistics


# ==================================================
# Monte-Carlo iterations (change if desired)
# ==================================================

ITER = 1000   # Increase for higher precision


# ==================================================
# Harmonic numbers (for STAR expectation)
# ==================================================

def harmonic(n):
    return sum(1.0 / i for i in range(1, n + 1))


def star_expectation(m, k):
    return m * (harmonic(m) - harmonic(k - 1))


# ==================================================
# One FULL run (adaptive)
# ==================================================

def simulate_full_once(m, k, rng):

    c = m - k
    n = comb(m, k)
    ell = comb(m - 1, k - 1)
    a = n - ell

    u = n
    t = 0

    # Decide regime once
    exact_possible = max(ell, a) < 2**63
    p = ell / n

    use_poisson = (u * p < 10)
    use_normal = (u > 50 and ell > 50 and a > 50)

    if not exact_possible and not use_normal and not use_poisson:
        raise RuntimeError(
            "FULL simulation refused: parameters too large and "
            "no justified approximation applies."
        )

    while u > 0:
        t += 1

        if exact_possible:
            u = rng.hypergeometric(u, n - u, a)

        elif use_poisson:
            lam = u * p
            u = max(0, u - rng.poisson(lam))

        else:
            mu = u * p
            var = u * p * (1 - p) * (n - ell) / (n - 1)
            x = int(round(rng.normal(mu, sqrt(var))))
            x = max(0, min(u, x))
            u -= x

    return t


# ==================================================
# Monte-Carlo FULL simulation
# ==================================================

def simulate_full(m, k, iters=ITER):

    rng = np.random.default_rng()
    times = []

    for i in range(iters):
        times.append(simulate_full_once(m, k, rng))
        if (i + 1) % max(1, iters // 5) == 0:
            print(f"[INFO] completed {i+1}/{iters} runs")

    return statistics.mean(times), statistics.pstdev(times)


# ==================================================
# Main program
# ==================================================

if __name__ == "__main__":

    print("\nFULL vs STAR comparison on Ω = binom([m], k)\n")

    m = int(input("Enter m: "))
    k = int(input("Enter k: "))

    n = comb(m, k)
    ell = comb(m - 1, k - 1)

    print("\nParameters:")
    print(f"  m = {m}")
    print(f"  k = {k}")
    print(f"  n = binom(m, k)")
    print(f"  ℓ = binom(m-1, k-1)")
    print(f"\nMonte-Carlo iterations: {ITER}\n")

    print("[INFO] Simulating FULL model...")

    try:
        full_mean, full_std = simulate_full(m, k)
        print("\nFULL model:")
        print(f"  Mean completion time = {full_mean:.6f}")
        print(f"  Std deviation        = {full_std:.6f}")
    except RuntimeError as e:
        print("\nFULL model:")
        print("  Simulation not performed.")
        print(" ", e)

    print("\nSTAR mechanism (exact):")
    print(f"  Expected time = {star_expectation(m, k):.6f}")
