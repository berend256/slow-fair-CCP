#!/usr/bin/env python3

# =====================================================================
# FULL vs STAR comparison on Ω = binom([m], k)
# =====================================================================
#
# PURPOSE OF THIS PROGRAM
# -----------------------
#
# This program numerically compares two FAIR coupon‑collection mechanisms
# defined on the SAME coupon space:
#
#       Ω = binom([m], k)
#
# where each coupon is a k‑subset of [m].
#
# The two mechanisms compared are:
#
#   (1) FULL  — simulated by Monte‑Carlo (exactly, via hypergeometric RVs)
#   (2) STAR  — not simulated; exact expectation is used
#
# The goal is to compare the expected coverage times of FULL and STAR
# for FIXED parameters (m, k).
#
#
# ---------------------------------------------------------------------
# COUPON SPACE AND PARAMETERS
# ---------------------------------------------------------------------
#
# Coupon space:
#     Ω = binom([m], k)
#     n = |Ω| = binom(m, k)
#
# Block size:
#     ℓ = binom(m-1, k-1)
#
# Both mechanisms are FAIR mechanisms with these parameters:
# every draw reveals exactly ℓ coupons, and all coupons are symmetric.
#
#
# ---------------------------------------------------------------------
# MECHANISM 1: STAR (NOT SIMULATED)
# ---------------------------------------------------------------------
#
# Definition of STAR:
#
#   In each draw:
#       • choose i ∈ [m] uniformly
#       • reveal all k‑subsets A ⊆ [m] such that i ∈ A
#
# The revealed set is a "star" of size
#       |F_i| = binom(m-1, k-1) = ℓ.
#
# Properties:
#   • Each coupon A ∈ Ω lies in exactly k stars
#   • The mechanism is fair
#
# Crucially:
#   The expected coverage time of STAR is known EXACTLY:
#
#       E[T_star(m,k)] = m ( H_m − H_{k−1} )
#
# where H_n is the nth harmonic number.
#
# Therefore:
#   → STAR is NOT simulated
#   → We simply compute and print its exact expectation
#
#
# ---------------------------------------------------------------------
# MECHANISM 2: FULL (SIMULATED)
# ---------------------------------------------------------------------
#
# Definition of FULL:
#
#   In each draw:
#       • choose uniformly at random an ℓ‑subset of Ω
#         (i.e., ℓ distinct k‑subsets of [m])
#
# This is the “maximally mixed” fair mechanism on Ω.
#
#
# CORRECT STATE COMPRESSION FOR FULL
# ---------------------------------
#
# We track:
#
#     U_t = number of unseen coupons in Ω after t draws
#
# Initially:
#     U_0 = n = binom(m, k)
#
#
# ONE‑STEP TRANSITION (EXACT)
# --------------------------
#
# Given U_t = u, a single FULL draw selects ℓ coupons uniformly from Ω.
# Among the n coupons, exactly u are unseen.
#
# Hence:
#
#     X_t ~ Hypergeom(n, u, ℓ)
#
# where X_t is the number of newly discovered coupons.
#
# The state update is:
#
#     U_{t+1} = U_t − X_t
#
# The process stops when U_t = 0.
#
#
# IMPLEMENTATION NOTE
# -------------------
#
# This program uses NumPy’s hypergeometric random‑variable generator.
# NumPy internally switches between exact, rejection, and normal‑type
# methods depending on the parameter regime.
#
# IMPORTANT RESTRICTION:
#   This program works ONLY when n = binom(m,k) and ℓ = binom(m-1,k-1)
#   fit into machine‑size integers (roughly n ≤ 10^18).
#
# Examples that WORK:
#   • small and moderate (m, k)
#   • sanity checks against exact formulas
#
# Examples that DO NOT work:
#   • m = 1000, k = 500  (n ~ 10^299)
#
# In such regimes, asymptotic analysis must replace simulation.
#
#
# ---------------------------------------------------------------------
# OUTPUT
# ---------------------------------------------------------------------
#
# The program prints:
#
#   • Monte‑Carlo mean and standard deviation of T_full
#   • Exact expectation E[T_star]
#
# A correct run should show:
#
#   • non‑zero variance for FULL
#   • sensible comparison between FULL and STAR
#
# =====================================================================

import statistics
from math import comb
import numpy as np

# ---------------------------------------------------------------------
# Global parameters
# ---------------------------------------------------------------------

# Number of Monte‑Carlo iterations for FULL
ITER = 100000

# NumPy random number generator
rng = np.random.default_rng()


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def harmonic(n):
    """Return the nth harmonic number H_n."""
    return sum(1.0 / i for i in range(1, n + 1))


def star_expectation(m, k):
    """
    Exact expected coverage time of the STAR mechanism.
    """
    return m * (harmonic(m) - harmonic(k - 1))


# ---------------------------------------------------------------------
# FULL simulation (exact hypergeometric)
# ---------------------------------------------------------------------

def simulate_full_once(m, k):
    """
    Simulate ONE run of the FULL mechanism.

    State variable:
        u = number of unseen coupons in Ω = binom([m],k)

    Transition:
        X ~ Hypergeom(n, u, ℓ)
        u <- u - X
    """
    n = comb(m, k)
    ell = comb(m - 1, k - 1)

    u = n
    t = 0

    while u > 0:
        t += 1
        x = rng.hypergeometric(u, n - u, ell)
        u -= x

    return t


def simulate_full(m, k, iters=ITER):
    """
    Monte‑Carlo simulation of the FULL mechanism.
    """
    times = []
    for i in range(iters):
        times.append(simulate_full_once(m, k))
        if (i + 1) % max(1, iters // 5) == 0:
            print(f"[INFO] completed {i+1}/{iters} runs")
    return statistics.mean(times), statistics.pstdev(times)


# ---------------------------------------------------------------------
# Main program
# ---------------------------------------------------------------------

if __name__ == "__main__":

    print("\nFULL vs STAR comparison on Ω = binom([m], k)\n")

    # Read input
    m, k = map(int, input("Enter m k: ").split())

    n = comb(m, k)
    ell = comb(m - 1, k - 1)

    print("\nParameters:")
    print(f"  m = {m}")
    print(f"  k = {k}")
    print(f"  n = binom(m,k) = {n}")
    print(f"  ℓ = binom(m-1,k-1) = {ell}")
    print(f"\nMonte‑Carlo iterations: {ITER}\n")

    # FULL simulation
    print("[INFO] Simulating FULL (exact hypergeometric)...")
    full_mean, full_std = simulate_full(m, k)

    print("\nFULL model (Monte‑Carlo):")
    print(f"  Mean completion time = {full_mean:.6f}")
    print(f"  Std deviation        = {full_std:.6f}")

    # STAR (exact)
    print("\nSTAR mechanism (exact):")
    print(f"  Expected time = {star_expectation(m, k):.6f}")
