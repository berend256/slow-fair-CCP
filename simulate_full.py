"""
Simulation of the full (n, ell) multi-coupon collector.

Model:
------
There are n coupon types.
Each draw selects ell DISTINCT coupons, uniformly at random, from all
binom(n, ell) possible subsets.

We do not track WHICH coupons were collected, only HOW MANY distinct
coupon types have been collected so far.

If at some stage n' coupon types have already been collected, then the
number of new coupon types obtained in the next draw is distributed as

    Hypergeometric(
        population size = n,
        number of successes = n - n',   (unseen coupons)
        sample size = ell
    ).

The process stops when all n coupon types have been collected.

Output:
-------
The program runs a Monte-Carlo simulation of this process and reports
the mean and standard deviation of the collection time.

For small (n, ell), where it is computationally feasible, the program
also computes the exact expected collection time using Pólya’s
inclusion–exclusion formula and reports it.

Precision:
----------
The number of Monte-Carlo iterations is set to ITER = 1000 by default.
If higher precision is required, increase ITER (at the cost of runtime).
"""

import numpy as np
import math
import statistics


# ==================================================
# User-adjustable parameter
# ==================================================
ITER = 1000   # Increase this value for higher precision


def polya_expectation(n, ell):
    """
    Compute E[T_full(n, ell)] using Pólya's inclusion–exclusion formula.

    WARNING:
    This computation costs Theta(n * ell) arithmetic operations and is
    feasible only for relatively small n (roughly n <= 300–500).
    """
    from math import comb

    denom = comb(n, ell)
    total = 0.0

    for k in range(1, n + 1):
        if n - k >= ell:
            miss_prob = comb(n - k, ell) / denom
        else:
            miss_prob = 0.0

        total += ((-1) ** (k - 1)) * comb(n, k) / (1.0 - miss_prob)

    return total


def simulate_full_multicoupon(n, ell, num_sims=ITER):
    """
    Monte-Carlo simulation of the full (n, ell) multi-coupon collector
    using hypergeometric transitions.
    """
    rng = np.random.default_rng()
    times = []

    for _ in range(num_sims):
        collected = 0
        t = 0

        while collected < n:
            t += 1
            new = rng.hypergeometric(
                ngood=n - collected,
                nbad=collected,
                nsample=ell
            )
            collected += new

        times.append(t)

    return {
        "mean": statistics.mean(times),
        "std": statistics.pstdev(times),
    }


# ==================================================
# Main program
# ==================================================
if __name__ == "__main__":

    n = int(input("Enter n (number of coupon types): "))
    ell = int(input("Enter ell (coupons per draw): "))

    if not (1 <= ell <= n):
        raise ValueError("ell must satisfy 1 <= ell <= n")

    print("\nRunning Monte-Carlo simulation...")
    print(f"Number of iterations: {ITER}")

    sim = simulate_full_multicoupon(n, ell)

    print("\nMonte-Carlo results:")
    print(f"Mean collection time: {sim['mean']:.6f}")
    print(f"Std deviation:        {sim['std']:.6f}")

    # Decide whether to attempt Pólya
    POLYA_THRESHOLD = 400

    print("\nExact calculation:")
    if n <= POLYA_THRESHOLD:
        print("Computing exact expected time using Pólya's formula...")
        exact = polya_expectation(n, ell)
        print(f"The exact expected time, by Pólya's formula, is {exact:.6f}")
    else:
        print(
            "An exact calculation using Pólya's formula would take too long "
            "for these values of (n, ell)."
        )
