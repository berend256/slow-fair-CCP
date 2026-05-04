#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
expected_time_full_model.py
===========================

Compute the expected coverage time for the FULL uniform batch
coupon collector with parameters (n, ell), using Pólya’s
inclusion–exclusion formula.

FULL MODEL:
Each step samples an ℓ-subset of [n] uniformly at random.

Correct Pólya formula:
    E[T_full(n, ell)]
    = sum_{s=1}^{n}
        (-1)^{s+1} * C(n, s)
        / (1 - C(n-s, ell) / C(n, ell))

Features:
    • exact rational computation (Fraction)
    • decimal approximation always available
    • rational displayed only if denominator is small enough
"""

from math import comb
from fractions import Fraction
import sys

# allow very large integer handling
sys.set_int_max_str_digits(1_000_000)


# ============================================================
# FULL MODEL EXPECTATION (PÓLYA)
# ============================================================

def expected_time_full(n, ell, max_digits=20):
    """
    Compute E[T_full(n, ell)] for the FULL uniform batch coupon collector.

    Returns dict with keys:
        'n'
        'ell'
        'fraction'            (Fraction or None)
        'decimal'             (float)
        'denominator_digits'  (int)
        'used_exact'          (bool)
    """

    C_n_ell = Fraction(comb(n, ell), 1)
    E = Fraction(0, 1)

    for s in range(1, n + 1):

        # probability that a batch avoids a fixed s-set
        if n - s >= ell:
            p_no = Fraction(comb(n - s, ell), 1) / C_n_ell
        else:
            p_no = Fraction(0, 1)

        denom = 1 - p_no
        term = Fraction(comb(n, s), 1) / denom

        if s % 2 == 1:
            E += term
        else:
            E -= term

    # count denominator digits safely
    denom_digits = int(E.denominator.bit_length() / 3.3219280948873626) + 1
    used_exact = denom_digits <= max_digits

    return {
        "n": n,
        "ell": ell,
        "fraction": E if used_exact else None,
        "decimal": float(E),
        "denominator_digits": denom_digits,
        "used_exact": used_exact
    }


# ============================================================
# STANDALONE EXECUTION
# ============================================================

if __name__ == "__main__":

    print("Compute full-model expected coverage time (Pólya).")

    n = int(input("Enter n (number of coupons): "))
    ell = int(input("Enter ell (batch size): "))

    res = expected_time_full(n, ell, max_digits=1000)

    print("\n--- Full Model Expectation ---")
    print(f"(n, ell) = ({res['n']}, {res['ell']})")

    if res["fraction"] is not None:
        print(f"Expected time = {res['fraction']} ≈ {res['decimal']:.12f}")
    else:
        print(f"Expected time ≈ {res['decimal']:.12f}")
        print(
            f"(Exact fraction omitted: denominator has "
            f"{res['denominator_digits']} digits)"
        )
