#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PG(2,q) Generator over GF(p^k) using the galois library
=======================================================

OVERVIEW
--------
Generates the finite projective plane PG(2,q), where q = p^k is a
prime power, using the `galois` library for finite field arithmetic.

KEY FEATURES
------------
- Correctly generates all points and lines of PG(2,q)
- Deduplicates points and lines using Python sets
- Computes incidence correctly via GF dot products
- Output lists points on each line, exactly q+1 per line
- Works for any prime power q supported by galois
"""

import numpy as np
import galois
from itertools import product

# ==================================================
# Projective plane generator
# ==================================================
def generate_pg2q(p, k):
    # Construct GF(q)
    q = p**k
    GF = galois.GF(q)
    elems = GF.elements
    zero = GF(0)
    one = GF(1)

    n = q*q + q + 1
    print(f"\nGenerating PG(2,{q}) over GF({p}^{k}) with {n} points")

    # -------- Points --------
    points_set = set()
    for x, y, z in product(elems, repeat=3):
        if not (x == zero and y == zero and z == zero):
            if x != 0:
                inv = x**-1
                pt = (int(one), int(y*inv), int(z*inv))
            elif y != 0:
                inv = y**-1
                pt = (0, int(one), int(z*inv))
            else:
                pt = (0, 0, 1)
            points_set.add(pt)

    points = [GF(list(p)) for p in sorted(points_set)]

    # -------- Lines --------
    lines_set = set()
    for a, b, c in product(elems, repeat=3):
        if not (a == zero and b == zero and c == zero):
            if a != 0:
                inv = a**-1
                ln = (int(one), int(b*inv), int(c*inv))
            elif b != 0:
                inv = b**-1
                ln = (0, int(one), int(c*inv))
            else:
                ln = (0, 0, 1)
            lines_set.add(ln)

    lines = [GF(list(l)) for l in sorted(lines_set)]
    assert len(lines) == len(points)

    # -------- Incidence --------
    for idx, line in enumerate(lines, 1):
        pts_on_line = []
        for i, point in enumerate(points, 1):
            if line @ point == 0:   # dot product in GF(q)
                pts_on_line.append(i)
        print(f"Line {idx}: {' '.join(map(str, pts_on_line))}")


# ==================================================
# Main
# ==================================================
if __name__ == "__main__":
    p = int(input("Enter prime p: "))
    k = int(input("Enter power k: "))
    generate_pg2q(p, k)
