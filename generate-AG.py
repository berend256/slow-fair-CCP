#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AG(2,q) Generator over GF(p^k) using galois
===========================================

OVERVIEW
--------
Generates the affine plane AG(2,q), where q = p^k is a prime power.
Each line is printed as a list of indices of points lying on it.

POINTS
------
All pairs (x, y) in GF(q)^2. Number of points = q^2.

LINES
-----
- Non-vertical: y = a x + b, a, b ∈ GF(q)
- Vertical: x = c, c ∈ GF(q)
Number of lines = q^2 + q

INCIDENCE
---------
Point lies on line if it satisfies the equation in GF(q).
"""

import numpy as np
import galois

def generate_ag2q(p, k):
    q = p**k
    GF = galois.GF(q)

    elems = GF.elements
    zero = GF(0)
    one = GF(1)

    # -------- Points --------
    points = []
    for x in elems:
        for y in elems:
            points.append(GF([x, y]))
    n_points = len(points)
    print(f"\nGenerating AG(2,{q}) over GF({p}^{k}) with {n_points} points")

    # -------- Lines --------
    lines = []

    # non-vertical: y = a x + b
    for a in elems:
        for b in elems:
            lines.append(("non-vertical", a, b))

    # vertical: x = c
    for c in elems:
        lines.append(("vertical", c))

    # -------- Print lines --------
    for idx, line in enumerate(lines, 1):
        pts_on_line = []
        if line[0] == "non-vertical":
            a, b = line[1], line[2]
            for i, pt in enumerate(points, 1):
                if pt[1] == a*pt[0] + b:
                    pts_on_line.append(i)
        else:
            c = line[1]
            for i, pt in enumerate(points, 1):
                if pt[0] == c:
                    pts_on_line.append(i)
        print(f"Line {idx}: {' '.join(map(str, pts_on_line))}")


# ==================================================
# Main
# ==================================================
if __name__ == "__main__":
    p = int(input("Enter prime p: "))
    k = int(input("Enter power k: "))
    generate_ag2q(p, k)
