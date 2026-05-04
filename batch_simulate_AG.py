
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
batch_simulate_AG.py
====================

Run Monte‑Carlo simulations for AG(2,q) for all prime powers q ≤ 17,
with user‑specified number of iterations per q.
"""

import datetime
from simulate_AG import simulate_coupon

# prime powers ≤ 17
prime_powers = []
primes = [2, 3, 5, 7, 11, 13, 17]

for p in primes:
    k = 1
    while p**k <= 17:
        prime_powers.append((p, k, p**k))
        k += 1

prime_powers.sort(key=lambda x: x[2])

# ==============================================================
# Ask user for number of iterations per q
# ==============================================================

print("Prime powers to be tested:", [q for (_,_,q) in prime_powers])
num_sims = int(input("How many simulations per q? (e.g., 1000, 100000, 1000000): "))

# ==============================================================
# Prepare output file
# ==============================================================

outfile = "AG_simulation_results.txt"

with open(outfile, "w") as f:

    f.write("Batch Monte‑Carlo Simulation for AG(2,q)\n")
    f.write("======================================\n")
    f.write(f"Start time: {datetime.datetime.now()}\n\n")

    for (p, k, q) in prime_powers:

        print(f"\nRunning AG(2,{q}) with {num_sims} simulations ...")

        results = simulate_coupon(p, k, num_sims=num_sims, verbose=False)

        mean = results["mean"]
        std = results["std"]
        trunc = results["trunc"]

        f.write(f"q = {q} (p={p}, k={k})\n")
        f.write(f"  Mean: {mean:.8f}\n")
        f.write(f"  Std: {std:.8f}\n")
        f.write(f"  Truncated: {trunc}/{num_sims}\n")
        f.write("\n")
        f.flush()

    f.write(f"End time: {datetime.datetime.now()}\n")
