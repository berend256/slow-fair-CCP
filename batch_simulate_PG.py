
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
batch_simulate_PG.py
====================

Run Monte‑Carlo simulations for PG(2,q) for any list of prime powers q,
with user‑specified number of iterations per q.

Results are written to a single text file:
    PG_simulation_results.txt

Requires:
    simulate_PG.py   (same directory)
"""

import datetime
from simulate_PG import simulate_coupon


# ============================================================
# Main batch driver
# ============================================================

def main():

    print("Batch Monte‑Carlo Simulation for PG(2,q)")
    print("========================================\n")

    # ----------------------------------------------------------
    # Ask user for q-values
    # ----------------------------------------------------------
    print("Enter q-values (prime powers), separated by spaces.")
    print("Example: 2 3 4 5 7 8 9 11 13")
    q_list = input("q-values: ").strip().split()

    try:
        q_list = [int(q) for q in q_list]
    except ValueError:
        print("Invalid q-value. Aborting.")
        return

    # ----------------------------------------------------------
    # Ask user for number of simulations
    # ----------------------------------------------------------
    num_sims = int(input("How many simulations per q? (e.g., 1000 or 1000000): "))
    print()

    # ----------------------------------------------------------
    # Output file (always overwritten)
    # ----------------------------------------------------------
    outfile = "PG_simulation_results.txt"

    with open(outfile, "w") as f:

        f.write("Monte‑Carlo Simulation for PG(2,q)\n")
        f.write("=================================\n")
        f.write(f"Start time: {datetime.datetime.now()}\n\n")

        # ------------------------------------------------------
        # Loop over q-values
        # ------------------------------------------------------
        for q in q_list:
            print(f"Running PG(2,{q}) with {num_sims} simulations...")

            # Determine p and k such that q = p^k
            # (We trust the user to enter prime powers, but we compute k.)
            # Naively: try all primes <= q
            found = False
            for p in range(2, q+1):
                k = 1
                while p**k <= q:
                    if p**k == q:
                        found = True
                        break
                    k += 1
                if found:
                    break

            if not found:
                print(f"  Warning: q={q} is not a prime power. Skipping.\n")
                f.write(f"q = {q}: not a prime power (skipped)\n\n")
                continue

            # --------------------------------------------------
            # Run simulation silently
            # --------------------------------------------------
            results = simulate_coupon(p, k, num_sims=num_sims, verbose=False)

            mean = results["mean"]
            std  = results["std"]
            trunc = results["trunc"]

            # --------------------------------------------------
            # Write block to file
            # --------------------------------------------------
            f.write(f"PG(2,{q}): (p={p}, k={k}),  n = {q*q + q + 1}\n")
            f.write(f"  Mean: {mean:.8f}\n")
            f.write(f"  Std:  {std:.8f}\n")
            f.write(f"  Truncated: {trunc}/{num_sims}\n\n")
            f.flush()

        f.write(f"End time: {datetime.datetime.now()}\n")

    print(f"\nResults written to {outfile}\n")


# ============================================================

if __name__ == "__main__":
    main()

