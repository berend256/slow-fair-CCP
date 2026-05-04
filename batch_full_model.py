
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
batch_full_model.py
===================

Compute full-model expected coverage times for the pairs (n, ell)
arising from AG(2,q) and/or PG(2,q), for any list of q-values chosen
by the user.

Writes all results to a single text file:
    full_model_results.txt

Requires:
    expected_time_full_model.py   (same directory)
"""

import datetime
from expected_time_full_model import expected_time_full


# ============================================================
# Utility: produce parameter pairs (n, ell)
# ============================================================

def AG_params(q):
    """Return (n, ell) for AG(2,q): n = q^2, ell = q."""
    return (q*q, q)

def PG_params(q):
    """Return (n, ell) for PG(2,q): n = q^2+q+1, ell = q+1."""
    return (q*q + q + 1, q+1)


# ============================================================
# Main batch driver
# ============================================================

def main():

    print("Batch computation of FULL-MODEL expected times (Pólya).")
    print("=======================================================")

    # Select which families to process
    mode = input(
        "Compute for (A) AG only, (P) PG only, or (B) Both?  [A/P/B]: "
    ).strip().upper()

    if mode not in {"A", "P", "B"}:
        print("Invalid choice. Aborting.")
        return

    # Select q-values
    print("\nEnter q-values (prime powers) separated by spaces.")
    print("Example: 2 3 4 5 7 8 9 11 13 17")
    q_list = input("q-values: ").strip().split()

    try:
        q_list = [int(q) for q in q_list]
    except ValueError:
        print("Invalid q-value detected. Aborting.")
        return

    # output file
    outfile = "full_model_results.txt"

    with open(outfile, "w") as f:

        f.write("Full-model expected times (Pólya inclusion–exclusion)\n")
        f.write("=====================================================\n")
        f.write(f"Computed on: {datetime.datetime.now()}\n\n")

        for q in q_list:

            if mode in {"A", "B"}:
                n, ell = AG_params(q)
                f.write(f"AG(2,{q}):  (n, ell) = ({n}, {ell})\n")
                res = expected_time_full(n, ell, max_digits=20)

                if res['fraction'] is not None:
                    f.write(f"  E[T_full] = {res['fraction']} ≈ {res['decimal']:.12f}\n")
                else:
                    f.write(f"  E[T_full] ≈ {res['decimal']:.12f}\n")
                    f.write(f"    (Exact denominator has {res['denominator_digits']} digits)\n")

                f.write("\n")

            if mode in {"P", "B"}:
                n, ell = PG_params(q)
                f.write(f"PG(2,{q}):  (n, ell) = ({n}, {ell})\n")
                res = expected_time_full(n, ell, max_digits=20)

                if res['fraction'] is not None:
                    f.write(f"  E[T_full] = {res['fraction']} ≈ {res['decimal']:.12f}\n")
                else:
                    f.write(f"  E[T_full] ≈ {res['decimal']:.12f}\n")
                    f.write(f"    (Exact denominator has {res['denominator_digits']} digits)\n")

                f.write("\n")

        f.write("=====================================================\n")
        f.write("End of computation.\n")

    print(f"\nResults written to {outfile}\n")


# ============================================================

if __name__ == "__main__":
    main()

