#!/usr/bin/env python3
"""Investigate where and why ordering diverges between two mzML files."""
import re, sys

def parse_spectra(path):
    with open(path) as f:
        text = f.read()
    spectra = []
    for m in re.finditer(r'<spectrum index="(\d+)".*?</spectrum>', text, re.DOTALL):
        block = m.group(0)
        def cv(name):
            p = re.search(rf'name="{re.escape(name)}" value="([^"]*)"', block)
            return p.group(1) if p else None
        spectra.append({
            'prec_mz': cv('selected ion m/z'),
            'charge': cv('charge state'),
            'rt': cv('scan start time'),
            'iim': cv('inverse reduced ion mobility'),
        })
    return spectra

a = parse_spectra(sys.argv[1])
b = parse_spectra(sys.argv[2])
n = min(len(a), len(b))

# Find first divergence point
first_diff = None
for i in range(n):
    if a[i]['charge'] != b[i]['charge']:
        first_diff = i
        break

print(f"First charge mismatch at index {first_diff}")
print(f"\nContext around divergence (idx {max(0,first_diff-3)} to {first_diff+8}):")
print(f"{'idx':>4} {'A_mz':>12} {'A_ch':>4} {'A_rt':>12} {'A_iim':>8} | {'B_mz':>12} {'B_ch':>4} {'B_rt':>12} {'B_iim':>8}")
print("-" * 100)
for i in range(max(0,first_diff-3), min(n, first_diff+8)):
    print(f"{i:4d} {a[i]['prec_mz']:>12} {a[i]['charge']:>4} {a[i]['rt']:>12} {a[i]['iim']:>8} | {b[i]['prec_mz']:>12} {b[i]['charge']:>4} {b[i]['rt']:>12} {b[i]['iim']:>8}")

# Check: do spectra around divergence match if we look at (mz, charge, iim) tuples?
print(f"\n--- Trying to find A[{first_diff}] in B nearby ---")
a_key = (a[first_diff]['prec_mz'][:6], a[first_diff]['charge'], a[first_diff]['iim'])
for j in range(max(0,first_diff-5), min(n, first_diff+10)):
    b_key = (b[j]['prec_mz'][:6], b[j]['charge'], b[j]['iim'])
    if a_key[0] == b_key[0]:
        print(f"  A[{first_diff}] mz~= B[{j}]: A=({a_key}) B=({b_key})")

# Count how many groups of same prec_mz exist around divergence
print(f"\n--- Groups of same prec_mz around divergence ---")
from collections import Counter
for label, spec in [("A (claude)", a), ("B (stefan)", b)]:
    mzs = [s['prec_mz'][:6] for s in spec[max(0,first_diff-2):first_diff+6]]
    groups = Counter(mzs)
    multi = {k:v for k,v in groups.items() if v > 1}
    if multi:
        print(f"  {label}: repeated m/z groups: {multi}")
        for k in multi:
            indices = [i+max(0,first_diff-2) for i,s in enumerate(spec[max(0,first_diff-2):first_diff+6]) if s['prec_mz'][:6]==k]
            for idx in indices:
                print(f"    [{idx}] mz={spec[idx]['prec_mz']} ch={spec[idx]['charge']} iim={spec[idx]['iim']}")
