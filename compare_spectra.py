#!/usr/bin/env python3
"""Compare defaultArrayLength of all spectra between two mzML files."""
import re, sys

def get_array_lengths(path):
    with open(path) as f:
        return [int(m) for m in re.findall(r'defaultArrayLength="(\d+)"', f.read())]

a = get_array_lengths(sys.argv[1])
b = get_array_lengths(sys.argv[2])

print(f"File A: {len(a)} spectra, File B: {len(b)} spectra")
if len(a) != len(b):
    print(f"COUNT MISMATCH: {len(a)} vs {len(b)}")

n = min(len(a), len(b))
diffs = 0
for i in range(n):
    if a[i] != b[i]:
        diffs += 1
        if diffs <= 20:
            print(f"  [{i}] A={a[i]}  B={b[i]}")

# Check for spectra only in the longer file
if len(a) > len(b):
    print(f"  A has {len(a)-len(b)} extra spectra at the end")
elif len(b) > len(a):
    print(f"  B has {len(b)-len(a)} extra spectra at the end")

print(f"\nSize mismatches: {diffs} / {n}")
if diffs == 0 and len(a) == len(b):
    print("ALL IDENTICAL!")
