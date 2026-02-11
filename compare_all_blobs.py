#!/usr/bin/env python3
"""Compare ALL m/z and intensity binary blobs between two mzML files."""
import base64, zlib, struct, re, sys

def extract_blobs(path):
    with open(path) as f:
        text = f.read()
    mz_blobs = []
    int_blobs = []
    for m in re.finditer(r'<binaryDataArrayList.*?</binaryDataArrayList>', text, re.DOTALL):
        block = m.group(0)
        binaries = re.findall(r'<binary>(.*?)</binary>', block)
        if len(binaries) == 2:
            mz_blobs.append(binaries[0])
            int_blobs.append(binaries[1])
    return mz_blobs, int_blobs

a_mz, a_int = extract_blobs(sys.argv[1])
b_mz, b_int = extract_blobs(sys.argv[2])

n = min(len(a_mz), len(b_mz))
mz_diffs = sum(1 for i in range(n) if a_mz[i] != b_mz[i])
int_diffs = sum(1 for i in range(n) if a_int[i] != b_int[i])

print(f"Spectra: A={len(a_mz)}, B={len(b_mz)}")
print(f"m/z blob mismatches:       {mz_diffs} / {n}")
print(f"intensity blob mismatches: {int_diffs} / {n}")
