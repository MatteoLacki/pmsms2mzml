#!/usr/bin/env python3
"""Compare two mzML files by scan number (order-independent)."""
import re, sys

def extract_spectra(path):
    with open(path) as f:
        text = f.read()
    spectra = {}
    for m in re.finditer(r'<spectrum [^>]*index="(\d+)"[^>]*>.*?</spectrum>', text, re.DOTALL):
        block = m.group(0)
        nid = re.search(r'NativeID:&quot;index=(\d+)&quot;', block)
        if nid:
            scan = int(nid.group(1))
        else:
            scan = int(m.group(1))
        # Remove title line for comparison (contains filename)
        cleaned = re.sub(r'<cvParam[^>]*name="spectrum title"[^>]*/>', '', block)
        spectra[scan] = cleaned
    return spectra

a = extract_spectra(sys.argv[1])
b = extract_spectra(sys.argv[2])

print(f"File A: {len(a)} spectra, File B: {len(b)} spectra")
only_a = set(a) - set(b)
only_b = set(b) - set(a)
common = set(a) & set(b)
if only_a: print(f"Only in A: {sorted(only_a)[:10]}")
if only_b: print(f"Only in B: {sorted(only_b)[:10]}")
print(f"Common scans: {len(common)}")

diffs = 0
for scan in sorted(common):
    if a[scan] != b[scan]:
        diffs += 1
        if diffs <= 3:
            la, lb = a[scan].split('\n'), b[scan].split('\n')
            for i, (x, y) in enumerate(zip(la, lb)):
                if x != y:
                    print(f"  Scan {scan} line {i}: A={x.strip()[:80]}  B={y.strip()[:80]}")
                    break

print(f"Mismatches: {diffs} / {len(common)}")
