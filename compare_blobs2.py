#!/usr/bin/env python3
"""Compare raw (unrounded) claude m/z vs stefan m/z, and show 3-decimal text."""
import base64, zlib, struct, sys

def extract_first_mz_blob(path):
    with open(path) as f:
        text = f.read()
    idx = text.find('name="m/z array"')
    idx = text.find('<binary>', idx)
    end = text.find('</binary>', idx)
    b64 = text[idx+8:end]
    raw = zlib.decompress(base64.b64decode(b64))
    n = len(raw) // 4
    return struct.unpack(f'<{n}f', raw)

raw = extract_first_mz_blob(sys.argv[1])    # no --decimals (raw floats)
stefan = extract_first_mz_blob(sys.argv[2])  # stefan

print(f"{'idx':>4} {'raw_float':>14} {'raw_3dec':>10} {'stefan_3dec':>12} {'match':>6}")
print("-" * 52)
diffs = 0
for i, (r, s) in enumerate(zip(raw, stefan)):
    r3 = f"{r:.3f}"
    s3 = f"{s:.3f}"
    match = "OK" if r3 == s3 else "DIFF"
    if match == "DIFF":
        diffs += 1
    if i < 30 or match == "DIFF" and diffs <= 30:
        print(f"{i:4d} {r:14.6f} {r3:>10} {s3:>12} {match:>6}")

print(f"\nTotal 3-decimal mismatches: {diffs} / {len(raw)}")
