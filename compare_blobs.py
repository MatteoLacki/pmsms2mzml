#!/usr/bin/env python3
"""Decode and compare the first spectrum m/z binary blobs from two mzML files."""
import base64, zlib, struct, sys, re

def extract_first_mz_blob(path):
    with open(path) as f:
        text = f.read()
    # Find first m/z array binary
    idx = text.find('name="m/z array"')
    idx = text.find('<binary>', idx)
    end = text.find('</binary>', idx)
    b64 = text[idx+8:end]
    raw = zlib.decompress(base64.b64decode(b64))
    n = len(raw) // 4
    return struct.unpack(f'<{n}f', raw)

claude = extract_first_mz_blob(sys.argv[1])
stefan = extract_first_mz_blob(sys.argv[2])

print(f"Claude: {len(claude)} values, Stefan: {len(stefan)} values")
if len(claude) != len(stefan):
    print("LENGTH MISMATCH!")
    sys.exit(1)

diffs = 0
for i, (c, s) in enumerate(zip(claude, stefan)):
    if c != s:
        diffs += 1
        if diffs <= 20:
            print(f"  [{i:3d}] claude={c:.6f}  stefan={s:.6f}  diff={c-s:+.9f}")

print(f"\nTotal differences: {diffs} / {len(claude)}")
if diffs == 0:
    print("IDENTICAL!")
