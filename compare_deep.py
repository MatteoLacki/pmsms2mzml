#!/usr/bin/env python3
"""Compare ordering: mmappet precursors vs MGF vs stefan mzML."""
import re, sys, struct

def parse_mzml(path):
    with open(path) as f:
        text = f.read()
    spectra = []
    for m in re.finditer(r'<spectrum index="(\d+)".*?</spectrum>', text, re.DOTALL):
        block = m.group(0)
        def cv(name):
            p = re.search(rf'name="{re.escape(name)}" value="([^"]*)"', block)
            return p.group(1) if p else None
        spectra.append((cv('selected ion m/z'), cv('charge state'), cv('inverse reduced ion mobility')))
    return spectra

def parse_mgf(path):
    spectra = []
    with open(path) as f:
        mz = ch = iim = None
        for line in f:
            line = line.strip()
            if line.startswith('PEPMASS='):
                mz = line.split('=')[1].split()[0]
            elif line.startswith('CHARGE='):
                ch = line.split('=')[1].rstrip('+').rstrip('-')
            elif line.startswith('ION_MOBILITY=') or line.startswith('INVERSE_ION_MOBILITY='):
                iim = line.split('=')[1]
            elif line == 'END IONS':
                spectra.append((mz, ch, iim))
                mz = ch = iim = None
    return spectra

def read_mmappet_col(path, dtype, count):
    import numpy as np
    return np.fromfile(path, dtype=dtype, count=count)

# Parse stefan and claude mzML
stefan = parse_mzml('reduced_input/stefan.mzML')
claude = parse_mzml('reduced_input/claude.mzML')

# Parse MGF
mgf = parse_mgf('reduced_input/mgf.mgf')

print(f"Claude mzML: {len(claude)} spectra")
print(f"Stefan mzML: {len(stefan)} spectra")
print(f"MGF: {len(mgf)} spectra")

# Check: does stefan match MGF order?
mgf_vs_stefan = 0
for i in range(min(len(mgf), len(stefan))):
    if mgf[i][0] != stefan[i][0] or mgf[i][1] != stefan[i][1]:
        mgf_vs_stefan += 1
        if mgf_vs_stefan <= 3:
            print(f"  MGF vs Stefan mismatch [{i}]: MGF=({mgf[i]}) Stefan=({stefan[i]})")
print(f"MGF vs Stefan mismatches: {mgf_vs_stefan}")

# Show MGF order around divergence point (254)
print(f"\nMGF around index 254:")
for i in range(251, min(262, len(mgf))):
    print(f"  [{i}] mz={mgf[i][0]} ch={mgf[i][1]} iim={mgf[i][2]}")

print(f"\nClaude around index 254:")
for i in range(251, min(262, len(claude))):
    print(f"  [{i}] mz={claude[i][0]} ch={claude[i][1]} iim={claude[i][2]}")
