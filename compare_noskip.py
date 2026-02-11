#!/usr/bin/env python3
"""Check: if we DON'T skip empty spectra, does mmappet order match MGF?"""
import numpy as np

base = 'reduced_input/pmsms.mmappet'
prec_dir = f'{base}/filtered_precursors_with_nontrivial_ms2.mmappet'

prec_mz = np.fromfile(f'{prec_dir}/6.bin', dtype=np.float64)
prec_charge = np.fromfile(f'{prec_dir}/8.bin', dtype=np.uint8)

def parse_mgf(path):
    spectra = []
    with open(path) as f:
        mz = ch = None
        for line in f:
            line = line.strip()
            if line.startswith('PEPMASS='):
                mz = float(line.split('=')[1].split()[0])
            elif line.startswith('CHARGE='):
                ch = int(line.split('=')[1].rstrip('+').rstrip('-'))
            elif line == 'END IONS':
                spectra.append((mz, ch))
                mz = ch = None
    return spectra

mgf = parse_mgf('reduced_input/mgf.mgf')

n = min(len(prec_mz), len(mgf))
diffs = 0
for i in range(n):
    mm_mz3 = f"{prec_mz[i]:.3f}"
    mgf_mz3 = f"{mgf[i][0]:.3f}"
    if mm_mz3 != mgf_mz3 or prec_charge[i] != mgf[i][1]:
        diffs += 1
        if diffs <= 5:
            print(f"  [{i}] mmappet=({prec_mz[i]:.3f}, ch={prec_charge[i]}) mgf=({mgf[i][0]:.3f}, ch={mgf[i][1]})")

print(f"Order mismatches (no skip): {diffs} / {n}")
