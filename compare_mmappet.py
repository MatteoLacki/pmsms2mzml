#!/usr/bin/env python3
"""Read mmappet precursor data and compare order with MGF."""
import numpy as np
import struct, re

base = 'reduced_input/pmsms.mmappet'
prec_dir = f'{base}/filtered_precursors_with_nontrivial_ms2.mmappet'
didx_dir = f'{base}/dataindex.mmappet'

# Read precursor schema
with open(f'{prec_dir}/schema.txt') as f:
    prec_schema = f.read().strip()
print(f"Precursor schema: {prec_schema}")

# Read precursor columns based on schema
# Schema: precursor_id_before_deisotoping(i64) frame(i32) scan(i32) tof(i32)
#         inv_ion_mobility(f64) intensity(u32) mz(f64) rt(f64) charge(u8)
prec_mz = np.fromfile(f'{prec_dir}/6.bin', dtype=np.float64)
prec_charge = np.fromfile(f'{prec_dir}/8.bin', dtype=np.uint8)
prec_iim = np.fromfile(f'{prec_dir}/4.bin', dtype=np.float64)
prec_rt = np.fromfile(f'{prec_dir}/7.bin', dtype=np.float64)

# Read dataindex
didx_size = np.fromfile(f'{didx_dir}/1.bin', dtype=np.uint64)

n_prec = len(prec_mz)
n_didx = len(didx_size)
print(f"Precursors: {n_prec}, Dataindex entries: {n_didx}")

# Which precursors have 0 fragments?
empty = [i for i in range(min(n_prec, n_didx)) if didx_size[i] == 0]
print(f"Empty precursors (0 fragments): {len(empty)}")
for i in empty:
    print(f"  [{i}] mz={prec_mz[i]:.6f} ch={prec_charge[i]} iim={prec_iim[i]:.4f}")

# Parse MGF for comparison
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
print(f"\nMGF spectra: {len(mgf)}")

# Compare mmappet order vs MGF order
print(f"\nmmappet vs MGF around index 254:")
print(f"{'idx':>4} {'mmappet_mz':>14} {'mm_ch':>5} {'mm_size':>8} | {'mgf_mz':>14} {'mgf_ch':>5} | match?")
print("-" * 80)
mm_i = 0
mgf_i = 0
for display_i in range(260):
    # Skip empty mmappet entries
    while mm_i < n_prec and mm_i < n_didx and didx_size[mm_i] == 0:
        mm_i += 1
    if mm_i >= n_prec or mgf_i >= len(mgf):
        break

    if display_i >= 251:
        mm_mz3 = f"{prec_mz[mm_i]:.3f}"
        mgf_mz3 = f"{mgf[mgf_i][0]:.3f}"
        match = "OK" if mm_mz3 == mgf_mz3 and prec_charge[mm_i] == mgf[mgf_i][1] else "DIFF"
        print(f"{display_i:4d} {prec_mz[mm_i]:14.6f} {prec_charge[mm_i]:5d} {didx_size[mm_i]:8d} | {mgf[mgf_i][0]:14.3f} {mgf[mgf_i][1]:5d} | {match}")

    mm_i += 1
    mgf_i += 1
