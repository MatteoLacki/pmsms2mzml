#!/usr/bin/env python3
"""Compare parquet fragment_spectrum_start/fragment_event_cnt vs dataindex idx/size."""
import numpy as np
import pyarrow.parquet as pq

base = 'reduced_input/pmsms.mmappet'
didx_dir = f'{base}/dataindex.mmappet'

# Dataindex
didx_idx = np.fromfile(f'{didx_dir}/2.bin', dtype=np.uint64)
didx_size = np.fromfile(f'{didx_dir}/1.bin', dtype=np.uint64)

# Parquet
pf = pq.read_table(f'{base}/filtered_precursors_with_nontrivial_ms2.parquet')
pq_start = pf.column('fragment_spectrum_start').to_numpy()
pq_cnt = pf.column('fragment_event_cnt').to_numpy()

print(f"Dataindex: {len(didx_idx)} entries")
print(f"Parquet: {len(pq_start)} entries")

# Check the 3 empty dataindex entries
empties = [i for i in range(min(len(didx_size), len(pq_cnt))) if didx_size[i] == 0]
print(f"\nEmpty dataindex entries ({len(empties)}):")
for i in empties:
    print(f"  [{i}] didx: idx={didx_idx[i]} size={didx_size[i]} | parquet: start={pq_start[i]} cnt={pq_cnt[i]}")

# General comparison around index 254
print(f"\nAround index 254:")
print(f"{'idx':>4} {'didx_idx':>10} {'didx_size':>10} {'pq_start':>10} {'pq_cnt':>10} {'match':>6}")
for i in range(251, 260):
    match = "OK" if didx_idx[i] == pq_start[i] and didx_size[i] == pq_cnt[i] else "DIFF"
    print(f"{i:4d} {didx_idx[i]:10d} {didx_size[i]:10d} {pq_start[i]:10d} {pq_cnt[i]:10d} {match:>6}")

# How many total mismatches?
n = min(len(didx_idx), len(pq_start))
idx_diff = np.sum(didx_idx[:n] != pq_start[:n])
size_diff = np.sum(didx_size[:n] != pq_cnt[:n])
print(f"\nTotal idx mismatches: {idx_diff} / {n}")
print(f"Total size mismatches: {size_diff} / {n}")
