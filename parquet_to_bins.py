#!/usr/bin/env python3
"""Add fragment_spectrum_start and fragment_event_cnt from parquet to precursor mmappet dir."""
import sys, os
import numpy as np
import pyarrow.parquet as pq

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <pmsms_dir>")
    sys.exit(1)

pmsms_dir = sys.argv[1]
parquet_path = os.path.join(pmsms_dir, "filtered_precursors_with_nontrivial_ms2.parquet")
mmappet_dir = os.path.join(pmsms_dir, "filtered_precursors_with_nontrivial_ms2.mmappet")

# Read parquet
table = pq.read_table(parquet_path)
frag_start = table.column('fragment_spectrum_start').to_numpy().astype(np.uint64)
frag_cnt = table.column('fragment_event_cnt').to_numpy().astype(np.uint64)

# Read existing schema to find next column index
schema_path = os.path.join(mmappet_dir, "schema.txt")
with open(schema_path) as f:
    lines = [l.strip() for l in f if l.strip()]
next_idx = len(lines)

# Write .bin files
frag_start.tofile(os.path.join(mmappet_dir, f"{next_idx}.bin"))
frag_cnt.tofile(os.path.join(mmappet_dir, f"{next_idx+1}.bin"))

# Append to schema
with open(schema_path) as f:
    content = f.read()
if not content.endswith('\n'):
    content += '\n'
with open(schema_path, 'w') as f:
    f.write(content)
    f.write(f"uint64 fragment_spectrum_start\n")
    f.write(f"uint64 fragment_event_cnt\n")

print(f"Added columns {next_idx} (fragment_spectrum_start) and {next_idx+1} (fragment_event_cnt)")
print(f"  {len(frag_start)} entries, total fragments: {frag_cnt.sum()}")
