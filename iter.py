#!/usr/bin/env python3

import subprocess

INPUT = "input/F9477_1000.mmappet"

for K in range(750, 1001):
    output = f"output/F9477_{K}_1thread.mzml"
    cmd = [
        "./myprog",
        INPUT,
        output,
        "--thread", "1",
        "--used-spectra-cnt", str(K),
    ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
