#!/usr/bin/env python3
"""Comprehensive diff of two mzML files for FragPipe-relevant fields."""
import re, sys

def parse_spectra(path):
    with open(path) as f:
        text = f.read()
    spectra = []
    for m in re.finditer(r'<spectrum index="(\d+)".*?</spectrum>', text, re.DOTALL):
        block = m.group(0)
        def cv(name):
            p = re.search(rf'name="{re.escape(name)}" value="([^"]*)"', block)
            return p.group(1) if p else None
        spectra.append({
            'index': int(m.group(1)),
            'arrayLen': re.search(r'defaultArrayLength="(\d+)"', block).group(1),
            'prec_mz': cv('selected ion m/z'),
            'charge': cv('charge state'),
            'rt': cv('scan start time'),
            'iim': cv('inverse reduced ion mobility'),
            'tic': cv('total ion current'),
            'base_mz': cv('base peak m/z'),
            'base_int': cv('base peak intensity'),
            'low_mz': cv('lowest observed m/z'),
            'high_mz': cv('highest observed m/z'),
            'title': cv('spectrum title'),
        })
    return spectra

a_spectra = parse_spectra(sys.argv[1])
b_spectra = parse_spectra(sys.argv[2])
n = min(len(a_spectra), len(b_spectra))

fields = ['prec_mz','charge','rt','iim','tic','base_mz','base_int','low_mz','high_mz','arrayLen','title']
for field in fields:
    diffs = 0
    examples = []
    for i in range(n):
        va, vb = a_spectra[i][field], b_spectra[i][field]
        if va != vb:
            diffs += 1
            if len(examples) < 3:
                examples.append((i, va, vb))
    status = f"{diffs}/{n} differ" if diffs else "ALL MATCH"
    print(f"{field:>12}: {status}")
    for idx, va, vb in examples:
        print(f"              [{idx}] A={va}  B={vb}")
