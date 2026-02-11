# pmsms2mzml

Converts pmsms.mmappet binary format directly to mzML, as FragPipe wants it.

## Dependencies

- **C++23 compiler** (clang++ or g++) (but C++17 worked too)
- **zlib** development headers (`sudo apt install zlib1g-dev` on Debian/Ubuntu)

## Build

```
make
```

Produces the `pmsms2mzml` binary.

## Usage

```
pmsms2mzml <pmsms_dir> <output> [options]
```

### Arguments

| Argument | Description |
|---|---|
| `pmsms_dir` | Path to the pmsms.mmappet directory |
| `output` | Output path (directory in default mode, file with `--indexed`) |

### Options

| Option | Description |
|---|---|
| `--threads N` | Number of threads (default: 1) |
| `--zlib-level N` | zlib compression level 1-9 (default: 1) |
| `--decimals N` | Truncate m/z values to N decimal places |
| `--indexed` | Produce indexed mzML (single file). Default is non-indexed |
| `--numpress` | Use MS-Numpress compression for m/z and intensity arrays |
| `--precursors-dir DIR` | Override precursor metadata directory |
| `--run-id NAME` | Override run ID (default: derived from output path) |
| `--check-mz-sorted` | Verify m/z values are sorted before writing; exit on failure |
| `--dry-run` | Compute output size without writing |
| `--used-spectra-cnt N` | Only convert the first N spectra |

### Output modes

**Default (non-indexed):** Output is a directory containing:
- `mzml.mzML` — the mzML file
- `idmap.mmappet/` — mapping from sequential mzML index to original precursor index (uint64 column)

The idmap allows remapping search engine results back to the original precursor indices.

**`--indexed`:** Output is a single indexed mzML file with spectrum offset index.

## License

MIT. See [LICENSE](LICENSE).

This project includes [MS-Numpress](https://github.com/ms-numpress/ms-numpress) (`src/MSNumpress.cpp`, `src/MSNumpress.hpp`), which is licensed under the Apache License 2.0. See the file headers for details.

## Input format

Expects a pmsms.mmappet directory with:
- `schema.txt` and columnar binary files (`0.bin`..`3.bin`): fragment-level data (tof, intensity, score, mz)
- `filtered_precursors_with_nontrivial_ms2.mmappet/`: precursor metadata (mz, rt, charge, inv_ion_mobility, fragment_spectrum_start, fragment_event_cnt)
