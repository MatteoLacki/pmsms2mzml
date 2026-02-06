# pmsms2mzml

Tool to translate pmsms format into mzml compatible with Fragpipe.

# Input:
* pmsms.mmappet
 matteo@pingu$ tree -h pmsms.mmappet/
[ 228]  pmsms.mmappet/
├── [502M]  0.bin # uint32 tof          this is what is inside 
├── [502M]  1.bin # uint32 intensity    this is what is inside     
├── [502M]  2.bin # float32 score       this is what is inside 
├── [502M]  3.bin # float32 mz          this is what is inside 
├── [ 110]  dataindex.mmappet
│   ├── [7.9M]  0.bin
│   ├── [7.9M]  1.bin
│   ├── [7.9M]  2.bin
│   ├── [3.9M]  3.bin
│   ├── [3.9M]  4.bin
│   └── [  80]  schema.txt
├── [ 58M]  filtered_precursors_with_nontrivial_ms2.parquet
└── [  52]  schema.txt # this file contains the schem for 0-3.bin.

This format contains raw fragment data:
* main folder overview: folder temp/F9477/cosine/pmsms.mmappet
```
              tof  intensity     score          mz
0           67109       2292  0.599361  234.146912
1           80687        944  0.399161  268.128693
...           ...        ...       ...         ...
131633098  255331        312  0.379669  910.427490
131633099  255539          0  0.464708  911.419312
[131633100 rows x 4 columns]
```
* index overview: folder temp/F9477/cosine/pmsms.mmappet/dataindex.mmappet
```
          ms1idx  size        idx  max_group_len  avg_group_len
0           5120   264          0             13       3.594697
1          13312   259        264              8       2.571429
...          ...   ...        ...            ...            ...
1032149  1030142   381  131632338              9       3.062992
1032150  1030143   381  131632719              9       3.062992

[1032151 rows x 5 columns]
```

* E.g. to get precursor's ms1idx tof spectrum look at 
```
idx = 0
max_group_len = 13
tof[idx:idx+max_group_len]
```
* corresponding output (without : example_output/mzml.mzML



# Example python writer

* src: example/mzMLwriter32.py 
* description:
    * python script that translates mgf into mzml.
* MAIN TASK:
    * replace mgf with pmsms directly and in another language.

## For Claude with love:
* update this file after each planning stage.
* update this file after each feature in the planning stage is implemented.
* use Makefile
    * ask when adding new entries into it.
* store outputs in output folder

## Deferred: Numpress compression
* example/mzML_numpress.py exists as reference but is NOT in scope right now.
* We will return to Numpress later. For now, use 32-bit float + zlib + base64 only.

## Input files (mmappet format, no parquet dependency)

### Fragment data: `pmsms.mmappet/`
Columnar binary: `0.bin`..`3.bin` + `schema.txt`
Schema: `uint32 tof, uint32 intensity, float32 score, float32 mz`

### Fragment index: `pmsms.mmappet/dataindex.mmappet/`
Schema: `uint64 ms1idx, uint64 size, uint64 idx, uint32 max_group_len, float32 avg_group_len`
For precursor i: fragments at `main[idx[i] .. idx[i]+size[i]]`

### Precursor metadata: `filtered_precursors_with_nontrivial_ms2.mmappet/`
Pre-converted from parquet to flat mmappet format (no parquet dependency needed).
Schema (9 columns):
| File | Type | Column |
|------|------|--------|
| 0.bin | int64 | precursor_id_before_deisotoping |
| 1.bin | int32 | frame |
| 2.bin | int32 | scan |
| 3.bin | int32 | tof |
| 4.bin | float64 | inv_ion_mobility |
| 5.bin | uint32 | intensity |
| 6.bin | float64 | mz |
| 7.bin | float64 | rt |
| 8.bin | uint8 | charge |

1,021,494 rows (filtered subset of dataindex's 1,032,151). This is the authoritative source for which precursors to emit. Row i here corresponds to dataindex row i.

## Precursor metadata → mzML mapping
| mmappet column              | mzML element                                  |
|-----------------------------|-----------------------------------------------|
| `mz`                        | `MS:1000744` selected ion m/z                 |
| `rt`                        | `MS:1000016` scan start time (seconds)        |
| `charge`                    | `MS:1000041` charge state                     |
| `inv_ion_mobility`          | `MS:1002815` inverse reduced ion mobility     |
| `precursor_id_before_deisotoping` | `MS:1000796` spectrum title              |
| dataindex `idx`             | fragment offset into main arrays              |
| dataindex `size`            | fragment count (== `fragment_event_cnt`)       |

## src/mmappet.hpp — header-only C++ mmappet reader/writer
* `MMappedData<T>` — mmap a single `.bin` column file, typed array access via `operator[]`, `.size()`, `.data()`
* `Dataset<T, Args...>` — variadic template opening all columns of an mmappet dir; validates schema types/names
* `Schema<T, Args...>` — defines a typed schema; `.open_dataset(path)` returns a `Dataset`
* `OpenDataset<T, Args...>(path, {"col1","col2",...})` — free function to open with schema+name validation
* Column access: `dataset.get_column<N>()` returns `MMappedData<T>&`
* Iterator: `for (auto [a,b,c] : dataset)` yields tuple per row
* Also: `DatasetWriter`, `IndexedDataset`, `IndexedWriter` (not needed for reading)

## Test data
* `input/testcutF9477.mmappet/` — 44 precursors, 9606 fragments. Cut from the full dataset.
* `example_output/tiny.mgf` + `example_output/tiny.mzml` — reference MGF/mzML for these 44 precursors.
* testcutF9477 has NO parquet file; metadata for these 44 comes from first 44 rows of the full precursors mmappet.

## Current implementation: main.cpp

Single-file C++ converter. Usage:
```
./myprog <pmsms_dir> <precursors_dir> <output.mzml> [--run-id NAME] [--zlib-level N]
```

Example:
```
./myprog input/testcutF9477.mmappet input/filtered_precursors_with_nontrivial_ms2_44.mmappet output/tiny_test.mzml --run-id tiny
```

Features implemented:
- Opens fragment data, dataindex, and precursor metadata via mmappet.hpp
- Writes indexed mzML with spectrum index and chromatogram index
- 32-bit float + zlib compression + base64 encoding for binary arrays
- Per-spectrum stats: lowest/highest mz, TIC, base peak mz/intensity
- Spectrum title format: `{run_id}.{idx}.{idx}.2 File:"", NativeID:"index={idx}"`
- CLI args: --run-id, --zlib-level

Verified: output matches reference (example_output/tiny.mzml) structurally — 44 spectra, 2074 lines, identical binary data. Minor value differences in metadata are expected (binary f32 precision vs MGF text rounding).

## General code quality patterns:

* avoid repetition but not at cost of too big abstraction.
* keep simple
* I will put tasks into AItodos folder. When we are done one task please move that text to AIdone folder and in Claude.nd do a summary of both the task and its final status.
