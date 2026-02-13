# jdist

GPU-accelerated exact Jaccard distances for genome-wide comparison.

- Installation: [install](install.md)
- Usage: [usage](usage.md)
- Paper excerpt: [paper](paper.md)
- Supplementary figures: [figures](figures.md)

# Installation

## Dependencies
- C++17 compiler
- OpenMP
- OpenCL headers + runtime (GPU optional; CPU OpenCL fallback supported)

## Build
```bash
make
```

This produces `./jdist`.

## Minimal test
```bash
make run
```

# Usage

## Input format
`jdist` expects a **tab-delimited binary presence/absence matrix**:

- **Rows**: features (unitigs, k-mers, gene families, or any user-defined features)
- **Columns**: samples/genomes
- **Values**: `0/1` (absence/presence)

The first column contains feature IDs; the first row contains sample IDs.

## Run
```bash
./jdist <input.tsv> <output.tsv> [num_threads]
```

Example:
```bash
./jdist examples/example_matrix.tsv examples/example_dist.tsv 4
```

## Output
A symmetric tab-delimited **distance matrix** with sample IDs as row/column headers. Values correspond to **Jaccard distance**:

\\[
d = 1 - \\frac{|A \\cap B|}{|A \\cup B|}
\\]

## Notes
- Empty samples (all zeros) are automatically filtered.
- If a GPU is unavailable, the program falls back to a CPU OpenCL device.

## Citation
If you use **jdist**, please cite:
- Torres RC, Meléndez-Sánchez D, Torres J, Almaguer D. *jdist: GPU-accelerated exact Jaccard distances for scalable genome-wide comparison.* (manuscript in preparation).
