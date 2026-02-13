# jdist

## Exact genome distance computation with GPU acceleration

![License](https://img.shields.io/badge/license-MIT-blue)
![Language](https://img.shields.io/badge/language-C++17-blue)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)
![Status](https://img.shields.io/badge/status-active-success)


**jdist** computes **exact Jaccard distances** between genomes (or any samples) from binary feature presence/absence matrices using hybrid **CPU + GPU (OpenCL)** acceleration.

The tool is designed for analyses requiring **genome-wide resolution** without relying on heuristic sketching methods, while remaining computationally scalable.

✔ Exact distances (no sketching)  
✔ Genome-wide feature support  
✔ GPU acceleration via OpenCL  
✔ Scales to thousands of genomes  

---

## Quick start

```bash
git clone https://github.com/BioinfoLab-UIMEIP/jdist
cd jdist
make
./jdist examples/example_matrix.tsv out.tsv 4

- Installation: [install](install.md)
- Usage: [usage](usage.md)
```
---

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

Example:

```
Feature   G1   G2   G3
f1        1    0    1
f2        0    1    1
...
```

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

```
      A      B      C
A     0    0.12   0.31
B    0.12    0    0.28
C    0.31  0.28     0
```

Distances correspond to **Jaccard distance**:

\[
d = 1 - \frac{|A \cap B|}{|A \cup B|}
\]

The matrix can be used directly in clustering, PCoA, phylogenetic or population analyses.

## Notes
- Empty samples (all zeros) are automatically filtered.
- If a GPU is unavailable, the program falls back to a CPU OpenCL device.

---

## When should I use jdist?

`jdist` is particularly useful when:

- accessory genome variation is large,
- core genome alignments are unstable,
- plasmids or modular genomes are analyzed,
- sketching methods lack resolution,
- exact distances are required,
- large genome collections must be compared efficiently.

---

## Method overview

`jdist` computes **exact Jaccard similarity** over binary genomic feature profiles.

Feature vectors are packed into 64-bit words and compared using bitwise AND/OR operations plus population counts. Pairwise comparisons are parallelized using a hybrid CPU multithreading + GPU OpenCL strategy, enabling efficient all-versus-all computations at population scale.

---

## Supplementary figures

Supplementary benchmarking figures are available:

- Fig. S1 — Rank-ordered distance overlay  
  [PDF](assets/pdf/10_rank_refdist_overlay.pdf)

- Fig. S2 — All-pairs scatter comparison  
  [PDF](assets/pdf/01_allpairs_scatter_vs_snp_raw.pdf)

- Fig. S3 — Distance distributions  
  [PDF](assets/pdf/05_distance_distributions.pdf)

- Fig. S4 — Procrustes: jdist vs Mash  
  [PDF](assets/pdf/08_procrustes_jaccard_vs_mash.pdf)

- Fig. S5 — Procrustes displacement histograms  
  [PDF](assets/pdf/09_procrustes_displacement_hists.pdf)

---

## Limitations

- Memory usage scales with number of features.
	Filtering extremely rare variants is recommended, as they often contribute little to global similarity estimates while substantially increasing computational cost.
- Generation of presence/absence matrices is external to the software.

---

## Citation

If you use **jdist**, please cite:

Torres RC, Meléndez-Sánchez D, Almaguer D, Torres J.  
*jdist: Exact Jaccard genome distances using GPU acceleration.*  
Manuscript in preparation.

---

## Contact

Roberto C. Torres, PhD.
Medical Research Unit on Infectious and Parasitic Diseases (UIMEIP)  
Instituto Mexicano del Seguro Social (IMSS)
Mexico City, Mexico
