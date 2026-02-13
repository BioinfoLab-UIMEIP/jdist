# jdist

🌐 Documentation / GitHub Pages: https://bioinfolab-uimeip.github.io/jdist/

**jdist** computes **exact Jaccard distances** between genomes (or any samples) from binary feature presence/absence matrices using hybrid **CPU + GPU (OpenCL)** acceleration.

The tool is intended for analyses requiring genome-wide resolution without relying on heuristic sketching methods, while remaining computationally scalable.

## Features

- Exact Jaccard distances (no sketching)
- Genome-wide feature support (unitigs, k-mers, genes, etc.)
- GPU acceleration via OpenCL
- CPU fallback when GPU is unavailable
- Scales to large genome collections

## Installation

### Requirements
- C++17 compiler
- OpenMP support
- OpenCL headers and runtime

### Build
```bash
git clone https://github.com/BioinfoLab-UIMEIP/jdist
cd jdist
make
```

Binary produced:
```
./jdist
```

## Usage

```bash
./jdist <input.tsv> <output.tsv> [threads]
```

Example:
```bash
./jdist example.tsv distances.tsv 8
```

Input must be a tab-delimited binary presence/absence matrix:

- Rows: genomic features
- Columns: genomes/samples
- Values: 0/1
- First column: feature IDs
- First row: sample IDs

Output is a symmetric distance matrix ready for clustering or phylogenetic workflows.

## Citation

If you use **jdist**, please cite:

Torres RC, Meléndez-Sánchez D, Almaguer D, Torres J.  
*jdist: Exact Jaccard genome distances using GPU acceleration.*  
Manuscript in preparation.

## License

Released under the MIT License.

## Contact

Roberto C. Torres  
Medical Research Unit on Infectious and Parasitic Diseases (UIMEIP)  
Instituto Mexicano del Seguro Social (IMSS)
Mexico City, Mexico

