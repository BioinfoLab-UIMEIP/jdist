 jdist

**jdist** computes exact Jaccard distances between genomes or samples using
binary presence/absence matrices of genomic features, with GPU acceleration
(OpenCL) for scalable all-versus-all comparisons.

## Features
- Exact Jaccard distance (no sketching)
- GPU-accelerated via OpenCL; CPU fallback if GPU is unavailable
- Works with any binary feature matrix (unitigs recommended)

## Input
Tab-delimited matrix: rows = features, columns = genomes/samples, values = 0/1.

## Output
Tab-delimited symmetric distance matrix with sample IDs as row/column headers.

## Compilation
g++ -O3 -fopenmp src/jdist.cpp -lOpenCL -o jdist

## Usage
./jdist input_matrix.tsv output_distances.tsv [num_threads]

## Citation
Torres RC et al. jdist: GPU-accelerated exact Jaccard distances for scalable genome-wide comparison. (in preparation)

## License
MIT License.
