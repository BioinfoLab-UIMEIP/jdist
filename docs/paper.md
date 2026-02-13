# Paper (selected sections)

This page includes short excerpts from the manuscript describing the motivation and method, intended for readers of the software repository. The full manuscript is in preparation.

## Motivation
Genome distance measures are widely used for clustering, population structure analysis, outbreak detection, and downstream phylogenetic or epidemiological workflows. Core-SNP approaches provide high resolution but depend on alignment/mapping and—by definition—exclude accessory variation. Alignment-free sketching tools (e.g., Mash) are extremely scalable but rely on subsampling (MinHash), which can reduce sensitivity in fine-scale comparisons.

In highly modular systems (e.g., some plasmid families), stable homologous backbones can be limited, making “core” definitions problematic. Exact set-based similarity over complete presence/absence profiles remains well-defined and captures shared structural features even under weak sequence-level homology.

## Method overview
`jdist` computes **exact Jaccard distances** over binary feature profiles. Feature vectors are packed into 64-bit words and compared using bitwise AND/OR operations plus population counts. The all-vs-all computation is accelerated using a hybrid CPU+GPU OpenCL strategy.

If you are interested in the manuscript or benchmarking details, please contact the corresponding author.
