## GarNet 2.0

This repository is a fast re-write of GarNet, originally written by [Sara Gosline](https://github.com/sgosline) as part of [OmicsIntegrator](https://github.com/fraenkel-lab/omicsintegrator).

The goal of GarNet is to use gene expression and epigenetic data to impute transcription factors (TFs) that played an important role in a biological system. Transcription factors bind in open chromatin regions to specific DNA sequences called "motifs," and affect the expression levels of nearby genes.
To determine which TFs were relevant to a biological system, users should supply epigenetic regions (peaks) of interest (i.e. open chromatin regions derived from ATAC-seq or DNase-seq on your cells or in a similar cell line) and differential gene expression data. GarNet:

1. Looks for known TF motifs (derived from [MotifMap](http://motifmap-rna.ics.uci.edu/) ) that occur within your epigenetic regions
2. Looks for known genes (derived from [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) ) that occur near your epigenetic regions
3. Maps the TFs and genes that were found near the same peaks to each other as those TFs potentially effect the expression of those genes
4. Uses linear regression to see if the change in expression level is dependent on the strength of the TF motif.

If motifs and genes are found near each other, inside relevant epigenetic regions, and expression level is significantly dependent on motif strength, we conclude that the TF is likely an important player in changing the gene expression in your system.

[Documentation about each of the functions can be found here.](https://fraenkel-lab.github.io/GarNet2/html/index.html)

TODO (Alex): Information about how to set up code.


---

#### Notes to ourselves:

Alternative IntervalTree implementations:

- https://gist.github.com/shoyer/c939325f509d7c027949 (keep an eye on https://github.com/pandas-dev/pandas/pull/8707 which has now moved to https://github.com/pandas-dev/pandas/pull/15309)
- https://github.com/ekg/intervaltree
- https://github.com/cpcloud/banyan


Things worth reading at some point:

- http://chrisalbon.com/python/pandas_apply_operations_to_dataframes.html
- http://fastinterval.readthedocs.io/en/latest/





