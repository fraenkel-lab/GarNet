# GarNet

Many biological labs commonly use RNA-seq or other high-throughput technologies to assess the gene expression changes in their experimental systems. In order to produce these gene expression changes, regulatory machinary in the cell must have caused changes in transcription. Assessing the effect of each factor of this regulatory machinary directly is often difficult. An open question, then, is how to use gene expression data to predict which cis-regulatory elements, such as transcription factor proteins, were at play in producing those gene expression changes.

The goal of GarNet is to use gene expression and epigenetic data to impute transcription factors (TFs) that played an important role in a biological system. Transcription factors bind in open chromatin regions to specific DNA sequences called "motifs," and affect the expression levels of nearby genes.
To determine which TFs were relevant to a biological system, users should supply epigenetic regions (peaks) of interest (i.e. open chromatin regions derived from ATAC-seq or DNase-seq on your cells or in a similar cell line) and differential gene expression data. GarNet:

1. Looks for known TF motifs (derived from [MotifMap](http://motifmap-rna.ics.uci.edu/) ) that occur within your epigenetic regions
![map TFs to peaks](figures/Picture1.png)
2. Looks for known genes (derived from [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) ) that occur near your epigenetic regions
![map genes to peaks](figures/Picture2.png)
3. Maps the TFs and genes that were found near the same peaks to each other as those TFs potentially effect the expression of those genes
4. For each TF, uses linear regression to see if the change in expression level is dependent on the strength of the Transcription factor binding motif.
![Regress motif strength on expression](figures/Picture3.png)

If a Transcription Factor binding motif is found near genes changing in expression, inside relevant epigenetic regions in this tissue type, and changes in gene expression are significantly dependent on the strength of that motif, we predict that the TF is likely an important player in the gene expression in your system. We assign a score to that TF based on the significance and slope of the regression of motif strength on gene expression.

[Documentation about each of the functions can be found here.](https://fraenkel-lab.github.io/GarNet2/html/index.html)

TODO (Alex): Information about how to set up code.


This repository is an updated version of Garnet, originally written by [Sara Gosline](https://github.com/sgosline) and Anthony Soltis as part of [OmicsIntegrator](https://github.com/fraenkel-lab/omicsintegrator).


---

#### Notes to ourselves:

Alternative IntervalTree implementations:

- https://gist.github.com/shoyer/c939325f509d7c027949 (keep an eye on https://github.com/pandas-dev/pandas/pull/8707 which has now moved to https://github.com/pandas-dev/pandas/pull/15309)
- https://github.com/ekg/intervaltree
- https://github.com/cpcloud/banyan


Things worth reading at some point:

- http://chrisalbon.com/python/pandas_apply_operations_to_dataframes.html
- http://fastinterval.readthedocs.io/en/latest/





