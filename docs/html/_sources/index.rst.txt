Welcome to GarNet's technical documentation!
============================================

.. _Readme: https://github.com/fraenkel-lab/GarNet
.. _Issues: https://github.com/fraenkel-lab/GarNet/issues

GarNet is an epigenomics data exploration tool. For more information about the scientific uses of GarNet,
please see the Readme_. To report an error, please refer to the Issues_.

Some important notes:

- GarNet uses tab-delimited files, otherwise known as tsv's. BED files are already tsv, so no worries there.

- GarNet saves intermediate files for all input datasets, which it can use instead of raw tsv's saving computation and time. Those are pickled IntervalTrees (files which end in .pickle)

- We've supplied such intermediate files for motifs and genes in the repo, leaving it to you to supply the peaks file and the expression file.


.. automodule:: garnet
	:members:




