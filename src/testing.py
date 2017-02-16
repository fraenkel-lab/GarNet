#!/usr/bin/env python3

from garnet import *


def output(dataframe, output_dir):
	dataframe.to_csv(output_dir + 'output', sep='\t', header=True, index=False)


peaks = "../example/wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.bed"
# reference = "../data/ucsc_hg19_knownGenes.tsv"
# motifs = "../data/HUMAN_hg19_BBLS_1_00_FDR_0_10.bed"

kgXref = "../data/ucsc_hg19_kgXref.tsv"
# expression = "../data/Tgfb_exp.tsv"
# output_file = "../src/output"

# peaks = "../data/peaks_IntervalTree_dictionary.pickle"
reference = "../data/reference_IntervalTree_dictionary.pickle"
motifs = "../data/motifs_IntervalTree_dictionary.pickle"

# peaks = "../example/A549_FOXA1_broadPeak.bed"
# motifs = "../example/foxa1.motifs"

output(map_known_genes_and_motifs_to_peaks(reference, motifs, peaks,
	{"upstream_window": 2000,
	 "downstream_window": 2000,
	 "tss": False,
	 "kgXref_file": kgXref,
	 "output_dir":'/Users/alex/Documents/GarNet2/src/'}), '/Users/alex/Documents/GarNet2/src/')

