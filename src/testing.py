#!/usr/bin/env python3

from garnet import *


def output(dataframe, output_dir):
	dataframe.to_csv(output_dir + 'output.tsv', sep='\t', header=True, index=False)

class Options:
	def __init__(self, options):
		self.__dict__.update(options)

peaks = "../example/wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.bed"
# reference = "../data/ucsc_hg19_knownGenes.tsv"
# motifs = "../data/HUMAN_hg19_BBLS_1_00_FDR_0_10.bed"

kgXref = "../data/ucsc_hg19_kgXref.tsv"
expression = "../data/Tgfb_exp.tsv"

# peaks = "../data/peaks_IntervalTree_dictionary.pickle"
reference = "../data/reference_IntervalTree_dictionary.pickle"
motifs = "../data/motifs_IntervalTree_dictionary.pickle"

# peaks = "../example/A549_FOXA1_broadPeak.bed"
# motifs = "../example/foxa1.motifs"

options = Options({"upstream_window": 2000,"downstream_window": 2000,"tss": False,"kgXref_file": kgXref,"output_dir":'/Users/alex/Documents/GarNet2/example/'})


output(map_known_genes_to_peaks(reference, peaks, options), options.output_dir)

# output(TF_regression(df, expression, options), options.output_dir)

