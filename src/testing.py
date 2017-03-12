#!/usr/bin/env python3

from garnet import *


def output(dataframe, output_dir):
	dataframe.to_csv(output_dir + 'output.tsv', sep='\t', header=True, index=False)

class Options:
	def __init__(self, options):
		self.__dict__.update(options)


# reference = "../data/ucsc_hg19_knownGenes.tsv" # was used to create
reference = "../data/reference_IntervalTree_dictionary.pickle"
kgXref = "../data/ucsc_hg19_kgXref.tsv"

# motifs = "../data/multiz46way_placental.txt" # was used to create
motifs = "../data/motifs_IntervalTree_dictionary.pickle"
# motifs = "../data/HUMAN_hg19_BBLS_1_00_FDR_0_10.bed"
# motifs = "../data/FOXA1_motifs.tsv"
# motifs = "../example/foxa1.motifs"

peaks = "../example/wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.bed"
# peaks = "../example/A549_FOXA1_broadPeak.bed"
# peaks = "../example/JUNB_ChIP_A549.bed"
# peaks = "../example/CTCF_ChIP_A549.bed"

expression = "../example/Tgfb_exp.tsv"

options = Options({"upstream_window": 2000,"downstream_window": 2000,"tss": False,"kgXref_file": kgXref,"output_dir":'/Users/alex/Documents/GarNet2/example/'})


output(map_motifs_to_peaks(motifs, peaks, options), options.output_dir)

# output(TF_regression(df, expression, options), options.output_dir)
