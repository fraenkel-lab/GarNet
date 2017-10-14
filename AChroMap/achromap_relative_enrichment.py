import sys, os

import numpy as np
import pandas as pd

from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests



open_matches = "/nfs/latdata/iamjli/ALS/network_analysis/iMNs_ALS_CTR_20171004/relative_enrichment_merged_genes/motif_counts.10kb_diffSites.txt"
all_matches = "/nfs/latdata/iamjli/ALS/network_analysis/iMNs_ALS_CTR_20171004/relative_enrichment_merged_genes/motif_counts.10kb_window.txt"


def hypergeom_cdf(row):
	N = row[0]+row[1]
	N_pos = row[2]+row[3]
	sampled = row[2]
	sampled_pos = row[3]

	rv = hypergeom(N, N_pos, sampled)
	p_val = sum(rv.pmf(np.arange(sampled_pos, sampled+1)))

	return p_val


def main():

	motif_matches = {}

	# Counts for windows around all genes
	for line in open(all_matches, 'r').readlines(): 

		count, gene_state, motif = line.strip().split()

		# initialize motif
		if motif not in motif_matches: motif_matches[motif] = {}

		if int(gene_state) == 0: # motifs near non-DEGs
			motif_matches[motif]["nonDEGs"] = int(count)

		if int(gene_state) == 1: # motifs near DEGs
			motif_matches[motif]["DEGs"] = int(count)


	# Counts for only differential open chromatin around all genes
	for line in open(open_matches, 'r').readlines():

		count, gene_state, motif = line.strip().split()

		# initialize motif
		if motif not in motif_matches: motif_matches[motif] = {}

		if int(gene_state) == 0: # motifs near non-DEGs
			motif_matches[motif]["open_nonDEGs"] = int(count)

		if int(gene_state) == 1: # motifs near DEGs
			motif_matches[motif]["open_DEGs"] = int(count)


	df = pd.DataFrame.from_dict(motif_matches, orient="index").fillna(0).astype(int)
	df["p_val"] = df.apply(hypergeom_cdf, axis=1)
	df["FDR"] = multipletests(df["p_val"])[1]
	df.sort_values(by="FDR", ascending=True, inplace=True)
	
	df.to_csv("/nfs/latdata/iamjli/ALS/network_analysis/iMNs_ALS_CTR_20171004/relative_enrichment_merged_genes/results.10kb_window.tsv", sep='\t')



main()