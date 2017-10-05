
import sys, os

import numpy as np
import pandas as pd

from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


binned_peaks = "/nfs/latdata/iamjli/ALS/network_analysis/iMNs_ALS_CTR_20171004/TF_prediction/ALS_CTR.diffPeaks.2kb_window.100bp_bins.tsv"
binned_motifs = "/nfs/latdata/iamjli/ALS/network_analysis/iMNs_ALS_CTR_20171004/TF_prediction/ALS_CTR.binned_TFs.2kb_window.100bp_bins.tsv"



def main():

	peaks_df = pd.read_csv(binned_peaks, sep='\t', names=["chrom", "start", "end", "gene", "isPeak"])
	motifs_df = pd.read_csv(binned_motifs, sep='\t', names=["chrom", "start", "end", "gene", "motif"]).drop_duplicates()

	TFs = motifs_df["motif"].unique()

	results = []

	for TF in TFs:
		TF_sites = motifs_df[motifs_df["motif"] == TF]
		merged = peaks_df.merge(TF_sites, on=["chrom", "start", "end", "gene"], how="outer")

		# There are N bins that are nearby DEG TSS, N_good of them contain open chromatin signature
		# There are n_samples predicted motifs, n_positives of them are also in open regions
		N, N_good = merged.shape[0], merged[merged["isPeak"]>0].shape[0]
		n_samples = merged.count()["motif"]
		n_positives = merged[merged["motif"].notnull() & merged["isPeak"]>0].shape[0]

		# For efficiency, calculate probability of seeing fewer than n_positives
		rv = hypergeom(N, N_good, n_samples)
		p_val = 1 - sum(rv.pmf(np.arange(0,n_positives)))

		results.append([TF, p_val, N, N_good, n_samples, n_positives])


	results_df = pd.DataFrame(results, columns=["TF_name", "p-value", "total_bins", "open_bins", "motifs_bins", "motifs-open_bins"])
	results_df["FDR"] = multipletests(results_df["p-value"])[1]
	results_df.sort_values(by="FDR", ascending=True, inplace=True)
	results_df = results_df[["TF_name", "p-value", "FDR", "total_bins", "open_bins", "motifs_bins", "motifs-open_bins"]]
	results_df.to_csv("/nfs/latdata/iamjli/ALS/network_analysis/iMNs_ALS_CTR_20171004/TF_prediction/results_motif_enrichment_in_open_chromatin.tsv", sep='\t', index=False)




main()