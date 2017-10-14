#!/usr/bin/env python3

# Core python modules
import sys
import os

# Peripheral python modules
import argparse
import logging

# Peripheral python external libraries
from pybedtools import BedTool

# python external libraries
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - MotifEnrich: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)


parser = argparse.ArgumentParser(description="Find enrichment of motifs in open chromatin near differentially expressed genes.")

class FullPaths(argparse.Action):
	"""Expand user- and relative-paths"""
	def __call__(self,parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def directory(dirname):
	if not os.path.isdir(dirname): raise argparse.ArgumentTypeError(dirname + " is not a directory")
	else: return dirname


# Input / Output parameters:
io_params = parser.add_argument_group("Input / Output Files")

io_params.add_argument("-g", "--DEG", dest='DEG_file', type=str, required=True,
	help ='(Required) File path to differentially expressed genes separated by line. Will ignore any columns but the first.')
io_params.add_argument("-s", "--DOS", dest='DOS_file', type=str, required=True,
	help ='(Required) File path to differentially open peaks with 3 columns: "chrom\tstart\tend')
io_params.add_argument('-o', '--output', dest='output_dir', action=FullPaths, type=directory, required=True,
	help='(Required) Output directory path')

# Command parameters 
search_params = parser.add_argument_group("Search Parameters")

search_params.add_argument("-d", dest="d", type=int, required=False, default=2000,
	help="Window around each gene's TSS to search for motifs up to 10kb. [default: 2kb]")
search_params.add_argument("-b", dest="b", type=int, required=False, default=100,
	help="Bin size. [default: 100bp]")
search_params.add_argument("-genome", dest="genome", type=str, required=False, default="hg19",
	help="Genome identifier. [default: 'hg19']")



def filter_TSS_file_by_DEGs(TSS_file, DEG_file):

	TSS_df = pd.read_csv(TSS_file, sep='\t', names=["chrom", "start", "end", "gene", "score", "strand"])
	DEG_list = [line.split()[0] for line in open(DEG_file, 'r').readlines()] # Get first item in each line

	filtered_TSS_df = TSS_df[TSS_df["gene"].isin(DEG_list)]
	logger.info("{} DEGs specified, {} found in reference file.".format(len(set(DEG_list)), len(set(filtered_TSS_df.gene.tolist()))))

	return filtered_TSS_df


def hypergeometric_cdf(N, N_pos, n, n_pos):

	# P-value does not change for large N, so we scale N to avoid floating point errors. 
	if N > 100000:
		N_pos = N_pos/N*100000
		N = 100000

	rv = hypergeom(N, N_pos, n)
	p_val = sum(rv.pmf(np.arange(n_pos, n+1)))

	return min(p_val, 1)


def TF_enrichment(TF, motif_df, DOS_df): 

	TF_sites = motif_df[motif_df["motif"] == TF]
	df = DOS_df.merge(TF_sites, on=["chrom", "start", "end", "gene"], how="outer")

	N, N_pos = df.shape[0], df[df["isPeak"]>0].shape[0]
	n = df.count()["motif"]
	n_pos = df[df["motif"].notnull() & df["isPeak"]>0].shape[0]

	p_val = hypergeometric_cdf(N, N_pos, n, n_pos)

	return [TF, p_val, N, N_pos, n, n_pos]


def main():

	args = parser.parse_args()

	TSS_file = "../motif_matches/reference/ucsc_reference.{}.tss.bed".format(args.genome)
	genome_size_file = "../motif_matches/genomes/{}.chrom.sizes".format(args.genome)
	motifs_file = "../motif_matches/motifs.cisBP.{}.10kb.1e-05.sorted.bed".format(args.genome)

	TSS_bed = BedTool.from_dataframe(filter_TSS_file_by_DEGs(TSS_file , args.DEG_file))

	logger.info("Binning search windows.")
	windows = TSS_bed.slop(b=args.d, genome=args.genome)
	binned_windows = windows.window_maker(b=windows, w=args.b, i="src").sort()

	# First filter removes instances where a motif is assigned to 2 windows by overlapping the boundary of a window.
	# Second filter ensures that there is only one motif per bin, but may remove legitemate matches.
	logger.info("Sorting motifs into bins.")
	binned_motifs = binned_windows.intersect(b=motifs_file, wo=True, sorted=True)
	binned_uniq_motifs_df = binned_motifs.to_dataframe(names=["chrom", "start", "end", "gene", "chrom1", "start1", "end1", "motif", "score", "strand", "overlap"])\
	                                     .sort_values(by=["chrom1", "start1", "end1", "motif", "overlap"]) \
	                                     .drop_duplicates(subset=["gene", "chrom1", "start1", "end1", "motif"], keep="last")
	binned_uniq_motifs_df = binned_uniq_motifs_df[["chrom", "start", "end", "gene", "motif"]].drop_duplicates()

	# Count open chromatin open in bins
	logger.info("Sorting DOS into bins.")
	binned_DOS = binned_windows.intersect(b=args.DOS_file, c=True)
	binned_DOS_df = binned_DOS.to_dataframe(names=["chrom", "start", "end", "gene", "isPeak"])


	#### REGRESSION #####
	logger.info("Performing enrichment on {} TFs.".format(len(binned_uniq_motifs_df["motif"].unique())))
	results = [TF_enrichment(TF, binned_uniq_motifs_df, binned_DOS_df) for TF in binned_uniq_motifs_df["motif"].unique()]
	results_df = pd.DataFrame(results, columns=["TF_name", "p-value", "total_bins", "open_bins", "motifs_bins", "motifs-open_bins"])
	results_df["FDR"] = multipletests(results_df["p-value"], method='b')[1]
	results_df.sort_values(by="FDR", inplace=True)
	results_df = results_df[["TF_name", "p-value", "FDR", "total_bins", "open_bins", "motifs_bins", "motifs-open_bins"]]


	results_df.to_csv(args.output_dir+"/results.{}_window.{}_bins.tsv".format(args.d, args.b), sep='\t', index=False)


main()