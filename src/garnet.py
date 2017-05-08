#!/usr/bin/env python3

# Core python modules
import sys
import os

# Peripheral python modules
import pickle
import logging
import sqlite3

# Core python external libraries
import numpy as np
import pandas as pd
from statsmodels.formula.api import ols as linear_regression
from statsmodels.graphics.regressionplots import abline_plot as plot_regression

# Peripheral python external libraries
import jinja2

templateLoader = jinja2.FileSystemLoader(searchpath=".")
templateEnv = jinja2.Environment(loader=templateLoader)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - GarNet: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)


######################################## File Parsing Logic #######################################

def parse_peaks_file(peaks_file):
	"""
	Parse a BED file with peaks from an epigenomics assay (e.g. ATAC) into a dataframe

	Arguments:
		peaks_file (string or FILE): BED file from epigenomics assay

	Returns:
		dataframe: peaks dataframe
	"""

	# if peaks file format is BED

	peaks_fieldnames = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]

	peaks_dataframe = pd.read_csv(peaks_file, delimiter='\t', names=peaks_fieldnames)

	peaks_dataframe.rename(index=str, columns={"chromStart":"peakStart", "chromEnd":"peakEnd", "name":"peakName", "score":"peakScore", "strand":"peakStrand"}, inplace=True)

	peaks_dataframe = peaks_dataframe[['peakName', 'chrom', 'peakStart', 'peakEnd', 'peakScore']]

	# if peaks file format is MACS

	# peaks_fieldnames = ["chr", "start", "end", "length", "summit", "tags", "-10*log10(pvalue)", "fold_enrichment FDR(%)"]

	# if peaks file format is GEM

	# peaks_fieldnames = ["Location", "IP binding strength", "Control binding strength", "Fold", "Expected binding strength", "Q_-lg10", "P_-lg10", "P_poiss", "IPvsEMP", "Noise", "KmerGroup", "KG_hgp", "Strand"]

	return peaks_dataframe


def parse_expression_file(expression_file):
	"""
	Parse gene expression scores from a transcriptomics assay (e.g. RNAseq) into a dataframe

	Arguments:
		expression_file (string or FILE): Two-column, tab-delimited file of gene / gene expression score

	Returns:
		dataframe: expression dataframe
	"""

	return pd.read_csv(expression_file, delimiter='\t', names=["name", "expression"])



def output(dataframe, output_dir, filename):
	dataframe.to_csv(os.path.join(output_dir, filename), sep='\t', header=True, index=False)



######################################### Public Functions #########################################

def map_peaks(garnet_file, peaks_file_or_list_of_peaks_files, options):
	"""
	Find motifs and associated genes local to peaks.

	This function searches for motifs "under" peaks from an epigenomics dataset and "around" peaks for genes.
	It then returns all pairs of motifs and genes which were found local to peaks.

	Arguments:
		garnet_file (str): filepath or file object for the garnet file.
		peaks_file_or_list_of_peaks_files (str or FILE or list): filepath or file object for the peaks file, or list of such paths or objects
		options (dict): {"upstream_window": int, "downstream_window": int, "tss": bool, "output_dir": string (optional)})

	Returns:
		dataframe: a dataframe with rows of transcription factor binding motifs and nearby genes
			with the restriction that these motifs and genes must have been found near a peak.
	"""

	GarNetDB = sqlite3.connect(garnet_file)

	# peaks_file_or_list_of_peaks_files is either a filepath or FILE, or a list of filepaths or FILEs.
	# Let's operate on a list in either case, so if it's a single string, put it in a list. #TODO, this will break if it's a single FILE.
	if isinstance(peaks_file_or_list_of_peaks_files, basestring): peaks_files = [peaks_file_or_list_of_peaks_files]
	else: peaks_files = peaks_file_or_list_of_peaks_files

	output = []

	for peaks_file in peaks_files:

		peaks = parse_peaks_file(peaks_file)

		peaks_with_associated_genes_and_motifs = query_db(peaks, GarNetDB)

		# Should probably map type_of_peak here

		output.append(motifs_and_genes)

	GarNetDB.close()

	# conversely, if this function was passed a single file, return a single dataframe
	if len(output) == 1: output = output[0]
	return output


def TF_regression(motifs_and_genes_dataframe, expression_file, options):
	"""
	Perform linear regression of the expression of genes versus the strength of the associated
	transcription factor binding motifs and report results.

	This function parses an expression file of two columns: gene symbol and expression value, and
	merges the expression profile into the motifs and genes file, resulting in information about
	transcription factor binding motifs local to genes, and those genes' expressions. We do linear
	regression, and if an output directory is provided, we output a plot for each TF and an html
	summary of the regressions.

	Arguments:
		motifs_and_genes_dataframe (dataframe): the outcome of map_known_genes_and_motifs_to_peaks
		expression_file (str or FILE): a tsv file of expression data, with geneSymbol, score columns
		options (dict): {"output_dir": string (optional)})

	Returns:
		dataframe: slope and pval of linear regfression for each transcription factor.
	"""

	expression_dataframe = parse_expression_file(expression_file)

	motifs_genes_and_expression_levels = motifs_and_genes_dataframe.merge(expression_dataframe, left_on='geneSymbol', right_on='name', how='inner')

	# the same geneSymbol might have different names but since the expression is geneSymbol-wise
	# these additional names cause bogus regression p-values. Get rid of them here.
	if 'geneSymbol' in motifs_genes_and_expression_levels.columns:
		motifs_genes_and_expression_levels.drop_duplicates(subset=['geneSymbol', 'motifID'], inplace=True)
	motifs_genes_and_expression_levels['motifScore'] = motifs_genes_and_expression_levels['motifScore'].astype(float)

	TFs_and_associated_expression_profiles = list(motifs_genes_and_expression_levels.groupby('motifName'))
	imputed_TF_features = []
	logger.info("Performing linear regression on "+str(len(TFs_and_associated_expression_profiles))+" transcription factor expression profiles...")

	for TF_name, expression_profile in TFs_and_associated_expression_profiles:

		# Occasionally there's only one gene associated with a TF, which we can't fit a line to.
		if len(expression_profile) < 5: continue

		# Ordinary Least Squares linear regression
		result = linear_regression(formula="expression ~ motifScore", data=expression_profile).fit()

		if options.get('output_dir'):
			plot = plot_regression(model_results=result, ax=expression_profile.plot(x="motifScore", y="expression", kind="scatter", grid=True))
			if not os.path.exists(options['output_dir']+'regression_plots/'): os.makedirs(options['output_dir']+'regression_plots/')
			plot.savefig(options['output_dir']+'regression_plots/' + TF_name + '.png')

		imputed_TF_features.append((TF_name, result.params['motifScore'], result.pvalues['motifScore']))

	imputed_TF_features_dataframe = pd.DataFrame(imputed_TF_features, columns=["Transcription Factor", "Slope", "P-Value"])

	# If we're supplied with an output_dir, we'll put a summary html file in there as well.
	if options.get('output_dir'):
		html_output = templateEnv.get_template("summary.jinja").render(images_dir=options['output_dir']+'regression_plots/', TFs=sorted(imputed_TF_features, key=lambda x: x[2]))
		with open(options['output_dir']+"summary.html", "w") as summary_output_file:
			summary_output_file.write(html_output)

	return imputed_TF_features_dataframes


######################################## Private Functions ########################################

def type_of_peak(row):
	"""
	Arguments:
		row (pd.Series): A row of data from a dataframe with peak and gene information

	Returns:
		str: a name for the relationship between the peak and the gene:
				- upstream if the start of the peak is more than 2kb above the start of the gene
				- promoter if the start of the peak is above the start of the gene
				- downstream if the start of the peak is below the start of the gene
	"""

	if row['geneStrand'] == '+':
		if -2000 >= row['peakStart'] - row['geneStart']: 	return 'upstream'
		if -2000 < row['peakStart'] - row['geneStart'] < 0: return 'promoter'
		if 0 <= row['peakStart'] - row['geneStart']: 		return 'downstream'  # a.k.a. row['peakStart'] < row['geneStart']
		return 'intergenic'
	if row['geneStrand'] == '-':
		if 2000 <= row['peakEnd'] - row['geneEnd']: 	return 'upstream'
		if 2000 > row['peakEnd'] - row['geneEnd'] > 0: 	return 'promoter'
		if 0 >= row['peakEnd'] - row['geneEnd']: 		return 'downstream'  # a.k.a. row['peakEnd'] < row['geneEnd']
		return 'intergenic'


def query_db(peaks, GarNetDB):
	"""

	Arguments:
		peaks (pd.DataFrame)
		GarNetDB (squlite3.connection):

	"""

	peaks.to_sql('peaks', GarNetDB, if_exists="replace")

	overlaps = pd.read_sql_query("""
		SELECT garnetdb.*
		FROM garnetdb
		JOIN peaks ON garnetdb.chr        == peaks.chr
				  AND garnetdb.motifStart <= peaks.peakEnd
				  AND garnetdb.motifEnd   >= peaks.peakStart;
		""", GarNetDB)

	return overlaps

