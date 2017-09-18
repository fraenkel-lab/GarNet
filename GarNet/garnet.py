#!/usr/bin/env python3

# Core python modules
import sys
import os

# Peripheral python modules
import pickle
import logging
import subprocess

# Core python external libraries
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from statsmodels.formula.api import ols as linear_regression
from statsmodels.graphics.regressionplots import abline_plot as plot_regression

# Peripheral python external libraries
from pybedtools import BedTool
import jinja2

# list of public methods:
__all__ = [ "construct_garnet_file", "map_peaks", "TF_regression" ]

templateLoader = jinja2.FileSystemLoader(os.path.dirname(os.path.abspath(__file__)))
templateEnv = jinja2.Environment(loader=templateLoader)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - GarNet: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)


######################################## File Parsing Logic #######################################

def parse_expression_file(expression_file):
	"""
	Parse gene expression scores from a transcriptomics assay (e.g. RNAseq) into a dataframe and ensures
	columns are of the correct type.

	Arguments:
		expression_file (string or FILE): Two-column, tab-delimited file of gene / gene expression score

	Returns:
		pd.dataframe: expression dataframe
	"""
	df = pd.read_csv(expression_file, delimiter='\t', names=["name", "expression"])
	df["expression"] = pd.to_numeric(df["expression"], errors="coerce")
	df.dropna(inplace=True)

	return df


def _parse_motifs_and_genes_file_or_dataframe(motifs_and_genes_file_or_dataframe):
	"""
	If the argument is a dataframe, return it. Otherwise if the argument is a string, try to read a dataframe from it, and return that
	"""

	if isinstance(motifs_and_genes_file_or_dataframe, str):
		motifs_and_genes_dataframe = pd.read_csv(motifs_and_genes_file_or_dataframe, delimiter='\t', header=0, index_col=False)

	elif isinstance(motifs_and_genes_file_or_dataframe, pd.DataFrame):
		motifs_and_genes_dataframe = motifs_and_genes_file_or_dataframe

	else: logger.critical('argument not recognized as a file or a dataframe, exiting...'); sys.exit(1)

	return motifs_and_genes_dataframe


######################################### Public Functions #########################################

def map_peaks(peaks_filepath_or_list_of_peaks_filepaths, garnet_filepath):
	"""
	Find motifs and associated genes local to peaks.

	This function intersects peaks from an epigenomics dataset with TF motifs.

	Arguments:
		peaks_filepath_or_list_of_peaks_filepaths (str or list): filepath of the peaks file, or list of such paths
		garnet_filepath (str): filepath to the garnet file.

	Returns:
		pd.DataFrame: a dataframe with rows of transcription factor binding motifs and nearby genes with
		the restriction that these motifs and genes must have been found near a peak.
	"""

	# peaks_filepath_or_list_of_peaks_filepaths is either a filepath or FILE, or a list of filepaths or FILEs.
	# Let's operate on a list in either case, so if it's a single string, put it in a list. #TODO, this will break if it's a single FILE.
	if isinstance(peaks_filepath_or_list_of_peaks_filepaths, list): list_of_peaks_filepaths = peaks_filepath_or_list_of_peaks_filepaths
	else: list_of_peaks_filepaths = [peaks_filepath_or_list_of_peaks_filepaths]
	assert all([os.path.isfile(peaks_filepath) for peaks_filepath in list_of_peaks_filepaths])

	output = []

	for peaks_filepath in list_of_peaks_filepaths:

		peaks = BedTool(peaks_filepath)
		motifs = BedTool(garnet_filepath)

		intersected = motifs.intersect(peaks, wa=True, f=1)
		intersected_df = intersected.to_dataframe(names=["chrom", "start", "end", "motifName", "motifScore", "motifStrand", "geneName", "geneStart", "geneEnd", "motif_gene_distance"])

		output.append(intersected_df)

	# If this function was passed a single file, return a single dataframe
	if len(output) == 1: output = output[0]
	return output


def TF_regression(motifs_and_genes_file_or_dataframe, expression_file, output_dir=None):
	"""
	Do linear regression of the expression of genes versus the strength of the assiciated transcription factor binding motifs and report results.

	This function parses an expression file of two columns: gene symbol and expression value, and
	merges the expression profile into the motifs and genes file, resulting in information about
	transcription factor binding motifs local to genes, and those genes' expressions. We do linear
	regression, and if an output directory is provided, we output a plot for each TF and an html
	summary of the regressions.

	Arguments:
		motifs_and_genes_file_or_dataframe (str or dataframe): the outcome of map_known_genes_and_motifs_to_peaks, either as a dataframe or a file
		expression_file (str or FILE): a tsv file of expression data, with geneName, score columns
		output_dir: (str): If you would like to output figures and a summary html page, supply an output directory

	Returns:
		pd.dataframe: slope, pval, and gene targets for each transcription factor.
	"""

	motifs_and_genes_dataframe = _parse_motifs_and_genes_file_or_dataframe(motifs_and_genes_file_or_dataframe)
	expression_dataframe = parse_expression_file(expression_file)

	motifs_genes_and_expression_levels = motifs_and_genes_dataframe.merge(expression_dataframe, left_on='geneName', right_on='name', how='inner')

	# the same geneName might have different names but since the expression is geneName-wise
	# keep the closest motif to gene in the case of duplicates
	# TODO: implement function to combine duplicates.
	if 'geneName' in motifs_genes_and_expression_levels.columns:
		motifs_genes_and_expression_levels["abs_distance"] = motifs_genes_and_expression_levels.motif_gene_distance.abs()
		motifs_genes_and_expression_levels = motifs_genes_and_expression_levels.sort_values("motifScore", ascending=True) \
		                                                                       .drop("abs_distance", axis=1) \
		                                                                       .drop_duplicates(subset=['geneName', 'motifName'], keep="first")

	motifs_genes_and_expression_levels['motifScore'] = motifs_genes_and_expression_levels['motifScore'].astype(float)

	TFs_and_associated_expression_profiles = list(motifs_genes_and_expression_levels.groupby('motifName'))
	imputed_TF_features = []

	logger.info("Performing linear regression for "+str(len(TFs_and_associated_expression_profiles))+" transcription factor expression profiles...")

	for TF_name, expression_profile in TFs_and_associated_expression_profiles:

		# Occasionally there's only one gene associated with a TF, which we can't fit a line to.
		if len(expression_profile) < 5: continue

		# This allows heavier points to be visualized on the top instead of being hidden by long distance points
		expression_profile = expression_profile.reindex(expression_profile.motif_gene_distance.abs().sort_values(inplace=False, ascending=False).index)

		# Ordinary Least Squares linear regression
		result = linear_regression(formula="expression ~ motifScore", data=expression_profile).fit()

		if output_dir:
			plot = plot_regression(model_results=result, ax=expression_profile.plot(x="motifScore", y="expression", kind="scatter", grid=True))

			# Add color to points based on distance to gene
			plt.scatter(expression_profile["motifScore"], expression_profile["expression"],
				c=[abs(v) for v in expression_profile["motif_gene_distance"].tolist()],
				norm=matplotlib.colors.LogNorm(vmin=1, vmax=100000, clip=True),
				cmap=matplotlib.cm.Blues_r)
			plt.title("%s, %0.4f" %(TF_name, result.pvalues['motifScore']))
			plt.colorbar()

			os.makedirs(os.path.join(output_dir, "regression_plots"), exist_ok=True)
			plot.savefig(os.path.join(output_dir, "regression_plots", TF_name.replace("/", "-") + '.png'))
			plt.close()

		# TODO: implement FDR calculation
		imputed_TF_features.append((TF_name, result.params['motifScore'], result.pvalues['motifScore'], ','.join(expression_profile['geneName'].tolist())))

	imputed_TF_features_dataframe = pd.DataFrame(imputed_TF_features, columns=["Transcription Factor", "Slope", "P-Value", "Targets"]).sort_values("P-Value")

	# If we're supplied with an output_dir, we'll put a summary html file in there as well.
	if output_dir:
		html_output = templateEnv.get_template("summary.jinja").render(images_dir=os.path.join(output_dir,"regression_plots",""), TFs=sorted(imputed_TF_features, key=lambda x: x[2]))
		with open(os.path.join(output_dir,"summary.html"), "w") as summary_output_file:
			summary_output_file.write(html_output)

	return imputed_TF_features_dataframe


######################################## Contruct Garnet File ########################################

def tss_from_bed(bed_file):
	"""
	Write a BED file defining the transcription start sites (TSS) inferred from the
	reference gene file. Sort the TSS and remove any duplicate entries.

	Arguments:
		bed_file (str): path to the reference gene BED file

	Returns:
		str: the path to the output TSS file.
	"""

	output_file = bed_file.replace(".bed", ".tss.bed")

	if os.path.isfile(output_file):
		logger.info('  - TSS file seems to already exists at ' + output_file)
		return output_file

	# These series of commands generates a BED file of TSS from the known genes file.
	# The `awk` command identifies whether the gene is transcribed on the forward or
	# reverse strand. If the gene is on the forward strand, it assumes the TSS is at the
	# start of the region, and if the gene is on the reverse strand, the TSS is placed at
	# the end. The resulting TSS file is sorted and uniquified, and written to output_file.

	# Code taken from https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md
	# TODO: To be replaced by pybedtools command, if possible
	bed_to_tss_cmd = ''' cat %s | awk 'BEGIN{OFS=FS="\t"} \
	       { if ($6 == "+") \
	         { print $1,$2,$2+1,$4,$5,$6 } \
	         else if ($6 == "-") \
	         { print $1,$3-1,$3,$4,$5,$6 } \
	       }' \
	| sort -k1,1 -k2,2n \
	| uniq \
	> %s''' %(bed_file, output_file)

	subprocess.call(bed_to_tss_cmd, shell=True)
	logger.info('  - Wrote TSS file to ' + output_file)

	return output_file


def construct_garnet_file(reference_file, motif_file_or_files, output_file, options=dict()):
	"""
	Construct the GarNet file by searching for any motifs that are within
	a certain window of a reference TSS. Generate a dataframe with these assocations,
	as well as the distance between the motif and gene TSS. Write the dataframe
	to the specified output file.

	Arguments:
		reference_file (str): path to the reference gene BED file
		motifs_file (str): path to the motifs BED file
		output_file (str): ouput GarNet file path
		options (dict): currently not used

	Returns:
		pd.DataFrame: motif-gene associations and distance between the two.
	"""

	# Check whether motif_file_or_files is single path or list of paths
	if isinstance(motif_file_or_files, list): motif_files = motif_file_or_files
	else: motif_files = [motif_file_or_files]
	assert all([os.path.isfile(motif_file) for motif_file in motif_files])

	# Generate BED file of TSS from reference gene file
	reference_tss_file = tss_from_bed(reference_file)
	reference_tss = BedTool(reference_tss_file)

	# Function that searches for all motifs within a window of 10kb from any gene in the
	# reference TSS file, and returns a dataframe.
	# TODO: window size should be generated via the options parameter.
	get_motif_genes_df = lambda motif_file: reference_tss.window(motif_file, w=10000) \
	                                                     .to_dataframe(names=["tssChrom", "tssStart", "tssEnd", "geneName", "tssScore", "tssStrand",
	                                                                          "motifChrom", "motifStart", "motifEnd", "motifName", "motifScore", "motifStrand"])

	# Perform motif-gene matching for each motif file, the concatenate them together.
	logger.info('  - Searching for motifs near genes. This may take a while...')
	motif_genes_df = pd.concat([get_motif_genes_df(motif_file) for motif_file in motif_files])


	# Calculate distance between motif and gene by first identifying the side of the motif
	# that is closest to the gene (this depends on which strand the motif is on), and next
	# finding the distance between the closest end to the gene TSS. A negative value indicates
	# that the motif is upstream of the TSS.
	motif_genes_df["motifClosestEnd"] = motif_genes_df.apply(lambda row: row["motifEnd"] if row["motifStrand"] == "+" else row["motifStart"], axis=1)
	motif_genes_df["motif_gene_distance"] = (motif_genes_df["motifClosestEnd"] - motif_genes_df["tssStart"]) * \
	                                         motif_genes_df.apply(lambda row: 1 if row["tssStrand"] == "+" else -1, axis=1)

	# Filter and reorder columns
	motif_genes_df = motif_genes_df[["motifChrom", "motifStart", "motifEnd", "motifName", "motifScore", "motifStrand",
	                                 "geneName", "tssStart", "tssEnd", "motif_gene_distance"]]

	motif_genes_df.to_csv(output_file, sep='\t', index=False, header=False)
	logger.info('  - %d motif-gene associations found and written to %s' %(motif_genes_df.shape[0], output_file))

	return motif_genes_df

