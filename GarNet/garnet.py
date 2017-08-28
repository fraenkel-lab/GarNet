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
from statsmodels.formula.api import ols as linear_regression
from statsmodels.graphics.regressionplots import abline_plot as plot_regression

# Peripheral python external libraries
from pybedtools import BedTool
import jinja2

# list of public methods:
__all__ = [ "map_peaks", "TF_regression" ]


templateLoader = jinja2.FileSystemLoader(searchpath=".")
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
		dataframe: expression dataframe
	"""
	df = pd.read_csv(expression_file, delimiter='\t', names=["name", "expression"])
	df["expression"] = pd.to_numeric(df["expression"], errors="coerce")
	df.dropna(inplace=True)

	return df


def parse_known_genes_file(known_genes_file, kgXref_file=None, organism="hg19"):
	"""
	Parse the RefSeq known genes file into a pandas dataframe

	The known genes file format is the following:
	http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.sql

	Arguments:
		known_genes_file (string or FILE): file procured from RefSeq with full list of genes in genome
		kgXref_file (string or FILE): (optional) additional "Cross Reference" file with more details on those genes

	Returns:
		dataframe: known genes dataframe
	"""

	known_genes_fieldnames = ["name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","proteinID","alignID"]

	known_genes_dataframe = pd.read_csv(known_genes_file, delimiter='\t', names=known_genes_fieldnames)

	known_genes_dataframe.rename(index=str, columns={"txStart":"geneStart", "txEnd":"geneEnd", "name":"geneName","strand":"geneStrand"}, inplace=True)

	known_genes_dataframe[['geneStart','geneEnd']] = known_genes_dataframe[['geneStart','geneEnd']].apply(pd.to_numeric)

	if kgXref_file:

		if organism == "hg19":
			kgXref_fieldnames = ["kgID","mRNA","spID","spDisplayID","geneSymbol","refseq","protAcc","description"]

		elif organism == "mm9":
			kgXref_fieldnames = ["kgID","mRNA","spID","spDisplayID","geneSymbol","refseq","protAcc","description"]

		elif organism == "mm10":
			kgXref_fieldnames = ["kgID","mRNA","spID","spDisplayID","geneSymbol","refseq","protAcc","description","rfamAcc","tRnaName"]

		else: logger.critical('organism name entered not recognized in parse_known_genes_file'); sys.exit(1)

		kgXref_dataframe = pd.read_csv(kgXref_file, delimiter='\t', names=kgXref_fieldnames)

		known_genes_dataframe = known_genes_dataframe.merge(kgXref_dataframe, left_on='geneName', right_on='kgID', how='left')
		known_genes_dataframe.rename(index=str, columns={"geneName":"ucID", "geneSymbol":"geneName"}, inplace=True)

	else: logger.info('Program was not supplied with a kgXref file, gene names will only be supplied as kgID')

	return known_genes_dataframe


def parse_motifs_file(motifs_file, organism="hg19"):
	"""
	Parse the MotifMap BED file listing Transcription Factor Binding Motifs in the genome
	Arguments:
		motifs_file (string or FILE): file procured from MotifMap with full list of TF binding sites in the genome
		organism (str): name of the organism, either "hg19", "mm9", or "mm10"
	Returns:
		dataframe: motif dataframe
	"""

	if organism is "mm9":
		motif_fieldnames = ["ZScore","BBLS","FDR","stop","FDR","strand","BLS","accession","FDR","cid","medianhits","start","name","orientation","chrom","stdevhits","LOD","NLOD","realhits"]

	elif organism is "mm10":
		motif_fieldnames = ["stdevhits","ZScore","BLS","name","chrom","FDR_lower","FDR","orientation","start","LOD","cid","strand","realhits","NLOD","BBLS","medianhits","stop","FDR_upper","accession"]

	elif organism is "hg19":
		motif_fieldnames = ["ZScore","FDR_lower","name","orientation","chrom","LOD","strand","start","realhits","cid","FDR","NLOD","BBLS","stop","medianhits","accession","FDR_upper","BLS","stdevhits"]

	else: logger.critical('organism name entered not recognized in parse_motifs_file'); sys.exit(1)

	motif_dataframe = pd.read_csv(motifs_file, delimiter='\t', names=motif_fieldnames)

	motif_dataframe.rename(index=str, columns={"start":"motifStart", "stop":"motifEnd", "FDR":"motifScore", "strand":"motifStrand", "name":"motifName"}, inplace=True)

	motif_dataframe[['motifStart','motifEnd']] = motif_dataframe[['motifStart','motifEnd']].apply(pd.to_numeric)

	return motif_dataframe


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


def save_as_pickled_object(obj, directory, filename):
	"""
	This is a defensive way to write pickle.write, allowing for very large files on all platforms
	"""
	filepath = os.path.join(directory, filename)
	max_bytes = 2**31 - 1
	bytes_out = pickle.dumps(obj)
	n_bytes = sys.getsizeof(bytes_out)
	with open(filepath, 'wb') as f_out:
		for idx in range(0, n_bytes, max_bytes):
			f_out.write(bytes_out[idx:idx+max_bytes])


def try_to_load_as_pickled_object_or_None(filepath):
	"""
	This is a defensive way to write pickle.load, allowing for very large files on all platforms
	"""
	max_bytes = 2**31 - 1
	try:
		input_size = os.path.getsize(filepath)
		bytes_in = bytearray(0)
		with open(filepath, 'rb') as f_in:
			for _ in range(0, input_size, max_bytes):
				bytes_in += f_in.read(max_bytes)
		obj = pickle.loads(bytes_in)
	except:
		return None
	return obj


######################################### Public Functions #########################################

def map_peaks(peaks_filepath_or_list_of_peaks_filepaths, garnet_filepath):
	"""
	Find motifs and associated genes local to peaks.

	This function intersects peaks from an epigenomics dataset with TF motifs. 

	Arguments:
		garnet_filepath (str): filepath to the garnet file.
		peaks_filepath_or_list_of_peaks_filepaths (str or list): filepath of the peaks file, or list of such paths

	Returns:
		ouput (pd.DataFrame): a dataframe with rows of transcription factor binding motifs and nearby genes with 
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

	# conversely, if this function was passed a single file, return a single dataframe
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
	# these additional names cause bogus regression p-values. Get rid of them here.
	if 'geneName' in motifs_genes_and_expression_levels.columns:
		motifs_genes_and_expression_levels.drop_duplicates(subset=['geneName', 'motifName'], inplace=True)
	motifs_genes_and_expression_levels['motifScore'] = motifs_genes_and_expression_levels['motifScore'].astype(float)

	TFs_and_associated_expression_profiles = list(motifs_genes_and_expression_levels.groupby('motifName'))
	imputed_TF_features = []

	logger.info("Performing linear regression for "+str(len(TFs_and_associated_expression_profiles))+" transcription factor expression profiles...")

	for TF_name, expression_profile in TFs_and_associated_expression_profiles:

		# Occasionally there's only one gene associated with a TF, which we can't fit a line to.
		if len(expression_profile) < 5: continue

		# Ordinary Least Squares linear regression
		result = linear_regression(formula="expression ~ motifScore", data=expression_profile).fit()

		if output_dir:
			plot = plot_regression(model_results=result, ax=expression_profile.plot(x="motifScore", y="expression", kind="scatter", grid=True))
			os.makedirs(os.path.join(output_dir, "regression_plots"), exist_ok=True)
			plot.savefig(os.path.join(output_dir, "regression_plots", TF_name + '.jpg'))

		imputed_TF_features.append((TF_name, result.params['motifScore'], result.pvalues['motifScore'], expression_profile['geneName'].tolist()))

	imputed_TF_features_dataframe = pd.DataFrame(imputed_TF_features, columns=["Transcription Factor", "Slope", "P-Value", "Targets"])

	# If we're supplied with an output_dir, we'll put a summary html file in there as well.
	if output_dir:
		html_output = templateEnv.get_template("summary.jinja").render(images_dir=os.path.join(output_dir,"regression_plots",""), TFs=sorted(imputed_TF_features, key=lambda x: x[2]))
		with open(os.path.join(output_dir,"summary.html"), "w") as summary_output_file:
			summary_output_file.write(html_output)

	return imputed_TF_features_dataframe


######################################## Contruct Garnet File ########################################

def tss_from_bed(bed_file): 
	"""
	This function writes a BED file defining the transcription start sites (TSS) inferred from the 
	reference gene file. It also sorts the TSS and removes any duplicate entries. 

	Arguments:
		bed_file (str): path to the reference gene BED file

	Returns:
		output_file (str): the path to the output TSS file. 
	"""

	output_file = bed_file.replace(".bed", ".tss.bed")

	""" These series of commands generates a BED file of TSS from the known genes file. 
	The `awk` command identifies whether the gene is transcribed on the forward or
	reverse strand. If the gene is on the forward strand, it assumes the TSS is at the 
	start of the region, and if the gene is on the reverse strand, the TSS is placed at 
	the end. The resulting TSS file is sorted and uniquified, and written to output_file. 

	Code taken from https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md
	To be replaced by pybedtools command, if possible """
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


def construct_garnet_file(reference_file, motifs_file, output_file, options): 
	"""
	This function constructs the GarNet file by searching for any motifs that are within 
	a certain window of a reference TSS. It generates a dataframe with these assocations, 
	as well as the distance between the motif and gene TSS. It also writes the dataframe
	to the specified output file. 

	Arguments:
		reference_file (str): path to the reference gene BED file
		motifs_file (str): path to the motifs BED file
		output_file (str): ouput GarNet file path

	Returns:
		motif_genes_df (pd.DataFrame): motif-gene associations and distance between the two. 
	"""

	# Generate BED file of TSS from reference gene file
	reference_tss_file = tss_from_bed(reference_file)

	motif = BedTool(motifs_file)
	reference_tss = BedTool(reference_tss_file)

	# Search for all motifs within a window of 10kb from any gene in the reference TSS file
	# TODO: window size should be generated via the options parameter.
	motif_genes = reference_tss.window(motif, w=10000)

	# Generate dataframe from pybedtools object. 
	motif_genes_df = motif_genes.to_dataframe(names=["tssChrom", "tssStart", "tssEnd", "geneName", "tssScore", "tssStrand", 
													 "motifChrom", "motifStart", "motifEnd", "motifName", "motifScore", "motifStrand"])

	""" Calculate distance between motif and gene by first identifying the side of the motif 
	that is closest to the gene (this depends on which strand the motif is on), and next
	finding the distance between the closest end to the gene TSS. A negative value indicates
	that the motif is upstream of the TSS. """
	motif_genes_df["motifClosestEnd"] = motif_genes_df.apply(lambda row: row["motifEnd"] if row["motifStrand"] == "+" else row["motifStart"], axis=1)
	motif_genes_df["motif_gene_distance"] = (motif_genes_df["motifClosestEnd"] - motif_genes_df["tssStart"]) * \
											motif_genes_df.apply(lambda row: 1 if row["tssStrand"] == "+" else -1, axis=1)

	# Filter and reorder columns
	motif_genes_df = motif_genes_df[["motifChrom", "motifStart", "motifEnd", "motifName", "motifScore", "motifStrand",
									 "geneName", "tssStart", "tssEnd", "motif_gene_distance"]]

	motif_genes_df.to_csv(output_file, sep='\t', index=False, header=False)
	logger.info('  - %d motif-gene associations found and written to %s' %(motif_genes_df.shape[0], output_file))

	return motif_genes_df

