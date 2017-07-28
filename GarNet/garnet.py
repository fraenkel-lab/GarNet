#!/usr/bin/env python3

# Core python modules
import sys
import os

# Peripheral python modules
import pickle
import logging

# Core python external libraries
import numpy as np
import pandas as pd
from statsmodels.formula.api import ols as linear_regression
from statsmodels.graphics.regressionplots import abline_plot as plot_regression

# Peripheral python external libraries
from intervaltree import IntervalTree
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

def parse_garnet_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or FILE object

	Returns:
		dict: {chr: IntervalTree of regions}
	"""

	garnet_file = try_to_load_as_pickled_object_or_None(filepath_or_file_object)

	if garnet_file == None: sys.exit('Unable to load garnet file')

	return garnet_file


def parse_peaks_file(peaks_file):
	"""
	Parse a BED file with peaks from an epigenomics assay (e.g. ATAC) into a dataframe

	Arguments:
		peaks_file (string or FILE): BED file from epigenomics assay

	Returns:
		dataframe: peaks dataframe
	"""

	peaks_fieldnames = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]

	peaks_dataframe = pd.read_csv(peaks_file, delimiter='\t', names=peaks_fieldnames)

	peaks_dataframe.rename(index=str, columns={"chromStart":"peakStart", "chromEnd":"peakEnd", "name":"peakName", "score":"peakScore", "strand":"peakStrand"}, inplace=True)

	peaks_dataframe = peaks_dataframe[['peakName', 'chrom', 'peakStart', 'peakEnd', 'peakScore']]

	peaks_dataframe[['peakStart','peakEnd']] = peaks_dataframe[['peakStart','peakEnd']].apply(pd.to_numeric)

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

def map_peaks(peaks_file_or_list_of_peaks_files, garnet_file):
	"""
	Find motifs and associated genes local to peaks.

	This function searches for motifs "under" peaks from an epigenomics dataset and "around" peaks for genes.
	It then returns all pairs of motifs and genes which were found local to peaks.

	Arguments:
		garnet_file (str): filepath to the garnet file.
		peaks_file_or_list_of_peaks_files (str or list): filepath of the peaks file, or list of such paths

	Returns:
		pd.dataframe: a dataframe with rows of transcription factor binding motifs and nearby genes with the restriction that these motifs and genes must have been found near a peak.
	"""
	logger.info("Mapping peaks against genome from garnet-file "+garnet_file)
	logger.info("Unpacking garnet file (this can take a while)... ")
	genome = parse_garnet_file(garnet_file)

	# peaks_file_or_list_of_peaks_files is either a filepath or FILE, or a list of filepaths or FILEs.
	# Let's operate on a list in either case, so if it's a single string, put it in a list. #TODO, this will break if it's a single FILE.
	if isinstance(peaks_file_or_list_of_peaks_files, list): list_of_peaks_files = peaks_file_or_list_of_peaks_files
	else: list_of_peaks_files = [peaks_file_or_list_of_peaks_files]

	output = []

	for peaks_file in list_of_peaks_files:

		logger.info("Constructing representation of peaks file "+peaks_file+"...")
		peaks = dict_of_IntervalTree_from_peak_file(peaks_file)

		peaks_with_associated_genomic_regions = intersection_of_dict_of_intervaltree(peaks, genome)

		peak_regions = [{**peak, **region} for peak, region in peaks_with_associated_genomic_regions]

		columns_to_output = ["chrom", "motifStart", "motifEnd", "motifName", "motifScore", "geneName", "geneStart", "geneEnd", "peakName"]
		peak_regions = pd.DataFrame.from_records(peak_regions, columns=columns_to_output)

		peak_regions = peak_regions.apply(type_of_peak, axis=1)

		output.append(peak_regions)

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


######################################## Private Functions ########################################

def construct_garnet_file(known_genes_file, motifs_file, options):
	"""
	Construct a representation of the genome to which to map peaks.

	This function searches for overlap of motifs and genes, and writes a file of motif / gene pairs.

	Arguments:
		known_genes_file (str or FILE): filepath or file object for the known genes file
		motifs_file (str or FILE): filepath or file object for the motifs file
		options (dict): {"upstream_window": int, "downstream_window": int, "tss": bool, "output_dir": string (optional)})
		kgXref_file (str or FILE): filepath or file object for the known genes reference file

	Returns:
		dataframe: A dataframe listing transcription factor binding motifs and nearby genes.
	"""

	reference = dict_of_IntervalTree_from_known_genes_file(known_genes_file, options)
	motifs = dict_of_IntervalTree_from_motifs_file(motifs_file)

	motifs_with_associated_genes = intersection_of_dict_of_intervaltree(motifs, reference)

	motifs_and_genes = [{**motif, **gene} for motif, gene in motifs_with_associated_genes]

	columns_to_output = ["chrom", "motifStart", "motifEnd", "motifName", "motifScore", "geneName", "geneStart", "geneEnd"]
	motifs_and_genes = pd.DataFrame.from_records(motifs_and_genes, columns=columns_to_output)
	motifs_and_genes['motif_to_gene_distance'] = motifs_and_genes['motifStart'] - motifs_and_genes['geneStart']

	return motifs_and_genes


def dict_of_IntervalTree_from_peak_file(peaks_file):
	"""
	Arguments:
		peaks_file (str or FILE): filepath or FILE object

	Returns:
		dict: dictionary of intervals in known genes to intervals in peaks.
	"""

	logger.info('  - Peaks file does not seem to have been generated by pickle, proceeding to parse...')
	peaks = parse_peaks_file(peaks_file)
	peaks = group_by_chromosome(peaks)
	logger.info('  - Parse complete, constructing IntervalTrees...')
	peaks = {chrom: IntervalTree_from_peaks(chromosome_peaks) for chrom, chromosome_peaks in peaks.items()}

	return peaks


def dict_of_IntervalTree_from_known_genes_file(known_genes_file, options):
	"""
	Arguments:
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		options (dict): options which may come from the argument parser.

	Returns:
		dict: dictionary of chromosome to IntervalTree of known genes
	"""

	logger.info('Parsing known genes file(s)...')
	reference = parse_known_genes_file(known_genes_file, kgXref_file=options.get('kgXref_file'))
	reference = group_by_chromosome(reference)
	logger.info('Parse complete, constructing IntervalTrees...')
	reference = {chrom: IntervalTree_from_reference(genes, options) for chrom, genes in reference.items()}

	return reference


def dict_of_IntervalTree_from_motifs_file(motifs_file):
	"""
	Arguments:
		motifs_file (str or FILE): filepath or file object for the motifs file

	Returns:
		dict: dictionary of chromosome to IntervalTree of TF binding motifs
	"""

	logger.info('Parsing motifs file...')
	motifs = parse_motifs_file(motifs_file)
	motifs = group_by_chromosome(motifs)
	logger.info('Parse complete, constructing IntervalTrees...')
	motifs = {chrom: IntervalTree_from_motifs(chromosome_motifs) for chrom, chromosome_motifs in motifs.items()}

	return motifs


def group_by_chromosome(dataframe):
	"""
	Arguments:
		dataframe (dataframe): Must be a dataframe with a chrom column

	Returns:
		dict: dictionary of chromosome names (e.g. 'chr1') to dataframes
	"""

	return dict(list(dataframe.groupby('chrom')))


def IntervalTree_from_peaks(peaks):
	"""
	Arguments:
		peaks (dataframe): Must be a dataframe with peakStart and peakEnd columns

	Returns:
		IntervalTree: of peaks
	"""

	intervals = zip(peaks.peakStart.values, peaks.peakEnd.values, peaks.to_dict(orient='records'))

	tree = IntervalTree.from_tuples(intervals)

	return tree


def IntervalTree_from_reference(reference, options):
	"""
	Arguments:
		reference (dataframe): Must be a dataframe with `strand`, `geneStart`, and `geneEnd` columns
		options (dict): {"upstream_window": int, "downstream_window": int, "tss": bool}

	Returns:
		IntervalTree: of genes from the reference
	"""

	upstream_window = options.get('upstream_window') or 2000
	downstream_window = options.get('downstream_window') or 2000
	window_ends_downstream_from_transcription_start_site_instead_of_transcription_end_site = options.get('tss') or False

	if window_ends_downstream_from_transcription_start_site_instead_of_transcription_end_site:
		starts = reference.apply(lambda x: x.geneStart - upstream_window if x.strand == '+' else x.geneEnd - upstream_window, axis=1)
		ends = reference.apply(lambda x: x.geneStart + downstream_window if x.strand == '+' else x.geneEnd + downstream_window, axis=1)

	else:
		starts = reference.apply(lambda x: x.geneStart - upstream_window, axis=1)
		ends = reference.apply(lambda x: x.geneEnd + downstream_window, axis=1)

	intervals = zip(starts.values, ends.values, reference.to_dict(orient='records'))

	tree = IntervalTree.from_tuples(intervals)

	return tree


def IntervalTree_from_motifs(motifs):
	"""
	Arguments:
		motifs (dataframe): Must be a dataframe with motifStart and motifEnd columns

	Returns:
		IntervalTree: of motifs
	"""

	motifs = motifs.dropna(subset=['motifStart', 'motifEnd'])

	intervals = zip(motifs.motifStart.astype(int).values, motifs.motifEnd.astype(int).values, motifs.to_dict(orient='records'))

	tree = IntervalTree.from_tuples(intervals)

	return tree


def intersection_of_dict_of_intervaltree(A, B):
	"""
	Arguments:
		A (dict): is a dictionary of {chrom (str): IntervalTree}
		B (dict): is a dictionary of {chrom (str): IntervalTree}

	Returns:
		dict: {keys shared between A and B: {intervals in A: [list of overlapping intervals in B]} }
	"""

	logger.info('Computing intersection operation of IntervalTrees for each chromosome...')

	# Keys are chromosomes. We only want to look through chromosomes where there is potential for overlap
	common_keys = set(A.keys()).intersection( set(B.keys()) )

	# In general, the below operation isn't perfectly elegant, due to the double-for:
	#   `for key in common_keys for a in A[key]`
	# The reason we use a double-for here is because we want to do one level of "flattening"
	# but we don't want to do it as a pre-processing or post-processing step. Specifically,
	# we're passed dictionaries of chrom: IntervalTree, and we'd like to return a single data
	# structure for the entire genome, not one per chromosome. The below can be read as:
	#
	#	for chromosome in chromosomes_found_in_both_datastructures:
	#		for a_interval in A[chromosome]
	#			for b_interval in B[chromosome].search(A_interval)
	#				Map a_interval to b_interval
	#
	# which can be expressed concisely in the following line of code (which is also maximally efficient)
	intersection = [(a.data, b.data) for key in common_keys for a in A[key] for b in B[key].search(a)]

	return intersection


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

