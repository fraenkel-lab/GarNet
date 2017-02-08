#!/usr/bin/env python3

# Core python modules
import sys
import os

# Peripheral python modules
import argparse
import pickle
import logging

# Core python external libraries
import numpy as np
import pandas as pd
import statsmodels.formula.api as sm

# Peripheral python external libraries
from intervaltree import IntervalTree

# list of public methods:
__all__ = [ "map_known_genes_and_motifs_to_peaks",
			"map_known_genes_to_peaks",
			"map_motifs_to_peaks",
			"map_known_genes_to_motifs",
			"motif_regression"]


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - GarNet: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)


parser = argparse.ArgumentParser(description="""
	Scans genome for nearby features within a given window size.
	If genes and peaks are provided, we map peaks to nearby genes.
	If genes and motif locations are provided, we map motifs to nearby genes.
	If peaks and motif locations are provided, we map motifs to nearby peaks.
	If all three are provided, we map genes to peaks, and map motifs to peaks.
""")

class FullPaths(argparse.Action):
	"""Expand user- and relative-paths"""
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def directory(dirname):
	if not os.path.isdir(dirname): raise argparse.ArgumentTypeError(dirname + " is not a directory")
	else: return dirname

parser.add_argument('-p', '--peaks', dest='peaks_file', type=argparse.FileType('r'),
	help='BED file containing epigenetic regions of interest')  # Add future file formats as we support them
parser.add_argument('-m', '--motifs', dest='motifs_file', type=argparse.FileType('r'),
	help='BED file containing locations and scores of TF motifs')
parser.add_argument('-g', '--genes', dest='known_genes_file', type=argparse.FileType('r'),
	help='file containing locations of known genes in the reference genome (i.e. from UCSC Annotation Database)')
parser.add_argument('-x', '--xref', dest='xref_file', type=argparse.FileType('r'),
	help='file containing information about known genes (i.e. from UCSC Annotation Database)')
parser.add_argument('-e', '--expression', dest='expression_file', type=argparse.FileType('r'),
	help='')

parser.add_argument('--up', dest='upstream_window', type=int, default=2000,
	help='window width in base pairs to consider upstream region [default: %default]')
parser.add_argument('--down', dest='downstream_window', type=int, default=2000,
	help='window width in base pairs to consider downstream region [default: %default]')
parser.add_argument('--tss', dest='tss', action='store_true',
	help='calculate downstream window from transcription start site instead of transcription end site')

parser.add_argument('-o', '--output', dest='output_dir', action=FullPaths, type=directory, required=True,
	help='output directory path')


# if __name__ == '__main__':

# 	args = parser.parse_args()
# 	options = {"upstream_window": args.upstream_window, "downstream_window": args.downstream_window, "tss": args.tss, "output_dir": args.output_dir}

# 	if args.peaks_file and args.motifs_file and args.known_genes_file:
# 		result_dataframe = map_known_genes_and_motifs_to_peaks(args.peaks_file, args.kgXref_file, args.motifs_file, args.known_genes_file, options)
# 		output(result_dataframe, args.output_dir)

# 		if args.expression_file:
# 			output_figs(motif_regression(result_dataframe, args.expression_file, options), args.output_dir)

# 	elif args.peaks_file and args.known_genes_file:
# 		output(map_known_genes_to_peaks(args.peaks_file, args.known_genes_file, options), args.output_dir)

# 	elif args.peaks_file and args.motifs_file:
# 		output(map_motifs_to_peaks(args.peaks_file, args.motifs_file, options), args.output_dir)

# 	elif args.known_genes_file and args.motifs_file:
# 		output(map_known_genes_to_motifs(args.motifs_file, args.known_genes_file, options), args.output_dir)

# 	else: raise InvalidCommandLineArgs()



######################################## File Parsing Logic #######################################

def parse_known_genes_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'r')`)

	The known genes file format is the following:
	http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.sql

	`name` varchar(255) NOT NULL DEFAULT '',
	`chrom` varchar(255) NOT NULL DEFAULT '',
	`strand` char(1) NOT NULL DEFAULT '',
	`txStart` int(10) unsigned NOT NULL DEFAULT '0',
	`txEnd` int(10) unsigned NOT NULL DEFAULT '0',
	`cdsStart` int(10) unsigned NOT NULL DEFAULT '0',
	`cdsEnd` int(10) unsigned NOT NULL DEFAULT '0',
	`exonCount` int(10) unsigned NOT NULL DEFAULT '0',
	`exonStarts` longblob NOT NULL,
	`exonEnds` longblob NOT NULL,
	`proteinID` varchar(40) NOT NULL DEFAULT '',
	`alignID` varchar(255) NOT NULL DEFAULT '',

	Returns:
		dataframe: representation of known genes file
	"""

	known_genes_fieldnames = ["name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","proteinID","alignID"]

	known_genes_dataframe = pd.read_csv(filepath_or_file_object, delimiter='\t', names=known_genes_fieldnames)

	known_genes_dataframe.rename(index=str, columns={"txStart":"geneStart", "txEnd":"geneEnd", "name":"geneName","strand":"geneStrand"}, inplace=True)

	return known_genes_dataframe


def parse_peaks_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'r')`)

	My contract is that I will return to you an instance of Peaks, independent of the filetype you supply me

	BED:

	chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
	chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
	chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.

	The 9 additional optional BED fields are:

	name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
	score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
	strand - Defines the strand. Either "." (=no strand) or "+" or "-".
	thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
	thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
	itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
	blockCount - The number of blocks (exons) in the BED line.
	blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
	blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.


	Returns:
		dataframe: representation of peaks file
	"""

	# if peaks file format is BED

	peaks_fieldnames = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]

	# if peaks file format is MACS

    # peaks_fieldnames = ["chr", "start", "end", "length", "summit", "tags", "-10*log10(pvalue)", "fold_enrichment FDR(%)"]

	# if peaks file format is GEM

	# Location:	the genome coordinate of this binding event
	# IP binding strength:	the number of IP reads associated with the event
	# Control binding strength:	the number of control reads in the corresponding region
	# Fold:	fold enrichment (IP/Control)
	# Expected binding strength:	the number of IP read counts expected in the binding region given its local context (defined by parameter W2 or W3), this is used as the Lambda parameter for the Poisson test
	# Q_-lg10:	-log10(q-value), the q-value after multiple-testing correction, using the larger p-value of Binomial test and Poisson test
	# P_-lg10:	-log10(p-value), the p-value is computed from the Binomial test given the IP and Control read counts (when there are control data)
	# P_poiss:	-log10(p-value), the p-value is computed from the Poission test given the IP and Expected read counts (without considering control data)
	# IPvsEMP:	Shape deviation, the KL divergence of the IP reads from the empirical read distribution (log10(KL)), this is used to filter predicted events given the --sd cutoff (default=-0.40).
	# Noise:	the fraction of the event read count estimated to be noise
	# KmerGroup:	the group of the k-mers associated with this binding event, only the most significant k-mer is shown, the n/n values are the total number of sequence hits of the k-mer group in the positive and negative training sequences (by default total 5000 of each), respectively
	# KG_hgp:	log10(hypergeometric p-value), the significance of enrichment of this k-mer group in the positive vs negative training sequences (by default total 5000 of each), it is the hypergeometric p-value computed using the pos/neg hit counts and total counts
	# Strand:	the sequence strand that contains the k-mer group match, the orientation of the motif is determined during the GEM motif discovery, '*' represents that no k-mer is found to associated with this event

	# peaks_fieldnames = ["Location", "IP binding strength", "Control binding strength", "Fold", "Expected binding strength", "Q_-lg10", "P_-lg10", "P_poiss", "IPvsEMP", "Noise", "KmerGroup", "KG_hgp", "Strand"]


	peaks_dataframe = pd.read_csv(filepath_or_file_object, delimiter='\t', names=peaks_fieldnames)

	peaks_dataframe.rename(index=str, columns={"chromStart":"peakStart", "chromEnd":"peakEnd", "name":"peakName", "score":"peakScore", "strand":"peakStrand"}, inplace=True)

	return peaks_dataframe


def parse_kgXref_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'r')`)

	Returns:
		dataframe: representation of xref file
	"""

	kgXref_fieldnames = ['kgID','mRNA','spID','spDisplayID','geneSymbol','refseq','protAcc','description']

	kgXref_dataframe = pd.read_csv(filepath_or_file_object, delimiter='\t', names=kgXref_fieldnames)

	return kgXref_dataframe


def parse_motifs_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'r')`)

	Returns:
		dataframe: representation of motif file
	"""

	motif_fieldnames = ["chrom", "start", "end", "name", "score", "strand"]

	motif_dataframe = pd.read_csv(filepath_or_file_object, delimiter='\t', names=motif_fieldnames)

	motif_dataframe['motifID'], motif_dataframe['motifName'] = motif_dataframe['name'].str.split('=', 1).str

	motif_dataframe.rename(index=str, columns={"start":"motifStart", "end":"motifEnd", "score":"motifScore", "strand":"motifStrand"}, inplace=True)

	return motif_dataframe


def parse_expression_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'r')`)

	Returns:
		dataframe: representation of expression file
	"""

	return pd.read_csv(filepath_or_file_object, delimiter='\t', names=["name", "expression"])


def save_as_pickled_object(object, filepath): return pickle.dump(object, open(filepath, "wb"))

def load_pickled_object(filepath): return pickle.load(open(filepath, "rb"))

def was_generated_by_pickle(filepath): return True


def output(dataframe, output_dir):
	"""
	Arguments:
		dataframe (dataframe): the principal result of the analysis we want to write out as a csv.
		output_dir (str): the fullpath of a directory we will write our output to.
	"""

	logger.info('Writing output file')

	dataframe.to_csv(output_dir + 'output', sep='\t')

	# finally, under the circumstances that motif_regression was called,
	# 	we'll make a plots dir, plot some shit.
	#	output another dataframe as a CSV with the relevant information


def output_figs(data, output_dir):
	"""
	Arguments:
		dataframe (dataframe): the principal result of the analysis we want to write out as a csv.
		output_dir (str): the fullpath of a directory we will write our output to.
	"""

	pass

	# plots!



######################################### Public Functions #########################################

def map_known_genes_and_motifs_to_peaks(known_genes_file, kgXref_file, motifs_file, peaks_file, options):
	"""
	Arguments:
		peaks_file (str or FILE): filepath or file object for the peaks file.
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		motifs_file (str or FILE): filepath or file object for the motifs file
		options (dict): options which may come from the argument parser.

	Returns:
		dataframe: dataframe
	"""

	peaks = dict_of_IntervalTree_from_peak_file(peaks_file, options['output_dir'])
	reference = dict_of_IntervalTree_from_reference_file(known_genes_file, kgXref_file, options, options['output_dir'])
	motifs = dict_of_IntervalTree_from_motifs_file(motifs_file, options['output_dir'])

	peaks_with_associated_genes_and_motifs = intersection_of_three_dicts_of_intervaltrees(peaks, reference, motifs)

	motifs_and_genes = [{**motif, **gene} for peak, genes, motifs in peaks_with_associated_genes_and_motifs for gene in genes for motif in motifs]

	columns_to_output = ["chrom", "motifStart", "motifEnd", "motifID", "motifName", "motifScore", "geneName", "geneSymbol", "geneStart", "geneEnd"]
	motifs_and_genes = pd.DataFrame.from_records(motifs_and_genes, columns=columns_to_output)

	## compute statistics & add fields

	return motifs_and_genes


def map_known_genes_to_peaks(known_genes_file, peaks_file, kgXref_file, options):
	"""
	Arguments:
		peaks_file (str or FILE): filepath or file object for the peaks file.
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		options (dict): options which may come from the argument parser.

	Returns:
		dataframe: dataframe
	"""

	peaks = dict_of_IntervalTree_from_peak_file(peaks_file, options['output_dir'])
	reference = dict_of_IntervalTree_from_reference_file(known_genes_file, kgXref_file, options, options['output_dir'])

	peaks_with_associated_genes = intersection_of_dict_of_intervaltree(peaks, reference)

	peaks_and_genes = [{**peak, **gene} for peak, gene in peaks_with_associated_genes]

	columns_to_output = ["chrom", "peakStart", "peakEnd", "peakName", "peakScore", "geneName", "geneStart", "geneEnd"]
	peaks_and_genes = pd.DataFrame.from_records(peaks_and_genes, columns=columns_to_output)
	# peaks_and_genes = pd.DataFrame.from_records(peaks_and_genes)

	## compute statistics & add fields

	return peaks_and_genes


def map_motifs_to_peaks(motifs_file, peaks_file, options):
	"""
	Arguments:
		peaks_file (str or FILE): filepath or file object for the peaks file.
		motifs_file (str or FILE): filepath or file object for the motifs file
		options (dict): options which may come from the argument parser.

	Returns:
		dataframe: dataframe
	"""

	peaks = dict_of_IntervalTree_from_peak_file(peaks_file, options['output_dir'])
	motifs = dict_of_IntervalTree_from_motifs_file(motifs_file, options['output_dir'])

	peaks_with_associated_motifs = intersection_of_dict_of_intervaltree(peaks, motifs)

	peaks_and_motifs = [{**peak, **motif} for peak, motif in peaks_with_associated_motifs]

	columns_to_output = ["chrom", "peakStart", "peakEnd", "peakName", "peakScore", "motifID", "motifName", "motifStart", "motifEnd", "motifScore"]
	peaks_and_motifs = pd.DataFrame.from_records(peaks_and_motifs, columns=columns_to_output)

	## compute statistics & add fields

	return peaks_and_motifs


def map_known_genes_to_motifs(known_genes_file, motifs_file, kgXref_file, options):
	"""
	Arguments:
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		motifs_file (str or FILE): filepath or file object for the motifs file
		options (dict): options which may come from the argument parser.

	Returns:
		dataframe: dataframe
	"""

	motifs = dict_of_IntervalTree_from_motifs_file(motifs_file, options['output_dir'])
	reference = dict_of_IntervalTree_from_reference_file(known_genes_file, kgXref_file, options, options['output_dir'])

	motifs_with_associated_genes = intersection_of_dict_of_intervaltree(motifs, reference)

	motifs_and_genes = [{**motif, **gene} for motif, gene in motifs_with_associated_genes]

	columns_to_output = ["chrom", "motifStart", "motifEnd", "motifID", "motifName", "motifScore", "geneName", "geneStart", "geneEnd"]
	motifs_and_genes = pd.DataFrame.from_records(motifs_and_genes, columns=columns_to_output)

	## compute statistics & add fields

	return motifs_and_genes


def motif_regression(motifs_and_genes_dataframe, expression_file):
	"""
	Arguments:

	Returns:

	"""

	expression_dataframe = parse_expression_file(expression_file)

	genes_and_motifs_matrix = construct_matrix(expression_dataframe, motifs_and_genes_dataframe)

	# http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.merge.html

	result = sm.ols(formula="A ~ B", data=genes_and_motifs_matrix).fit()

	print(result.params)
	print(result.summary())



######################################## Private Functions ########################################

def dict_of_IntervalTree_from_peak_file(peaks_file, output_dir):
	"""
	Arguments:
		peaks_file (str or FILE): filepath or file object for the peaks file

	Returns:
		dict: dictionary of intervals in known genes to intervals in peaks.
	"""

	logger.info("Checking if peaks file was generated by pickle, which would imply it's been parsed to IntervalTree already")

	if was_generated_by_pickle(peaks_file):
		logger.info('  - Peaks file seems to have been generated by pickle, assuming IntervalTree format and proceeding...')
		return load_pickled_object(peaks_file)

	logger.info('  - Peaks file does not seem to have been generated by pickle, proceeding to parse...')
	peaks = parse_peaks_file(peaks_file)
	peaks = group_by_chromosome(peaks)
	logger.info('  - Parse complete, constructing IntervalTrees...')
	peaks = {chrom: IntervalTree_from_peaks(chromosome_peaks) for chrom, chromosome_peaks in peaks.items()}

	logger.info('  - IntervalTree construction complete, saving pickle file for next time.')
	save_as_pickled_object(peaks, output_dir + 'peaks_IntervalTree_dictionary.pickle')

	return peaks


def dict_of_IntervalTree_from_reference_file(known_genes_file, kgXref_file, options, output_dir):
	"""
	Arguments:
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		options (dict): options which may come from the argument parser.

	Returns:
		dict: dictionary of chromosome to IntervalTree of known genes
	"""

	logger.info("Checking if known genes file was generated by pickle, which would imply it's been parsed to IntervalTree already")

	if was_generated_by_pickle(known_genes_file):
		logger.info('  - Known Genes file seems to have been generated by pickle, assuming IntervalTree format and proceeding...')
		return load_pickled_object(known_genes_file)

	logger.info('  - Known Genes file does not seem to have been generated by pickle, proceeding to parse...')
	reference = parse_known_genes_file(known_genes_file)

	if kgXref_file:
		logger.info('  - Program was supplied with a kgXref file, merging with reference...')
		kgXref = parse_kgXref_file(kgXref_file)
		reference.merge(kgXref, left_on='geneName', right_on='kgID', how='left')
	else: logger.info('  - Program was not supplied with a kgXref file, gene names will only be supplied as kgID')

	reference = group_by_chromosome(reference)
	logger.info('  - Parse complete, constructing IntervalTrees...')
	reference = {chrom: IntervalTree_from_reference(genes, options) for chrom, genes in reference.items()}

	logger.info('  - IntervalTree construction complete, saving pickle file for next time.')
	save_as_pickled_object(reference, output_dir + 'reference_IntervalTree_dictionary.pickle')

	return reference


def dict_of_IntervalTree_from_motifs_file(motifs_file, output_dir):
	"""
	Arguments:
		motifs_file (str or FILE): filepath or file object for the motifs file

	Returns:
		dict: dictionary of chromosome to IntervalTree of TF binding motifs
	"""

	logger.info("Checking if motifs file was generated by pickle, which would imply it's been parsed to IntervalTree already")

	if was_generated_by_pickle(motifs_file):
		logger.info('  - Motifs file seems to have been generated by pickle, assuming IntervalTree format and proceeding...')
		return load_pickled_object(motifs_file)

	logger.info('  - Motifs file does not seem to have been generated by pickle, proceeding to parse...')
	motifs = parse_motifs_file(motifs_file)
	motifs = group_by_chromosome(motifs)
	logger.info('  - Parse complete, constructing IntervalTrees...')
	motifs = {chrom: IntervalTree_from_motifs(chromosome_motifs) for chrom, chromosome_motifs in motifs.items()}

	logger.info('  - IntervalTree construction complete, saving pickle file for next time.')
	save_as_pickled_object(motifs, output_dir + 'motifs_IntervalTree_dictionary.pickle')

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
		options (dict): options which shall be unpacked here

	Returns:
		IntervalTree: of genes from the reference
	"""

	upstream_window = options['upstream_window']
	downstream_window = options['downstream_window']
	window_ends_downstream_from_transcription_start_site_instead_of_transcription_end_site = options['tss']


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
		motifs (dataframe): Must be a dataframe with start and end columns

	Returns:
		IntervalTree: of motifs
	"""

	intervals = zip(motifs.motifStart.values, motifs.motifEnd.values, motifs.to_dict(orient='records'))

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



def intersection_of_three_dicts_of_intervaltrees(A, B, C):
	"""
	Arguments:
		A (dict): is a dictionary of {chrom (str): IntervalTree}
		B (dict): is a dictionary of {chrom (str): IntervalTree}
		C (dict): is a dictionary of {chrom (str): IntervalTree}

	Returns:
		dict: {keys shared between A, B and C: {intervals in A: [[list of overlapping intervals in B], [list of overlapping intervals in C]]} }
	"""

	logger.info('Computing intersection operation of three IntervalTrees for each chromosome...')

	# See `intersection_of_dict_of_intervaltree` (above) for programmer's notes. This function is nearly identical.
	common_keys = set(A.keys()).intersection( set(B.keys()) ).intersection( set(C.keys()) )

	intersection = [(a.data, [b.data for b in B[key].search(a)],
							 [c.data for c in C[key].search(a)] ) for key in common_keys for a in A[key]]

	return intersection


########################################### Error Logic ###########################################


class Error(Exception):
	def __init__(self, message):
		self.message = message
	def __str__(self):
		return self.message


class InvalidCommandLineArgs(Error):
	def __init__(self, args):
		pass
		# initialize super with a generic message


########################################## Testing Logic ##########################################


peaks = "/Users/alex/Documents/GarNet2/data/A549_FOXA1_broadPeak.bed"
reference = "/Users/alex/Documents/GarNet2/data/ucsc_hg19_knownGenes.tsv"
kgXref = "/Users/alex/Documents/GarNet2/data/ucsc_hg19_kgXref.tsv"
motifs = "/Users/alex/Documents/GarNet2/data/HUMAN_hg19_BBLS_1_00_FDR_0_10.bed"

# peaks = "/Users/alex/Documents/GarNet2/src/peaks_IntervalTree_dictionary.pickle"
# reference = "/Users/alex/Documents/GarNet2/src/reference_IntervalTree_dictionary.pickle"
# motifs = "/Users/alex/Documents/GarNet2/src/motifs_IntervalTree_dictionary.pickle"

output(map_known_genes_and_motifs_to_peaks(reference, kgXref, motifs, peaks, {"upstream_window":2000, "downstream_window":2000, "tss":False, "output_dir":'/Users/alex/Documents/GarNet2/src/'}), '/Users/alex/Documents/GarNet2/src/')


