#!/usr/bin/env python3

# Core python modules
import sys

# Peripheral python modules
import argparse
import pickle

# Core python external libraries
import numpy as np
import pandas as pd

# Peripheral python external libraries
from intervaltree import Interval, IntervalTree


parser = argparse.ArgumentParser(description="""
	if genes and peaks are provided, we map peaks to known genes
	if genes and motifs are provided, we map motifs to known genes
	if peaks and motifs are provided, we map motifs to peaks
	if all three are provided, we map known genes to peaks, and map motifs to peaks
""")

parser.add_argument('-p', '--peaks', dest='peaks_file', type=argparse.FileType('r'),
	help='')
parser.add_argument('-m', '--motifs', dest='motifs_file', type=argparse.FileType('r'),
	help='')
parser.add_argument('-g', '--genes', dest='known_genes_file', type=argparse.FileType('r'),
	help='')
parser.add_argument('-x', '--xref', dest='xref_file', type=argparse.FileType('r'),
	help='')

parser.add_argument('--up', dest='upstream_window', type=int, default=2000,
	help='window width in base pairs to consider promoter region [default: %default]')
parser.add_argument('--down', dest='downstream_window', type=int, default=2000,
	help='window width in base pairs to consider downstream region [default: %default]')
parser.add_argument('--tss', dest='tss', action='store_true',
	help='calculate downstream window from transcription start site instead of transcription end site')

parser.add_argument('-o', '--output', dest='output_file', type=argparse.FileType('w'),
	help='')


if __name__ == '__main__':

	args = parser.parse_args(sys.argv[1:])

	if args.peaks_file and args.motifs_file and args.known_genes_file:
		map_known_genes_and_TF_binding_motifs_to_peaks(args.peaks_file, args.motifs_file, args.known_genes_file, args)

	elif args.peaks_file and args.known_genes_file:
		map_known_genes_to_peaks(args.peaks_file, args.known_genes_file, args)

	elif args.peaks_file and args.motifs_file:
		map_motifs_to_peaks(args.peaks_file, args.motifs_file, args)

	elif args.known_genes_file and args.motifs_file:
		map_motifs_to_known_genes(args.known_genes_file, args.motifs_file, args)

	else: raise InvalidUsage("invalid usage")

	if args.ouput_file:
		output(result, output_file)


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

	return known_genes_dataframe


def parse_peaks_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'rb')`)

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

	peaks_dataframe = pd.read_csv(filepath_or_file_object, delimiter='\t', names=peaks_fieldnames)

	# if peaks file format is MACS


	# if peaks file format is GPS


	return peaks_dataframe


def parse_kgXref_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'rb')`)

	Returns:
		dataframe: representation of xref file
	"""

	kgXref_fieldnames = ['kgID','mRNA','spID','spDisplayID','geneSymbol','refseq','protAcc','description']

	kgXref_dataframe = pd.read_csv(filepath_or_file_object, delimiter='\t', names=kgXref_fieldnames)

	return kgXref_dataframe


def parse_motif_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'rb')`)

	Returns:
		dataframe: representation of motif file
	"""

	motif_fieldnames = ["chrom", "start", "end", "name", "score", "strand"]

	motif_dataframe = pd.read_csv(filepath_or_file_object, delimiter='\t', names=motif_fieldnames)

	motif_dataframe["name"] = motif_dataframe["name"].split('=', expand=True)

	# split motif_dataframe.name around = and keep the second half.

	return motif_dataframe


def save(object, filepath): return pickle.dump(object, open(filepath, "wb"))

def load_pickled_object(filepath): return pickle.load(open(filepath, "rb"))

def was_generated_by_pickle(filepath): return False


def output(data, output_file):
	"""
	Arguments:
		options (dict): a filepath might be in here?

	Returns:
		None
	"""

	output_file.write('...')


######################################### Public Functions #########################################

def map_known_genes_and_TF_binding_motifs_to_peaks(peaks_file, motifs_file, known_genes_file, options):
	"""
	Arguments:
		peaks_file (str or FILE): filepath or file object for the peaks file.
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		motifs_file (str or FILE): filepath or file object for the mnotifs file
		options (dict): options which may come from the argument parser.

	Returns:
		dict: dictionary of intervals in peaks to intervals in known genes and motifs.
	"""

	reference = dict_of_IntervalTree_from_reference_file(known_genes_file, options)
	peaks = dict_of_IntervalTree_from_peak_file(peaks_file)
	# kgXref = parse_kgXref_file(kgXref_file)
	motifs = dict_of_IntervalTree_from_motif_file(motifs_file)

	peaks_with_associated_genes_and_motifs = intersection_of_three_dicts_of_intervaltrees(peaks, reference, motifs, options)

	return peaks_with_associated_genes_and_motifs


def map_known_genes_to_peaks(peaks_file, known_genes_file, options): # kgXref_file too?
	"""
	Arguments:
		peaks_file (str or FILE): filepath or file object for the peaks file.
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		options (dict): options which may come from the argument parser.

	Returns:
		dict: dictionary of intervals in peaks to intervals in known genes.
	"""

	reference = dict_of_IntervalTree_from_reference_file(known_genes_file, options)
	peaks = dict_of_IntervalTree_from_peak_file(peaks_file)
	# kgXref = parse_kgXref_file(kgXref_file)

	peaks_with_associated_genes = intersection_of_dict_of_intervaltree(peaks, reference, options)

	return peaks_with_associated_genes


def map_motifs_to_peaks(peaks_file, motifs_file, options):
	"""
	Arguments:
		peaks_file (str or FILE): filepath or file object for the peaks file.
		motifs_file (str or FILE): filepath or file object for the mnotifs file
		options (dict): options which may come from the argument parser.

	Returns:
		dict: dictionary of intervals in peaks to intervals in known genes.
	"""


	peaks = dict_of_IntervalTree_from_peak_file(peaks_file)
	motifs = dict_of_IntervalTree_from_motifs_file(motifs_file)

	peaks_with_associated_motifs = intersection_of_dict_of_intervaltree(peaks, motifs, options)

	return peaks_with_associated_genes


def map_motifs_to_known_genes(known_genes_file, motifs_file, options):  # kgXref_file too?
	"""
	Arguments:
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		motifs_file (str or FILE): filepath or file object for the mnotifs file
		options (dict): options which may come from the argument parser.

	Returns:
		dict: dictionary of intervals in known genes to intervals in motifs.
	"""

	reference = dict_of_IntervalTree_from_reference_file(known_genes_file, options)
	motifs = dict_of_IntervalTree_from_motifs_file(motifs_file)
	# kgXref = parse_kgXref_file(kgXref_file)

	genes_with_associated_motifs = intersection_of_dict_of_intervaltree(reference, motifs, options)

	return genes_with_associated_motifs


######################################## Private Functions ########################################

def dict_of_IntervalTree_from_peak_file(peaks_file):
	"""
	Arguments:
		peaks_file (str or FILE): filepath or file object for the mnotifs file

	Returns:
		dict: dictionary of intervals in known genes to intervals in motifs.
	"""

	if was_generated_by_pickle(peaks_file): return load_pickled_object(peaks_file)

	peaks = parse_peaks_file(peaks_file)
	peaks = group_by_chromosome(peaks)
	peaks = {chrom: IntervalTree_from_peaks(chromosome_peaks) for chrom, chromosome_peaks in peaks}

	return peaks


def dict_of_IntervalTree_from_reference_file(known_genes_file, options):
	"""
	Arguments:
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		options (dict): options which may come from the argument parser.

	Returns:
		dict: dictionary of intervals in known genes to intervals in motifs.
	"""

	if was_generated_by_pickle(known_genes_file): return load_pickled_object(known_genes_file)

	reference = parse_known_genes_file(known_genes_file)
	reference = group_by_chromosome(reference)
	reference = {chrom: IntervalTree_from_peaks(genes, options) for chrom, genes in reference}

	return peaks


def dict_of_IntervalTree_from_motifs_file(motifs_file):
	"""
	Arguments:
		motifs_file (str or FILE): filepath or file object for the mnotifs file

	Returns:
		dict: dictionary of intervals in known genes to intervals in motifs.
	"""

	if was_generated_by_pickle(motifs_file): return load_pickled_object(motifs_file)

	motifs = parse_motifs_file(motifs_file)
	motifs = group_by_chromosome(motifs)
	motifs = {chrom: IntervalTree_from_motifs(chromosome_motifs) for chrom, chromosome_motifs in motifs}

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
		peaks (dataframe): Must be a dataframe with chromStart and chromEnd columns

	Returns:
		IntervalTree: of peaks
	"""

	intervals = zip(peaks.chromStart.values, peaks.chromEnd.values, peaks.values.tolist())

	tree = IntervalTree.from_tuples(intervals)

	return tree


def IntervalTree_from_reference(reference, options):
	"""
	Arguments:
		reference (dataframe): Must be a dataframe with `strand`, `txStart`, and `txEnd` columns
		options (dict): options which shall be unpacked here

	Returns:
		IntervalTree: of genes from the reference
	"""

	upstream_window = options.upstream_window
	downstream_window = options.downstream_window
	window_ends_downstream_from_transcription_start_site_instead_of_transcription_end_site = options.tss


	if window_ends_downstream_from_transcription_start_site_instead_of_transcription_end_site:
		starts = reference.apply(lambda x: x.txStart - upstream_window if x.strand == '+' else x.txEnd - upstream_window, axis=1)
		ends = reference.apply(lambda x: x.txStart + downstream_window if x.strand == '+' else x.txEnd + downstream_window, axis=1)

	else:
		starts = reference.apply(lambda x: x.txStart - upstream_window, axis=1)
		ends = reference.apply(lambda x: x.txEnd + downstream_window, axis=1)


	intervals = zip(starts.values, ends.values, reference.values.tolist())

	tree = IntervalTree.from_tuples(intervals)

	return tree


def IntervalTree_from_motifs(motifs):
	"""
	Arguments:
		motifs (dataframe): Must be a dataframe with start and end columns

	Returns:
		IntervalTree: of motifs
	"""

	intervals = zip(motifs.start.values, motifs.end.values, motifs.values.tolist())

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

	intersection = {}

	for common_key in set(A.keys()).intersection( set(B.keys()) ):

		intersection[common_key] = {a.data: [b.data for b in B[common_key].search(a)] for a in A[common_key]}

	return intersection



def intersection_of_three_dicts_of_intervaltrees(A, B, C):
	"""
	Arguments:
		A (dict): is a dictionary of {chrom (str): IntervalTree}
		B (dict): is a dictionary of {chrom (str): IntervalTree}

	Returns:
		dict: {keys shared between A and B: {intervals in A: [list of overlapping intervals in B]} }
	"""

	intersection = {}

	for common_key in set(A.keys()).intersection( set(B.keys()) ):

		intersection[common_key] = {a.data: [b.data for b in B[common_key].search(a)] for a in A[common_key]}

	return intersection


########################################### Error Logic ###########################################



# peaks = parse_peaks_file("/Users/alex/Documents/GarNet2/data/A549_FOXA1_broadPeak.bed")

# reference = parse_known_genes_file("/Users/alex/Documents/GarNet2/data/ucsc_hg19_knownGenes.txt")

# motifs = parse_motif_file("/Users/alex/Documents/GarNet2/data/HUMAN_hg19_BBLS_1_00_FDR_0_10.bed")




