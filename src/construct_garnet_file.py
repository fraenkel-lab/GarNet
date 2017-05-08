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

# Peripheral python external libraries
from intervaltree import IntervalTree

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - GarNet: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)


######################################## File Parsing Logic #######################################

def parse_known_genes_file(known_genes_file, kgXref_file=None):
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

	if kgXref_file:

		kgXref_fieldnames = ["kgID","mRNA","spID","spDisplayID","geneSymbol","refseq","protAcc","description"]
		kgXref_dataframe = pd.read_csv(kgXref_file, delimiter='\t', names=kgXref_fieldnames)

		known_genes_dataframe = known_genes_dataframe.merge(kgXref_dataframe, left_on='geneName', right_on='kgID', how='left')
		known_genes_dataframe.rename(index=str, columns={"geneName":"ucID", "geneSymbol":"geneName"}, inplace=True)

	else: logger.info('Program was not supplied with a kgXref file, gene names will only be supplied as kgID')

	return known_genes_dataframe


def parse_motifs_file(motifs_file):
	"""
	Parse the MotifMap BED file listing Transcription Factor Binding Motifs in the genome

	Arguments:
		motifs_file (string or FILE): file procured from MotifMap with full list of TF binding sites in the genome

	Returns:
		dataframe: motif dataframe
	"""

	motif_fieldnames = ["ZScore","FDR_lower","name","orientation","chrom","LOD","strand","start","realhits","cid","FDR","NLOD","BBLS","stop","medianhits","accession","FDR_upper","BLS","stdevhits"]
	# motif_fieldnames = ["chrom", "start", "end", "name", "score", "strand"]
	# motif_fieldnames = ["motifName", "chrom", "motifStrand", "motifScore", "motifStart", "motifEnd"]

	motif_dataframe = pd.read_csv(motifs_file, delimiter='\t', names=motif_fieldnames)

	# motif_dataframe['motifID'], motif_dataframe['motifName'] = motif_dataframe['name'].str.split('=', 1).str

	# motif_dataframe.rename(index=str, columns={"start":"motifStart", "end":"motifEnd", "score":"motifScore", "strand":"motifStrand"}, inplace=True)
	motif_dataframe.rename(index=str, columns={"start":"motifStart", "stop":"motifEnd", "FDR":"motifScore", "strand":"motifStrand", "name":"motifName"}, inplace=True)

	return motif_dataframe


######################################### Public Functions #########################################

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

	columns_to_output = ["chrom", "motifStart", "motifEnd", "motifID", "motifName", "motifScore", "geneName", "geneStart", "geneEnd"]
	motifs_and_genes = pd.DataFrame.from_records(motifs_and_genes, columns=columns_to_output)
	motifs_and_genes['motif_to_gene_distance'] = motifs_and_genes['motifStart'] - motifs_and_genes['geneStart']

	logger.info('Writing GarNetDB file to disk...')

	GarNetDB = sqlite3.connect(os.path.join(options['output_dir'], 'garnetDB.sql'))
	motifs_and_genes.to_sql('garnetdb', GarNetDB, if_exists="replace")

	GarNetDB.execute("CREATE INDEX chr_start_stop on garnetdb(chrom, motifStart, motifEnd);")

	GarNetDB.commit()
	GarNetDB.close()

	return motifs_and_genes


######################################## Private Functions ########################################

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
