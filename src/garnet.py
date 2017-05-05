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
from intervaltree import IntervalTree
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

		motifs_and_genes = [{**peak, **motif_gene_pair} for peak, motif_gene_pairs in peaks_with_associated_genes_and_motifs for motif_gene_pair in motif_gene_pairs]

		columns_to_output = ["chrom", "peakName", "motifStart", "motifEnd", "motifID", "motifName", "motifScore", "geneName", "geneSymbol", "geneStart", "geneEnd"]
		motifs_and_genes = pd.DataFrame.from_records(motifs_and_genes, columns=columns_to_output).set_index("peakName")

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


def query_db(peaks, GarNetDB):
	"""

	Arguments:
		peaks (pd.DataFrame)
		GarNetDB (squlite3.connection):

	"""
	overlaps = []

	df['motifStart']
	df['motifEnd']
	df['chrom']

	peaks.to_sql('peaks', GarNetDB, if_exists="replace")

	overlaps[chrom] = pd.read_sql_query("""
		SELECT garnetdb.*
		FROM garnetdb
		JOIN peaks ON garnetdb.chr        == peaks.chr
				  AND garnetdb.motifStart <= peaks.peakEnd
				  AND garnetdb.motifEnd   >= peaks.peakStart;
		""", GarNetDB)

	return overlaps

