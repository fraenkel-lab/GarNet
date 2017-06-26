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
__all__ = [ "map_known_genes_and_motifs_to_peaks",
			"map_known_genes_to_peaks",
			"map_motifs_to_peaks",
			"map_known_genes_to_motifs",
			"TF_regression",
			"batch_scan_epigenomics_files" ]


templateLoader = jinja2.FileSystemLoader(searchpath=".")
templateEnv = jinja2.Environment(loader=templateLoader)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - GarNet: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)


######################################## File Parsing Logic #######################################

def parse_known_genes_file(known_genes_filepath_or_file_object, kgXref_filepath_or_file_object_or_None):
	"""
	Arguments:
		known_genes_filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'r')`)

	The known genes file format is the following:
	http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.sql

	Returns:
		dataframe: known genes dataframe
	"""

	known_genes_fieldnames = ["name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","proteinID","alignID"]

	known_genes_dataframe = pd.read_csv(known_genes_filepath_or_file_object, delimiter='\t', names=known_genes_fieldnames)

	known_genes_dataframe.rename(index=str, columns={"txStart":"geneStart", "txEnd":"geneEnd", "name":"geneName","strand":"geneStrand"}, inplace=True)

	if kgXref_filepath_or_file_object_or_None == None:
		logger.info('  - Program was not supplied with a kgXref file, gene names will only be supplied as kgID')

	else:
		kgXref_dataframe = parse_kgXref_file(kgXref_filepath_or_file_object_or_None)

		known_genes_dataframe = known_genes_dataframe.merge(kgXref_dataframe, left_on='geneName', right_on='kgID', how='left')

	return known_genes_dataframe


def parse_peaks_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'r')`)

	Returns:
		dataframe: peaks dataframe
	"""

	# if peaks file format is BED

	peaks_fieldnames = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]

	peaks_dataframe = pd.read_csv(filepath_or_file_object, delimiter='\t', names=peaks_fieldnames)

	peaks_dataframe.rename(index=str, columns={"chromStart":"peakStart", "chromEnd":"peakEnd", "name":"peakName", "score":"peakScore", "strand":"peakStrand"}, inplace=True)

	# if peaks file format is MACS

	# peaks_fieldnames = ["chr", "start", "end", "length", "summit", "tags", "-10*log10(pvalue)", "fold_enrichment FDR(%)"]

	# if peaks file format is GEM

	# peaks_fieldnames = ["Location", "IP binding strength", "Control binding strength", "Fold", "Expected binding strength", "Q_-lg10", "P_-lg10", "P_poiss", "IPvsEMP", "Noise", "KmerGroup", "KG_hgp", "Strand"]

	return peaks_dataframe


def parse_kgXref_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'r')`)

	Returns:
		dataframe: additional known genes information as a dataframe
	"""

	kgXref_fieldnames = ["kgID","mRNA","spID","spDisplayID","geneSymbol","refseq","protAcc","description"]

	kgXref_dataframe = pd.read_csv(filepath_or_file_object, delimiter='\t', names=kgXref_fieldnames)

	return kgXref_dataframe


def parse_motifs_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'r')`)

	Returns:
		dataframe: motif dataframe
	"""

	motif_fieldnames = ["ZScore","FDR_lower","name","orientation","chrom","LOD","strand","start","realhits","cid","FDR","NLOD","BBLS","stop","medianhits","accession","FDR_upper","BLS","stdevhits"]
	# motif_fieldnames = ["chrom", "start", "end", "name", "score", "strand"]
	# motif_fieldnames = ["motifName", "chrom", "motifStrand", "motifScore", "motifStart", "motifEnd"]

	motif_dataframe = pd.read_csv(filepath_or_file_object, delimiter='\t', names=motif_fieldnames)

	# motif_dataframe['motifID'], motif_dataframe['motifName'] = motif_dataframe['name'].str.split('=', 1).str

	# motif_dataframe.rename(index=str, columns={"start":"motifStart", "end":"motifEnd", "score":"motifScore", "strand":"motifStrand"}, inplace=True)
	motif_dataframe.rename(index=str, columns={"start":"motifStart", "stop":"motifEnd", "FDR":"motifScore", "strand":"motifStrand", "name":"motifName"}, inplace=True)

	return motif_dataframe


def parse_expression_file(filepath_or_file_object):
	"""
	Arguments:
		filepath_or_file_object (string or FILE): A filepath or file object (conventionally the result of a call to `open(filepath, 'r')`)

	Returns:
		dataframe: expression dataframe
	"""

	return pd.read_csv(filepath_or_file_object, delimiter='\t', names=["name", "expression"])


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


def output(dataframe, output_dir, filename):
	"""
	Arguments:
		dataframe (dataframe): the principal result of the analysis we want to write out as a csv.
		output_dir (str): the fullpath of a directory we will write our output to.
	"""

	logger.info('Writing output file '+filename)

	dataframe.to_csv(os.path.join(output_dir, filename), sep='\t', header=True, index=False)


######################################### Public Functions #########################################

def map_known_genes_and_motifs_to_peaks(known_genes_file, motifs_file, peaks_file, options):
	"""
	Find motifs and associated genes local to peaks.

	This function searches for motifs "under" peaks from an epigenomics dataset and "around" peaks for genes.
	It then returns all pairs of motifs and genes which were found local to peaks.

	Arguments:
		peaks_file (str or FILE): filepath or file object for the peaks file.
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		motifs_file (str or FILE): filepath or file object for the motifs file
		options (dict): {"upstream_window": int, "downstream_window": int, "tss": bool, "output_dir": string (optional)})

	Returns:
		dataframe: a dataframe with rows of transcription factor binding motifs and nearby genes
			with the restriction that these motifs and genes must have been found near a peak.
	"""

	peaks = dict_of_IntervalTree_from_peak_file(peaks_file, options.get('output_dir'))
	reference = dict_of_IntervalTree_from_reference_file(known_genes_file, options, options.get('output_dir'))
	motifs = dict_of_IntervalTree_from_motifs_file(motifs_file, options.get('output_dir'))

	peaks_with_associated_genes_and_motifs = intersection_of_three_dicts_of_intervaltrees(peaks, reference, motifs)

	motifs_and_genes = [{**motif, **gene, **peak} for peak, genes, motifs in peaks_with_associated_genes_and_motifs for gene in genes for motif in motifs]

	columns_to_output = ["chrom", "motifStart", "motifEnd", "motifID", "motifName", "motifScore", "geneName", "geneSymbol", "geneStart", "geneEnd", "peakName"]
	motifs_and_genes = pd.DataFrame.from_records(motifs_and_genes, columns=columns_to_output)

	return motifs_and_genes


def map_known_genes_to_peaks(known_genes_file, peaks_file, options):
	"""
	Find all genes nearby to peaks.

	This function searches in the neighborhood of peaks for genes and returns each peak / gene pair
	which were found to be local to one another.

	Arguments:
		peaks_file (str or FILE): filepath or file object for the peaks file.
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		options (dict): {"upstream_window": int, "downstream_window": int, "tss": bool, "output_dir": string (optional)})

	Returns:
		dataframe: A dataframe listing peaks and nearby genes
	"""

	peaks = dict_of_IntervalTree_from_peak_file(peaks_file, options.get('output_dir'))
	reference = dict_of_IntervalTree_from_reference_file(known_genes_file, options, options.get('output_dir'))

	peaks_with_associated_genes = intersection_of_dict_of_intervaltree(peaks, reference)

	peaks_and_genes = pd.DataFrame.from_records([{**peak, **gene} for peak, gene in peaks_with_associated_genes])

	peaks_and_genes['distance'] = abs((peaks_and_genes['peakStart'] + peaks_and_genes['peakEnd'])/2 - peaks_and_genes['geneStart'])
	peaks_and_genes['type'] = peaks_and_genes.apply(type_of_peak, axis=1)  # upstream/promoter/downstream/intergenic

	columns_to_output = ["chrom", "peakStart", "peakEnd", "peakName", "peakScore", "geneName", "geneStart", "geneEnd", "distance", "type"]

	return peaks_and_genes[columns_to_output]


def map_motifs_to_peaks(motifs_file, peaks_file, options):
	"""
	Find known transcription factor binding motifs motifs "below" input epigenetic peaks.

	This function searches for overlap of motifs and peaks, and returns peak / motif pairs.

	Arguments:
		peaks_file (str or FILE): filepath or file object for the peaks file.
		motifs_file (str or FILE): filepath or file object for the motifs file
		options (dict): {"output_dir": string (optional)})

	Returns:
		dataframe: A dataframe listing peaks and nearby transcription factor binding motifs
	"""

	peaks = dict_of_IntervalTree_from_peak_file(peaks_file, options.get('output_dir'))
	motifs = dict_of_IntervalTree_from_motifs_file(motifs_file, options.get('output_dir'))

	peaks_with_associated_motifs = intersection_of_dict_of_intervaltree(peaks, motifs)

	peaks_and_motifs = [{**peak, **motif} for peak, motif in peaks_with_associated_motifs]

	columns_to_output = ["chrom", "peakStart", "peakEnd", "peakName", "peakScore", "motifID", "motifName", "motifStart", "motifEnd", "motifScore"]
	peaks_and_motifs = pd.DataFrame.from_records(peaks_and_motifs, columns=columns_to_output)

	return peaks_and_motifs


def map_known_genes_to_motifs(known_genes_file, motifs_file, options):
	"""
	Associate genes local to motifs with those motifs, without peak information.

	This function searches for overlap of motifs and genes, and returns motif / gene pairs.

	Arguments:
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		motifs_file (str or FILE): filepath or file object for the motifs file
		options (dict): {"upstream_window": int, "downstream_window": int, "tss": bool, "output_dir": string (optional)})

	Returns:
		dataframe: A dataframe listing transcription factor binding motifs and nearby genes.
	"""

	motifs = dict_of_IntervalTree_from_motifs_file(motifs_file, options.get('output_dir'))
	reference = dict_of_IntervalTree_from_reference_file(known_genes_file, options, options.get('output_dir'))

	motifs_with_associated_genes = intersection_of_dict_of_intervaltree(motifs, reference)

	motifs_and_genes = [{**motif, **gene} for motif, gene in motifs_with_associated_genes]

	columns_to_output = ["chrom", "motifStart", "motifEnd", "motifID", "motifName", "motifScore", "geneName", "geneStart", "geneEnd"]
	motifs_and_genes = pd.DataFrame.from_records(motifs_and_genes, columns=columns_to_output)

	peaks_and_genes['distance'] = peaks_and_genes['motifStart'] - peaks_and_genes['geneStart']

	return motifs_and_genes


def TF_regression(motifs_and_genes_dataframe, expression_file, options):
	"""
	Do linear regression of the expression of genes versus the strength of the assiciated transcription factor binding motifs and report results.

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

	return imputed_TF_features_dataframe


def batch_scan_epigenomics_files(list_of_peaks_files, known_genes_file, motifs_file, options):
	"""
	Scan each peak file for nearby motifs and genes, in the same manner as map_known_genes_and_motifs_to_peaks

	This function principally removes the overhead of loading and unloading and re-loading
	motifs intervaltrees and reference intervaltrees from memory, which is the most time-intensive
	operation in this package. Each of the other functions in this package is targeted towards
	individual samples, but if you receive many peaks files at once, it makes sense to analyze them
	all in one fell swoop.
	This function cannot be called from the CLI defined by the ArgParser above.
	This function writes resulting files to the specified output_dir.

	Arguments:
		list_of_peaks_files (list): a list of filepaths associated with epigenomics files
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		motifs_file (str or FILE): filepath or file object for the motifs file
		options (dict): {"upstream_window": int, "downstream_window": int, "tss": bool, "output_dir": string (optional)})
	"""

	reference = dict_of_IntervalTree_from_reference_file(known_genes_file, options, options.get('output_dir'))
	motifs = dict_of_IntervalTree_from_motifs_file(motifs_file, options.get('output_dir'))

	for peaks_file in list_of_peaks_files:

		logger.info(peaks_file)
		peaks = dict_of_IntervalTree_from_peak_file(peaks_file, None)

		peaks_with_associated_genes_and_motifs = intersection_of_three_dicts_of_intervaltrees(peaks, reference, motifs)

		motifs_and_genes = [{**motif, **gene, **peak} for peak, genes, motifs in peaks_with_associated_genes_and_motifs for gene in genes for motif in motifs]

		columns_to_output = ["chrom", "motifStart", "motifEnd", "motifID", "motifName", "motifScore", "geneName", "geneSymbol", "geneStart", "geneEnd", "peakName"]
		motifs_and_genes = pd.DataFrame.from_records(motifs_and_genes, columns=columns_to_output)

		output(motifs_and_genes, options.get('output_dir'), peaks_file + '.garnet')


######################################## Private Functions ########################################

def dict_of_IntervalTree_from_peak_file(peaks_file, output_dir):
	"""
	Arguments:
		peaks_file (str or FILE): filepath or file object for the peaks file

	Returns:
		dict: dictionary of intervals in known genes to intervals in peaks.
	"""

	logger.info("Checking if peaks file was generated by pickle and trying to load it...")
	peaks = try_to_load_as_pickled_object_or_None(peaks_file)
	if peaks:
		logger.info('  - Peaks file seems to have been generated by pickle, assuming IntervalTree format and proceeding...')
		return peaks

	logger.info('  - Peaks file does not seem to have been generated by pickle, proceeding to parse...')
	peaks = parse_peaks_file(peaks_file)
	peaks = group_by_chromosome(peaks)
	logger.info('  - Parse complete, constructing IntervalTrees...')
	peaks = {chrom: IntervalTree_from_peaks(chromosome_peaks) for chrom, chromosome_peaks in peaks.items()}

	if output_dir:
		logger.info('  - IntervalTree construction complete, saving pickle file for next time.')
		save_as_pickled_object(peaks, output_dir, 'peaks_IntervalTree_dictionary.pickle')

	return peaks


def dict_of_IntervalTree_from_reference_file(known_genes_file, options, output_dir):
	"""
	Arguments:
		known_genes_file (str or FILE): filepath or file object for the known_genes file
		options (dict): options which may come from the argument parser.

	Returns:
		dict: dictionary of chromosome to IntervalTree of known genes
	"""

	logger.info("Checking if known genes file was generated by pickle...")
	reference = try_to_load_as_pickled_object_or_None(known_genes_file)
	if reference:
		logger.info('  - Known Genes file seems to have been generated by pickle, assuming IntervalTree format and proceeding...')
		return reference

	logger.info('  - Known Genes file does not seem to have been generated by pickle, proceeding to parse...')
	reference = parse_known_genes_file(known_genes_file, options.get('kgXref_file'))
	reference = group_by_chromosome(reference)
	logger.info('  - Parse complete, constructing IntervalTrees...')
	reference = {chrom: IntervalTree_from_reference(genes, options) for chrom, genes in reference.items()}

	if output_dir:
		logger.info('  - IntervalTree construction complete, saving pickle file for next time.')
		save_as_pickled_object(reference, output_dir, 'reference_IntervalTree_dictionary.pickle')

	return reference


def dict_of_IntervalTree_from_motifs_file(motifs_file, output_dir):
	"""
	Arguments:
		motifs_file (str or FILE): filepath or file object for the motifs file

	Returns:
		dict: dictionary of chromosome to IntervalTree of TF binding motifs
	"""

	logger.info("Checking if motifs file was generated by pickle...")
	motifs = try_to_load_as_pickled_object_or_None(motifs_file)
	if motifs:
		logger.info('  - Motifs file seems to have been generated by pickle, assuming IntervalTree format and proceeding...')
		return motifs

	logger.info('  - Motifs file does not seem to have been generated by pickle, proceeding to parse...')
	motifs = parse_motifs_file(motifs_file)
	motifs = group_by_chromosome(motifs)
	logger.info('  - Parse complete, constructing IntervalTrees...')
	motifs = {chrom: IntervalTree_from_motifs(chromosome_motifs) for chrom, chromosome_motifs in motifs.items()}

	if output_dir:
		logger.info('  - IntervalTree construction complete, saving pickle file for next time.')
		save_as_pickled_object(motifs, output_dir, 'motifs_IntervalTree_dictionary.pickle')

	return motifs


def group_by_chromosome(dataframe):
	"""
	Arguments:
		dataframe (dataframe): Must be a dataframe with a chrom column

	Returns:
		dict: dictionary of chromosome names (e.g. 'chr1') to dataframes
	"""

	return dict(list(dataframe.groupby('chrom')))


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


