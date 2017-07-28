#!/usr/bin/env python3

# Core python modules
import sys

# Peripheral python modules
import argparse

# import this module
from GarNet import map_peaks, TF_regression


parser = argparse.ArgumentParser(prog="GarNet", description="""
	Scans genome for features in the garnet file local to peaks from the peaks file.
	Performs univariate regression of TF binding affinity to RNA expression if RNA data is supplied.
	If intermediate file and expression are supplied, only performs TF regression step, assumes peaks already mapped.
""", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

class FullPaths(argparse.Action):
	"""Expand user- and relative-paths"""
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def directory(dirname):
	if not os.path.isdir(dirname): raise argparse.ArgumentTypeError(dirname + " is not a directory")
	else: return dirname

parser.add_argument('-p', '--peaks', dest='peaks_file', type=argparse.FileType('r'), required=False,
	help='BED file containing genomic regions of interest')  # Add future file formats as we support them
parser.add_argument('-g', '--garnet', dest='garnet_file', type=argparse.FileType('r'), required=False,
	help='A garnet-generated file of of the genome to search against')
parser.add_argument('-e', '--expression', dest='expression_file', type=argparse.FileType('r'), required=False,
	help='a two-column tab-delimited file of genes and (relative) expression levels')
parser.add_argument('-i', '--intermediate', dest='intermediate_file', type=argparse.FileType('r'), required=False,
	help='a .garnet.tsv file output by this program with mapped peaks')

parser.add_argument('-o', '--output', dest='output_dir', action=FullPaths, type=directory, required=True,
	help='output directory path')


def output(dataframe, output_dir, filename):
	"""
	Arguments:
		dataframe (dataframe): the principal result of the analysis we want to write out as a csv.
		output_dir (str): the fullpath of a directory we will write our output to.
	"""
	dataframe.to_csv(os.path.join(output_dir, filename), sep='\t', header=True, index=False)


def main():

	args = parser.parse_args()

	# "Map Peaks"
	if args.peaks_file and args.garnet_file:
		result_dataframe = map_peaks(args.peaks_file, args.garnet_file)
		output(result_dataframe, args.output_dir, args.peaks_file+'.garnet.tsv')

		# "Map Peaks + TF Regression"
		if args.expression_file:
			output(TF_regression(result_dataframe, args.expression_file), args.output_dir, args.expression_file+'.prizes.tsv')

	# "TF Regression"
	elif args.expression_file and args.intermediate_file:
		output(TF_regression(args.intermediate_file, args.expression_file), args.output_dir, args.expression_file+'.prizes.tsv')

	else: raise Exception('Improper parameters. See GarNet --help for more.')


if __name__ == '__main__':
	main()
