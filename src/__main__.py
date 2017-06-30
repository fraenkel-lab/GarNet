#!/usr/bin/env python3

# Core python modules
import sys

# Peripheral python modules
import argparse

# import this module
from . import map_peaks, TF_regression


parser = argparse.ArgumentParser(prog="GarNet", description="""
	Scans genome for features in the garnet file local to peaks from the peaks file.
""", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

class FullPaths(argparse.Action):
	"""Expand user- and relative-paths"""
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def directory(dirname):
	if not os.path.isdir(dirname): raise argparse.ArgumentTypeError(dirname + " is not a directory")
	else: return dirname

parser.add_argument('-p', '--peaks', dest='peaks_file', type=argparse.FileType('r'),
	help='BED file containing genomic regions of interest')  # Add future file formats as we support them
parser.add_argument('-g', '--garnet', dest='garnet_file', type=argparse.FileType('r'),
	help='A garnet-generated file of of the genome to search against')
parser.add_argument('-e', '--expression', dest='expression_file', type=argparse.FileType('r'),
	help='a two-column tab-delimited file of genes and (relative) expression levels')

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

	if args.peaks_file and args.garnet_file:
		result_dataframe = map_peaks(args.peaks_file, args.garnet_file)
		output(result_dataframe, args.output_dir, args.peaks_file+'.garnet.tsv')

		if args.expression_file:
			output(TF_regression(result_dataframe, args.expression_file, options), args.output_dir, args.expression_file+'.prizes.tsv')

	else: raise Exception('GarNet requires a peaks (foreground) file and a garnet (background) file. GarNet --help for more.')



if __name__ == '__main__':
	main()
