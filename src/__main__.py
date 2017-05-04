#!/usr/bin/env python3

# Core python modules
import sys

# Peripheral python modules
import argparse

# import this module
from . import construct_garnet_file, map_peaks, TF_regression


parser = argparse.ArgumentParser(prog="GarNet", description="""


""", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

class FullPaths(argparse.Action):
	"""Expand user- and relative-paths"""
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def directory(dirname):
	if not os.path.isdir(dirname): raise argparse.ArgumentTypeError(dirname + " is not a directory")
	else: return dirname

parser.add_argument('-p', '--peaks', dest='peaks_file', type=argparse.FileType('r'),
	help='BED file containing epigenetic regions of interest')  # Add future file formats as we support them
parser.add_argument('-db', '--garnetdb', dest='garnet_file', type=argparse.FileType('r'),
	help='GarNetDB file containing Reference Genome and Transcription Factor Binding Motifs')
parser.add_argument('-e', '--expression', dest='expression_file', type=argparse.FileType('r'),
	help='Gene Expression file (tsv format) for TF regression')

parser.add_argument('--up', dest='upstream_window', type=int, default=2000,
	help='window width in base pairs to consider upstream region')
parser.add_argument('--down', dest='downstream_window', type=int, default=2000,
	help='window width in base pairs to consider downstream region')
parser.add_argument('--tss', dest='tss', action='store_true',
	help='calculate downstream window from transcription start site instead of transcription end site')

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
	options = {"upstream_window": args.upstream_window, "downstream_window": args.downstream_window, "tss": args.tss, "kgXref_file": args.kgXref_file, "output_dir": args.output_dir}

	if args.peaks_file and args.garnet_file:
		result_dataframe = map_peaks(args.peaks_file, args.garnet_file, options)
		output(result_dataframe, args.output_dir, args.peaks_file+'.garnet')

		if args.expression_file:
			output(TF_regression(result_dataframe, args.expression_file, options), args.output_dir, args.expression_file+'.prizes')

	else: raise Exception('GarNet requires at least two files, some combination of [known genes, motifs, peaks] or all three. GarNet --help for more.')


if __name__ == '__main__':
	main()
