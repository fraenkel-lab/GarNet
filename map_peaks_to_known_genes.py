

from csv import DictReader, DictWriter
from optparse import OptionParser



usage = '%prog [options] <knownGene file> <peaks file>'
description = """
Map the peaks in <peaks file> to genes in <knownGene file>.  <knownGene file> is\
format is as specified in http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.sql, though BED format is also accepted.\
<peaks file> format is as produced by GPS, MACS or BED.  If *auto* is chosen (default) file extension \
is examined for *.xls* for default MACS format, *.txt* for GPS, or *.bed* for BED format.  
"""
parser = OptionParser(usage=usage,description=description,epilog='')#,formatter=MultiLineHelpFormatter())
parser.add_option('--upstream-window',dest='upst_win',type='int',default=100000,help='window width in base pairs to consider promoter region [default: %default]')
parser.add_option('--downstream-window',dest='dnst_win',type='int',default=0,help='window width in base pairs to consider downstream region [default: %default]')
parser.add_option('--tss',dest='tss',action='store_true',default=False, help='calculate downstream window from transcription start site instead of transcription end site')
parser.add_option('--map-output',dest='peak_output',default=None,help='filename to output mapped peaks to [default: stdout]')
parser.add_option('--stats-output',dest='stats_output',default=sys.stderr,help='filename to output summary stats in conversion [default: stderr]')
parser.add_option('--peaks-format',dest='peaks_fmt',default='auto',type='choice',choices=['auto','MACS','BED'],help='format of peaks input file [default: %default]')

parser.add_option('--intergenic',dest='intergenic',action='store_true',help='write intergenic peaks to the gene file as well with None as gene ID')
parser.add_option('--utilpath',default='../src/',dest='addpath',help='Destination of chipsequtil library')
parser.add_option('--symbol-xref',dest='symbol_xref',default=None,help='Provide kgXref table file supplied to find a gene symbol and add as second column of output')



def parse_known_genes_file(file):
	"""



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
	"""

	known_genes_file_fieldnames = ["name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","proteinID","alignID"]

	reader = csv.DictReader(file, delimiter='\t', fieldnames=column_names)


def parse_peaks_file(file):
	"""
	My contract is that I will return to you an instance of Peaks, independent of the filetype you supply me
	"""






def parse_kgXref_file(file):
	"""
	"""

	kgXref_file_fieldnames = ['kgID','mRNA','spID','spDisplayID','geneSymbol','refseq','protAcc','description']

	reader = csv.DictReader(file, delimiter='\t', fieldnames=column_names)



def main():
	"""
	"""

	opts, args = parser.parse_args(sys.argv[1:])
	if len(args) != 2: parser.error('Must provide two filename arguments')


	with open(args[0]) as known_genes_file:
		reference = parse_known_genes_file(known_genes_file)

	with open(args[1]) as peaks_file:
		peaks = parse_peaks_file(peaks_file)

	with open(args[2]) as kgXref_file:
		kgXref = parse_kgXref_file(kgXref_file)






if __name__ == '__main__' :
	main()





