

# Core python modules
import sys
import os

# Peripheral python external libraries
from pybedtools import BedTool

sys.path.append("../GarNet/")
import garnet




def binned_windows_from_bed(bed_file, window_size, bin_size, genome="hg19"): 

	# Get TSS from reference file
	reference_tss_file = garnet.tss_from_bed(bed_file)

	# Get fixed window around TSS
	tss_windows = BedTool(reference_tss_file).slop(genome=genome, b=window_size)

	# Bin windows
	binned_tss_windows = BedTool().window_maker(b=tss_windows.fn, w=bin_size, i="src")

	return binned_tss_windows



def main():

	binned_tss_windows = binned_windows_from_bed("../garnet_data/reference/ucsc_reference.hg19.bed", 2000, 100)

	intersected = binned_tss_windows.intersect(b=["/nfs/latdata/iamjli/ALS/data/iMNs/ATAC-seq/EdgeR_iMNs_CTR_v_ALS.bed", "/nfs/latdata/iamjli/packages/GarNetMotifs/data_old/cisBP_mapped_motifs_h19/ENSG00000251493_LINE2070_FOXD1_D.bed"],
		names=["peak", "motif"], wa=True, wb=True)

	peaks = binned_tss_windows.intersect(b="/nfs/latdata/iamjli/ALS/data/iMNs/ATAC-seq/EdgeR_iMNs_CTR_v_ALS.bed", c=True)

	print(peaks.head())




main()