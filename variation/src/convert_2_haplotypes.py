"""
Convert a csv binary datafile to a haplotype format.

Usage: python convert_2_haplotype.py data_file haplotype_data_file snps_window

	-w ... 			Window size in SNPs 
"""

import sys
import dataParsers as dp

if __name__ == '__main__':
	#Load Data..
	sd = dp.parse_numerical_snp_data(sys.argv[1])#), filter=0.01)
	#Convert data..
	sd.haplotize(snp_window=int(sys.argv[3]))
	#Save data..
	sd.writeToFile(sys.argv[2])
