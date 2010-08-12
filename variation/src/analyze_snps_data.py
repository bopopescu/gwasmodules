"""
Contains basic functions to analyze SNP data, e.g. LD levels, polymorphism rates, Fst, recombination, etc.
"""

def calc_r2_levels( graph_file=None):
	"""
	Returns statistics on LD levels, and can plot them.
	"""
	import dataParsers
	snps_data = dataParsers.parse_snp_data("/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t43_192.csv")
	snps_data.convert_2_binary()
	print "Converted to binary"
	snps_data.snpsDataList[0].calc_r2_levels()


def calc_fst_levels(snps_data, bin_size=100000, graph_file=None):
	pass


if __name__ == "__main__":
	calc_r2_levels()