"""
A file which contains hdf5 interface for the phenotypes and results.

Overall data structure is: 

One hdf5 file per phenotype.
Three types of tables.
	- A info table, one for each transformation/accession subset.
	- phenotype table, one for each transformation/accession subset
	- result table, one for each transformation/accession subset and one for each analysis method.
	
	The minimum contains an info record, and a raw phenotype table.
"""
import tables


class PhenotypeValue(tables.IsDescription):
	"""
	Phenotype value wrapper
	"""
	ecotype = tables.Int32Col()
	accession_name = tables.StringCol(16)
	measurement_info = tables.StringCol(256)
	mean_value = tables.Float32Col()
	stdev = tables.Float32Col()
	
	
class ResultRecord(tables.IsDescription):
	"""
	"""
	chromosome = tables.Int32Col()
	position = tables.Int32Col()
	score = tables.Float32Col() #Perhaps 64 bits?? 
	

def init_file(filename):
	print 'Setting up file %s'%filename
	# Open a file in "w"rite mode
	h5file = tables.openFile(filename, mode = "w", title = "Phenotype_results_file")
	# Create a new group under "/" (root)
	group = h5file.createGroup("/", 'phenotype', 'Detector information')
	# Create one table on it
	table = h5file.createTable(group, 'readout', PhenotypeValue, "Readout example")
	# Fill the table with 10 particles
	h5file.close()
