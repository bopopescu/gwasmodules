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
import phenotypeData as pd


class PhenotypeInfo(tables.IsDescription):
	"""
	Phenotype info container
	"""
	name = tables.StringCol(256)
	num_values = tables.Int32Col()
	std_dev = tables.Float32Col()
	growth_conditions = tables.StringCol(256)
	phenotype_scoring = tables.StringCol(256)
	method_description = tables.StringCol(256)
	measurement_scale = tables.StringCol(256)



class PhenotypeValue(tables.IsDescription):
	"""
	Phenotype value wrapper
	"""
	ecotype = tables.Int32Col()
	accession_name = tables.StringCol(16)
	mean_value = tables.Float32Col()
	std_dev = tables.Float32Col()
	comment = tables.StringCol(256)


class ResultRecord(tables.IsDescription):
	"""
	A general result record structure
	"""
	chromosome = tables.Int32Col()
	position = tables.Int32Col()
	score = tables.Float32Col() #Perhaps 64 bits?? 
	maf = tables.Float32Col()
	mac = tables.Int32Col()

class ResultRecordLinearModel(ResultRecord):
	"""
	Linear model, mixed models, etc.
	"""
	genotype_var_perc = tables.Float32Col()
	beta0 = tables.Float32Col()
	beta1 = tables.Float32Col()
	correlation = tables.Float32Col()


class ResultRecordKW(ResultRecord):
	"""
	Kruskal Wallis
	"""
	statistic = tables.Float32Col()


class ResultRecordFT(ResultRecord):
	"""
	Fisher's exact test
	"""
	odds_ratio = tables.Float32Col()


def init_file(hdf5_filename):
	print 'Setting up file %s' % hdf5_filename
	# Open a file in "w"rite mode
	h5file = tables.openFile(hdf5_filename, mode="w", title="Phenotype_results_file")
	# Create a new group under "/" (root)
	group = h5file.createGroup("/", 'phenotypes', 'Basic phenotype folder')
	h5file.close()


def add_new_phenotype_file(hdf5_file_name, phenotype_file, phen_name, growth_conditions='', phenotype_scoring='',
			method_description='', measurement_scale='', is_binary=False):
	"""
	Initializes the phenotype group for this phenotype and inserts it into the file object.
	"""
	#Now parsing the phenotype file
	phend = pd.readPhenotypeFile(phenotype_file)
	_init_phenotype_(hdf5_file_name, phen_name, growth_conditions=growth_conditions, phenotype_scoring=phenotype_scoring,
			method_description=method_description, measurement_scale=measurement_scale, is_binary=is_binary)
	add_phenotype_values(hdf5_filename, phen_name, phend.accessions, phend.getPhenVals(1), transformation='raw',
			accessions=phend.accessionNames, std_dev_values=None, value_comments=None)


def _init_phenotype_(hdf5_file_name, phen_name, growth_conditions='', phenotype_scoring='',
			method_description='', measurement_scale='', is_binary=False):
	"""
	Insert a new phenotype into the DB
	"""
	group = h5file.createGroup("/", 'phenotypes', 'Basic phenotype folder')



def add_phenotype_values(hdf5_filename, phen_name, ecotypes, values, transformation='raw', accessions=None,
			std_dev_values=None, value_comments=None):
	"""
	Adds phenotype values.
	"""
	pass



def get_phenotype_values(hdf5_filename, phen_name, ecotypes, transformation='raw'):
	"""
	Returns the phenotype values
	"""
	pass


def get_phenotype_names(hdf5_filename):
	"""
	Returns the phenotype values
	"""
	pass


def get_phenotype_transformations(hdf5_filename, phen_name):
	"""
	Returns the phenotype values
	"""
	pass


def add_results(hdf5_filename, phen_name, chromosomes, positions, scores, mafs, macs, analysis_type, transformation='raw', **kwargs):
	"""
	Add a result to the hdf5 file.
	"""
	pass


#def get_a	


