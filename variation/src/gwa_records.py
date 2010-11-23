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
import itertools
import scipy as sp
import util
import linear_models as lm
import dataParsers as dp
import math
import time
import pdb

chromosome_ends = [30429953, 19701870, 23467451, 18578708, 26992130]

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
	is_binary = tables.BoolCol()



class TransformationInfo(tables.IsDescription):
	"""
	Info on the transformation.
	"""
	name = tables.StringCol(256)
	description = tables.StringCol(256)




class PhenotypeValue(tables.IsDescription):
	"""
	Phenotype value wrapper
	"""
	ecotype = tables.Int32Col()
	accession_name = tables.StringCol(16)
	mean_value = tables.Float32Col()
	std_dev = tables.Float32Col()
	comment = tables.StringCol(256)


class ResultInfo(tables.IsDescription):
	name = tables.StringCol(256)
	comment = tables.StringCol(256)
	max_score = tables.Float64Col()


class ResultRecordLM(tables.IsDescription):
	"""
	Linear model, mixed models, etc.
	"""
	chromosome = tables.Int32Col()
	position = tables.Int32Col()
	score = tables.Float64Col()
	maf = tables.Float32Col()
	mac = tables.Int32Col()
	genotype_var_perc = tables.Float32Col()
	beta0 = tables.Float32Col()
	beta1 = tables.Float32Col()
	correlation = tables.Float32Col()


class ResultRecordKW(tables.IsDescription):
	"""
	Kruskal Wallis
	"""
	chromosome = tables.Int32Col()
	position = tables.Int32Col()
	score = tables.Float64Col()
	maf = tables.Float32Col()
	mac = tables.Int32Col()
	statistic = tables.Float32Col()


class ResultRecordFT(tables.IsDescription):
	"""
	Fisher's exact test
	"""
	chromosome = tables.Int32Col()
	position = tables.Int32Col()
	score = tables.Float64Col()
	maf = tables.Float32Col()
	mac = tables.Int32Col()
	odds_ratio = tables.Float32Col()

class ResultRecordLMStep(tables.IsDescription):
	"""
	Linear model, mixed models, etc. for stepwise analysis.
	"""
	chromosome = tables.Int32Col()
	position = tables.Int32Col()
	score = tables.Float64Col()
	maf = tables.Float32Col()
	mac = tables.Int32Col()
	genotype_var_perc = tables.Float32Col()
	beta0 = tables.Float32Col()
	beta1 = tables.Float32Col()
	correlation = tables.Float32Col()

_file_counter = {}

class GWASRecord():

	def __init__(self, hdf5_file_name):
		self.filename = hdf5_file_name

	def _open(self, mode, *args):
		h5file = tables.openFile(self.filename, mode, *args)
		if self.filename not in _file_counter:
			_file_counter[self.filename] = 1
		else:
			_file_counter[self.filename] += 1
		return h5file

	def _close(self, h5file):
		if self.filename in _file_counter:
			_file_counter[self.filename] -= 1
			if _file_counter[self.filename] == 0:
				print "closing hdf5 file"
				h5file.close()


	def init_file(self):
		h5file = self_open(mode="w", title="Phenotype_results_file")
		g = h5file.createGroup("/", 'phenotypes', 'Basic phenotype folder')
		h5file.createTable(g, 'info', PhenotypeInfo, "Phenotyping information")
		h5file.flush()
		self._close(h5file)

	def add_new_phenotype(self, phen_name, phenotype_values, ecotypes, accession_names=None, growth_conditions='',
				phenotype_scoring='', method_description='', measurement_scale='', is_binary=False):
		"""
		Initializes the phenotype group for this phenotype and inserts it into the file object.
		"""
		#Now parsing the phenotype file
		self.h5file = self._open(mode="r+")
		self._init_phenotype_(phen_name, num_vals=len(phenotype_values), std_dev=sp.std(phenotype_values),
				growth_conditions=growth_conditions, phenotype_scoring=phenotype_scoring,
				method_description=method_description, measurement_scale=measurement_scale,
				is_binary=is_binary)

		self._add_phenotype_values_(phen_name, ecotypes, phenotype_values, transformation='raw',
				accessions=accession_names, std_dev_values=None, value_comments=None)
		self.h5file.flush()
		self._close(self.h5file)


	def _init_phenotype_(self, phen_name, num_vals=0.0, std_dev=0.0, growth_conditions='', phenotype_scoring='',
				method_description='', measurement_scale='', is_binary=False):
		"""
		Insert a new phenotype into the DB
		"""
		group = self.h5file.createGroup("/phenotypes", phen_name, 'Phenotype folder for ' + phen_name)
		table = self.h5file.createTable(group, 'transformation_info', TransformationInfo, "Transformation information")
		table = self.h5file.getNode('/phenotypes/info')
		info = table.row
		info['name'] = phen_name
		info['num_values'] = num_vals
		info['std_dev'] = std_dev
		info['growth_conditions'] = growth_conditions
		info['phenotype_scoring'] = phenotype_scoring
		info['method_description'] = method_description
		info['measurement_scale'] = measurement_scale
		info['is_binary'] = is_binary
		info.append()
		table.flush()



	def add_phenotype_file(self, phen_file_name=None, file_object=None, transformation='raw', transformation_description=None):
		"""
		Adds phenotype values, to an existing phenotype, e.g. when applying different transformations.
		"""
		phen_file_name
		phed = pd.readPhenotypeFile(phen_file_name, file_object, with_db_ids=False)
		phed.filter_na_ecotypes()
		for pid in phed.phenIds:
			(phen_vals, ecotypes) = phed.get_vals_accs(pid)
			self.add_phenotype_values(phed.getPhenotypeName(pid), ecotypes, phen_vals,
						transformation=transformation,
						transformation_description=transformation_description)


	def add_phenotype_values(self, phen_name, ecotypes, values, transformation='raw', transformation_description=None,
				accessions=None, std_dev_values=None, value_comments=None):
		"""
		Adds phenotype values, to an existing phenotype, e.g. when applying different transformations.
		"""
		self.h5file = self._open(mode="r+")
		self._add_phenotype_values_(phen_name, ecotypes, values, transformation, transformation_description,
				accessions, std_dev_values, value_comments)
		self.h5file.flush()
		self._close(self.h5file)


	def delete_transformation(self, phen_name, transformation='raw'):
		"""
		Deletes transformation from an existing phenotype.
		"""
		self.h5file = self._open(mode="r+")
		try:
			try:
				self._delete_transformation_(phen_name, transformation)
			except Exception, err:
				raise err
		finally:
			self.h5file.flush()
			self._close(self.h5file)

	def _delete_transformation_(self, phen_name, transformation='raw'):
		table = self.h5file.getNode('/phenotypes/%s/transformation_info' % phen_name)
		trans_index = -1
		for x in table.iterrows():
			if x['name'] == transformation:
				trans_index = x.nrow
		if trans_index == -1:
			raise Exception('Transformation not found')
		self.h5file.removeNode('/phenotypes/%s' % phen_name, transformation, True)
		table.removeRows(trans_index)
		table.flush()


	def _add_phenotype_values_(self, phen_name, ecotypes, values, transformation='raw', transformation_description=None,
				accessions=None, std_dev_values=None, value_comments=None):
		"""
		"""

		phen_group = self.h5file.getNode('/phenotypes/%s' % phen_name)
		table = self.h5file.getNode('/phenotypes/%s/transformation_info' % phen_name)
		info = table.row
		info['name'] = transformation
		if transformation_description: info['description'] = transformation_description
		info.append()
		table.flush()

		trans_group = self.h5file.createGroup(phen_group, transformation, 'Transformation: ' + transformation)
		table = self.h5file.createTable(trans_group, 'values', PhenotypeValue, "Phenotype values")
		value = table.row
		for i, (ei, v) in enumerate(itertools.izip(ecotypes, values)):
			value['ecotype'] = ei
			value['mean_value'] = v
			if accessions: value['accession_name'] = accessions[i]
			if std_dev_values: value['std_dev'] = std_dev_values[i]
			if value_comments: value['comment'] = value_comments[i]
			value.append()
		table.flush()




	def get_phenotype_values(self, phen_name, transformation='raw'):
		"""
		Returns the phenotype values
		"""
		h5file = self._open(mode="r")
		table = h5file.getNode('/phenotypes/%s/%s/values' % (phen_name, transformation))
		d = {'ecotype' : [], 'mean_value' : [], 'accession_name': [], 'std_dev': [], 'comment':[]}
		for x in table.iterrows():
			for k in d:
				d[k].append(x[k])
		self._close(h5file)
		return d



	def get_phenotype_info(self, phen_name=None):
		"""
		Returns the phenotype meta data in a dict.
		"""
		dict_list = []
		self.h5file = self._open(mode="r")
		table = self.h5file.getNode('/phenotypes/info')
		if not phen_name:
			for x in table.iterrows():
				d = {'name': '', 'num_values': 0, 'std_dev': 0.0, 'growth_conditions': '',
					'phenotype_scoring': '', 'method_description': '', 'measurement_scale': '',
					'is_binary': False}
				for k in d:
					d[k] = x[k]
				d['transformation'] = self._get_phenotype_transformations_(x['name'])
				#for

				dict_list.append(d)
		else:
			for x in table.where('name=="%s"' % phen_name):
				d = {'name': '', 'num_values': 0, 'std_dev': 0.0, 'growth_conditions': '',
					'phenotype_scoring': '', 'method_description': '', 'measurement_scale': '',
					'is_binary': False}
				for k in d:
					d[k] = x[k]
				d['transformation'] = self._get_phenotype_transformations_(x['name'])
				dict_list.append(d)
		self._close(self.h5file)
		return dict_list



	def _get_phenotype_transformations_(self, phen_name):
		dict_list = []
		table = self.h5file.getNode('/phenotypes/%s/transformation_info' % phen_name)
		for x in table.iterrows():
			d = {'name': '', 'description': ''}
			for k in d:
				d[k] = x[k]
			d['phenotype'] = phen_name
			d['analysis_method'] = self._get_analysis_methods_(phen_name, x['name'])
			dict_list.append(d)
		return dict_list



	def get_phenotype_transformations(self, phen_name):
		"""
		Returns the phenotype values
		"""
		self.h5file = self._open(mode="r")
		d = self._get_phenotype_transformations_(phen_name)
		self._close(self.h5file)
		return d



	def _get_analysis_methods_(self, phen_name, transformation):
		dict_list = []
		try:
			table = self.h5file.getNode('/phenotypes/%s/%s/result_info' % (phen_name, transformation))
			for x in table.iterrows():
				d = {'name': '', 'comment': ''}
				for k in d:
					d[k] = x[k]
				d['phenotype'] = phen_name
				d['transformation'] = transformation
				dict_list.append(d)
		except Exception, err_str:
			print "No results found:", err_str
		return dict_list



	def get_analysis_methods(self, phen_name, transformation):
		"""
		Returns the phenotype values
		"""
		self.h5file = self._open(mode="r")
		d = _get_analysis_methods_(phen_name, transformation)
		self._close(self.h5file)
		return d



	def add_results(self, phen_name, analysis_method, chromosomes, positions, scores, mafs, macs,
			analysis_comment='', transformation='raw', **kwargs):
		"""
		Add a result to the hdf5 file.
		"""
		h5file = self._open(mode="r+")
		try:
			trans_group = h5file.getNode('/phenotypes/%s/%s' % (phen_name, transformation))
			table = h5file.createTable(trans_group, 'result_info', ResultInfo, "Result information")
			info = table.row
			info['name'] = analysis_method
			info['max_score'] = max(scores)
			if analysis_comment: info['comment'] = analysis_comment
			info.append()
			table.flush()

			analysis_group = h5file.createGroup(trans_group, analysis_method, 'Analysis method: ' + analysis_method)
			if analysis_method in ['emmax', 'lm']:
				table = h5file.createTable(analysis_group, 'results', ResultRecordLM, "Regression result")
			elif analysis_method == 'kw':
				table = h5file.createTable(analysis_group, 'results', ResultRecordKW, "Regression result")
			else:
				raise Exception('Not implemented for analysis method %s' % analysis_method)
			result = table.row


			for i, cpsmm in enumerate(itertools.izip(chromosomes, positions, scores, mafs, macs)):
				(result['chromosome'], result['position'], result['score'], result['maf'], result['mac']) = cpsmm
				if analysis_method == 'kw':
					result['statistic'] = kwargs['statistics'][i]
				else: #EMMAX or LM
					result['beta0'] = kwargs['beta0'][i]
					result['beta1'] = kwargs['beta1'][i]
					result['correlation'] = kwargs['correlation'][i]
					result['genotype_var_perc'] = kwargs['genotype_var_perc'][i]
				result.append()
		except Exception, err:
			raise err
		finally:
			table.flush()
#		table.cols.chromosome.createIndex()
#		table.cols.score.createIndex()
#		table.cols.mac.createIndex()
			self._close(h5file)




	def get_results(self, phen_name, analysis_method, transformation='raw', min_mac=0, max_pval=1.0):
		"""
		Return results..
		"""
		d = {'chromosome': [], 'position': [], 'score': [], 'maf': [], 'mac': []}
		if analysis_method == 'kw':
			d['statistic'] = []
		else:
			d['beta0'] = []
			d['beta1'] = []
			d['correlation'] = []
			d['genotype_var_perc'] = []

		h5file = self.open(mode="r")
		table = h5file.getNode('/phenotypes/%s/%s/%s/results' % (phen_name, transformation, analysis_method))
		for x in table.where('(score<=%f) & (mac>=%d)' % (max_pval, min_mac)):
			for k in d:
				d[k].append(x[k])
		self._close(h5file)
		return d



	def get_results_by_chromosome(self, phen_name, analysis_method, transformation='raw', min_mac=15, min_score=0.0, \
				top_fraction=0.05, chromosomes=[1, 2, 3, 4, 5], log_transform=True):
		"""
		Return results..
		"""
		cd = {}
		h5file = self._open(mode="r")
		table = h5file.getNode('/phenotypes/%s/%s/%s/results' % (phen_name, transformation, analysis_method))
		info_table = h5file.getNode('/phenotypes/%s/%s/result_info' % (phen_name, transformation))
		for x in info_table.where('name=="%s"' % analysis_method):
			max_score = x['max_score']

		d_keys = ['score', 'position', 'maf', 'mac']
		if analysis_method == 'kw':
			d_keys.append('statistic')
		else:
			d_keys.extend(['beta0', 'beta1', 'correlation', 'genotype_var_perc'])

		for chromosome in chromosomes:
			sort_list = []
			test_chromosomes = []
			for x in table.where('(chromosome==%d) &(score>=%f) & (mac>=%d)' % (chromosome, min_score, min_mac)):
				sort_list.append(tuple([x[k] for k in d_keys]))
				test_chromosomes.append(x['chromosome'])
			print len(sort_list)
			sort_list.sort(reverse=True)
			sort_list = sort_list[:int(top_fraction * len(sort_list))]
			transp_list = map(list, zip(*sort_list))
			d = {}
			for i, k in enumerate(d_keys):
				d[k] = transp_list[i]
			cd[chromosome] = d
		cd['chromosome_ends'] = chromosome_ends
		print '/phenotypes/%s/%s/result_info' % (phen_name, transformation)
		cd['max_score'] = max_score
		self._close(h5file)
		return cd


	def get_phenotype_bins(self, phen_name, transformation='raw', bin_number=20):
		phen_vals = self.get_phenotype_values(phen_name, transformation)['mean_value']
		return self._get_phenotype_bins(phen_vals, bin_number)

	def _get_phenotype_bins(self, phen_vals, bin_number=20):
		"""
		Returns the 
		"""
		#Get the phenotype

		min_phen_val = min(phen_vals)
		max_phen_val = max(phen_vals)
		chunk_size = ((max_phen_val - min_phen_val) / bin_number) * (1 + 1e-5 / bin_number)
		pa = (sp.array(phen_vals) - min_phen_val) / chunk_size
		bin_counts = sp.bincount(sp.array(pa , dtype='int32'))
		keys = ['x_axis', 'frequency']
		bin_list = []
		for i, bin_count in enumerate(bin_counts):
			x1 = min_phen_val + i * chunk_size
			x2 = x1 + chunk_size
			bin_list.append({'x_axis':'%0.2f-%0.2f' % (x1, x2), 'frequency':bin_count})

		return bin_list


	def preview_transform_phenotype(self, phen_name, transformation, original_transformation='raw'):
		new_phen_vals = self.transform_phenotype(phen_name, transformation, original_transformation)
		return self._get_phenotype_bins(new_phen_vals)

	def transform_phenotype(self, phen_name, transformation, original_transformation='raw', store=False):
		"""
		Apply a transformation to a phenotype.
		"""
		phen_data = self.get_phenotype_values(phen_name, original_transformation)
		phen_vals = sp.array(phen_data['mean_value'])
		if transformation == 'raw':
			return
		elif transformation == 'log':
			new_phen_vals = sp.log(phen_vals - min(phen_vals) + sp.var(phen_vals) * 0.1)
		elif transformation == 'sqrt':
			new_phen_vals = sp.sqrt(phen_vals - min(phen_vals) + sp.var(phen_vals) * 0.1)
		if store:
			self.add_phenotype_values(phen_name, phen_data['ecotype'], new_phen_vals, transformation)
		return new_phen_vals

	def exists_transformation(self, phen_name, transformation):
		transformations = self.get_phenotype_transformations(phen_name)
		for trans in transformations:
			if trans['name'] == transformation:
				return True
		return False



	def perform_gwas(self, phen_name, transformation='raw', analysis_method='kw', \
			snps_data_file='/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t54.csv.binary',
			kinship_file='/Users/bjarnivilhjalmsson/Projects/Data/250k/kinship_matrix_cm54.pickled'):

		"""
		Performs GWAS and updates the datastructure.
		"""
		import bisect
		import gwa
		if analysis_method not in ['lm', 'emmax', 'kw']:
			raise Exception('analysis method %s not supported' % analysis_method)
		phen_data = self.get_phenotype_values(phen_name, transformation)
		sd = dp.parse_numerical_snp_data(snps_data_file)
		sd.filter_accessions(map(str, phen_data['ecotype']))
		phen_vals = []
		for ei in sd.accessions:
			i = bisect.bisect(phen_data['ecotype'], int(ei)) - 1
			phen_vals.append(phen_data['mean_value'][i])
		snps = sd.getSnps()
		positions = sd.getPositions()
		chromosomes = []
		for s, c in itertools.izip(sd.snpsDataList, sd.chromosomes):
			chromosomes.extend([c] * len(s.snps))
		maf_dict = sd.get_mafs()
		correlations = gwa.calc_correlations(snps, phen_vals)


		kwargs = {}
		if analysis_method == 'emmax':
			k = lm.load_kinship_from_file(kinship_file, sd.accessions)
			res = lm.emmax(snps, phen_vals, k)
		elif analysis_method == 'lm':
			res = lm.linear_model(snps, phen_vals)
		elif analysis_method == 'kw':
			kw_res = util.kruskal_wallis(snps, phen_vals)
			scores = map(lambda x:-math.log10(x), kw_res['ps'])
			self.add_results(phen_name, analysis_method, chromosomes, positions, scores, maf_dict['marfs'],
					maf_dict['mafs'], transformation=transformation, statistics=kw_res['ds'],
					correlation=correlations)
		else:
			raise Exception('analysis method %s not supported' % analysis_method)

		if analysis_method in ['lm', 'emmax']:
			betas = map(list, zip(*res['betas']))
			scores = map(lambda x:-math.log10(x), res['ps'])
			self.add_results(phen_name, analysis_method, chromosomes, positions, scores, maf_dict['marfs'],
					maf_dict['mafs'], transformation=transformation,
					genotype_var_perc=res['var_perc'], beta0=betas[0], beta1=betas[1],
					correlation=correlations)
		print 'Done!'


	def conditional_analysis(self):
		"""
		Conditions on the selected SNP, for either EMMAX, or LM.
		"""



#	def add_new_phenotype_file(hdf5_file_name, phenotype_file, phen_name, growth_conditions='', phenotype_scoring='',
#				method_description='', measurement_scale='', is_binary=False):
#		"""
#		Initializes the phenotype group for this phenotype and inserts it into the file object.
#		"""
#		#Now parsing the phenotype file
#		h5file = tables.openFile(self.filename, mode="r+")
#		print h5file
#		phend = pd.readPhenotypeFile(phenotype_file)
#		_init_phenotype_(h5file, phen_name, growth_conditions=growth_conditions, phenotype_scoring=phenotype_scoring,
#				method_description=method_description, measurement_scale=measurement_scale, is_binary=is_binary)
#		add_phenotype_values(h5file, phen_name, phend.accessions, phend.getPhenVals(1), transformation='raw',
#				accessions=phend.accessionNames, std_dev_values=None, value_comments=None)
#		h5file.flush()
#		h5file.close()
#


def _test_():
	#Load phenotype data..
	import phenotypeData as pd
	import gwaResults as gr
	phed = pd.readPhenotypeFile('/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_raw_092910.tsv')
	pid1 = 1
	phed.filter_accessions_w_missing_data(pid1)
	phen_name = phed.getPhenotypeName(pid1)
	phen_vals = phed.getPhenVals(pid1)
	ecotypes = phed.accessions
	is_binary = phed.isBinary(pid1)

	#Creating the first hdf5 file
	hdf5_file_name_1 = '/Users/bjarni.vilhjalmsson/tmp/test1.hdf5'
	gwa_record = GWASRecord(hdf5_file_name_1)
	gwa_record.init_file()
	gwa_record.add_new_phenotype(phen_name, phen_vals, ecotypes, is_binary=is_binary)
	print "First file is constructed"

	print "Now testing it"
	r = gwa_record.get_phenotype_values(phen_name, 'raw')
	#print r

	phed = pd.readPhenotypeFile('/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_raw_092910.tsv')
	pid2 = 5
	phed.filter_accessions_w_missing_data(pid2)
	phen_name = phed.getPhenotypeName(pid2)
	phen_vals = phed.getPhenVals(pid2)
	ecotypes = phed.accessions
	is_binary = phed.isBinary(pid2)
	gwa_record.add_new_phenotype(phen_name, phen_vals, ecotypes, is_binary=is_binary)

	print "Now testing it"
	r = gwa_record.get_phenotype_values(phen_name, 'raw')
	#print r
	r = gwa_record.get_phenotype_info(phen_name)
	print r

	gwa_record.transform_phenotype('FT10', transformation='sqrt')
	print "Now testing it"
	r = gwa_record.get_phenotype_values(phen_name, 'raw')
	#print r
	r = gwa_record.get_phenotype_info(phen_name)
	print r

	result_file = '/Users/bjarnivilhjalmsson/tmp/pi1_pid5_FT10_emmax_none.pvals'
	res = gr.Result(result_file=result_file, name='FT10')
	res.neg_log_trans()

#	for c in ['chromosomes', 'positions', 'scores', 'marfs', 'mafs', 'genotype_var_perc', 'beta0', \
#		'beta1', 'correlations']:
#		print c, res.snp_results[c][:10]


	gwa_record.add_results(phen_name, 'emmax', res.snp_results['chromosomes'], res.snp_results['positions'],
			res.scores, res.snp_results['marfs'], res.snp_results['mafs'],
			transformation='raw', genotype_var_perc=res.snp_results['genotype_var_perc'],
			beta0=res.snp_results['beta0'], beta1=res.snp_results['beta1'],
			correlation=res.snp_results['correlations'])


	print "Result added."

	print "Now fetching a result."
	res = gwa_record.get_results(phen_name, 'emmax')#, min_mac=15, max_pval=0.01)
	print "Result loaded"
#	for c in ['chromosome', 'position', 'score', 'maf', 'mac', 'genotype_var_perc', 'beta0', \
#		'beta1', 'correlation']:
#		print c, res[c][:10]
	r = gwa_record.get_phenotype_info()
	print r
	s1 = time.time()
 	res = gwa_record.get_results_by_chromosome(phen_name, 'emmax')
 	print "Result re-loaded"
       	secs = time.time() - s1
        if secs > 60:
                mins = int(secs) / 60
                secs = secs - mins * 60
                print 'Took %d mins and %f seconds.' % (mins, secs)
        else:
                print 'Took %f seconds.' % (secs)
	for chromosome in [1, 2, 3, 4, 5]:
		for c in ['position', 'score', 'maf', 'mac', 'genotype_var_perc', 'beta0', \
			'beta1', 'correlation']:
			print c, res[chromosome][c][:10]
	print res['chromosome_ends']
	print res['max_score']
	print gwa_record.get_phenotype_bins(phen_name)
	s1 = time.time()
	gwa_record.perform_gwas('LD', analysis_method='kw')
      	secs = time.time() - s1
        if secs > 60:
                mins = int(secs) / 60
                secs = secs - mins * 60
                print 'Took %d mins and %f seconds.' % (mins, secs)
        else:
                print 'Took %f seconds.' % (secs)

	gwa_record.transform_phenotype('LD', transformation='log')
	gwa_record.perform_gwas('LD', analysis_method='emmax', transformation='log')



if __name__ == '__main__':
	_test_()



