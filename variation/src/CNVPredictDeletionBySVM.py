#!/usr/bin/env python
"""

Examples:
	#2010-7-1 use CNV method 8 data as training to predict
	CNVPredictDeletionBySVM.py -u yh
	-i ~/mnt/hpc-cmb/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_2Flanking_GADA_A0.5T8.0M5.tsv -t 8 -m 10 -y 1 -c
	
	# 2010-7-16 for TAIR9 in banyan's mysql db
	CNVPredictDeletionBySVM.py
	-i /usr/local/home_ubuntu/crocea/mnt/hpc-cmb/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_TAIR9_4Flanking_GADA_A0.5T12.0M5.tsv
	-t 7 -m 8 -y 1 -c -u yh -z banyan.usc.edu
	
	# 2010-7-25 for TAIR9 in banyan's mysql db. 
	CNVPredictDeletionBySVM.py
	-i /usr/local/home_ubuntu/crocea/mnt/hpc-cmb/script/variation/data/CNV/call_53_WABQN_b100_j50_lts_TAIR9_4Flanking_GADA_A0.5T12.0M5.tsv
	-t 10 -m 11 -y 1 -c -u yh -z banyan.usc.edu -f 0.6 -l 2 -x 100
	
Description:
	2010-7-1 This program predicts deletions based on training data fetched from table CNVCall via a SVM model.
		SVM model involves 3 features: amplitude (mean intensity of a segment), log10(no_of_probes), probeDensity per kb.
		
		Final data gets inserted into db, CNVCall under a different cnv_method_id.
	
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, numpy, getopt
from pymodule import figureOutDelimiter, getListOutOfStr, getColName2IndexFromHeader, PassingData
import Stock_250kDB

class CNVPredictDeletionBySVM(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
						('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('input_fname', 1, ): ['', 'i', 1, 'input file. output by MpiGADA.py', ],\
						('training_cnv_method_id', 1, int): [8, 't', 1, 'CNV method id of the training data in CNVCall'],\
						('cnv_method_id', 1, int): [10, 'm', 1, 'CNV method id for the predicted deletions to be stored in CNVCall.'],\
						('cnv_type_id', 1, int): [1, 'y', 1, 'CNV type id. table CNVType', ],\
						('max_amplitude', 1, float): [-0.095, 'a', 1, 'only segments whose amplitude below this gets into the prediction'],\
						('minPercUnCoveredByLerContig', 1,  float): [0.6, 'f', 1, 'minimum LerContig-uncovered fraction for a segment in the training data to be regarded as true deletion'],\
						('max_median_intensity_dist', 1, float): [250, 'x', 1, 'maximum distance for a 2nd model to be included'],\
						('min_array_median_intensity', 1, float): [250, 'n', 1, 'minimum median intensity for one array to be included'],\
						('minNoOfModelArrays', 1, int): [5, '', 1, 'minimum number of model arrays for one array to ensure consensus'],\
						('SVM_C', 1, float): [10, 'C', 1, 'The C parameter for SVM'],\
						('SVM_gamma', 0, float): [0., 'g', 1, 'gamm used in SVM'],\
						('SVM_eps', 0, float): [1e-2, 'e', 1, 'epsilon, when to stop the iteration..', ],\
						('deletedFractionType', 1, float): [1, 'l', 1, 'deleted fraction type: 1 (percUnCoveredByLerContig), 2 (fractionDeletedInPECoverageData)'],\
						('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self, **keywords):
		"""
		2009-10-28
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
		#self.arrays_to_form_model = getListOutOfStr(self.arrays_to_form_model, data_type=int)
	
	def get_array_id2median_intensity(self, min_array_median_intensity=250):
		"""
		2010-7-1
		"""
		sys.stderr.write("Getting  array median intensity ...")
		array_id2median_intensity = {}
		query = Stock_250kDB.ArrayInfo.query.filter(Stock_250kDB.ArrayInfo.median_intensity>=min_array_median_intensity)
		for row in query:
			array_id = row.id
			array_id2median_intensity[array_id] = row.median_intensity
		sys.stderr.write("%s arrays.\n"%(len(array_id2median_intensity)))
		return array_id2median_intensity
	
	def constructSVMModels(self, db_250k, arrays_to_form_model, array_id2median_intensity,\
						minPercUnCoveredByLerContig=0.6, cnv_method_id=6, kernel_type=None, C=10, gamma=0., \
						eps=1e-2, deletedFractionType=1):
		"""
		2010-7-25
			add argument deletedFractionType
				1: CNVCall.percUnCoveredByLerContig
				2: CNVCall.fractionDeletedInPECoverageData
		2010-7-1
		"""
		sys.stderr.write("Constructing SVM models for %s arrays ...\n"%(len(arrays_to_form_model)))
		from svm import svm_problem, svm_parameter, svm_model, cross_validation, LINEAR, POLY, RBF
		if kernel_type is None:
			kernel_type = RBF
		param = svm_parameter(C = C, eps=eps, probability = 1, gamma=gamma, kernel_type = kernel_type)
		array_id2model = {}
		for array_id in arrays_to_form_model:
			if array_id not in array_id2median_intensity:	#model array has to be in array_id2median_intensity
				continue
			cnvFeatureData = self.getCNVFeatureData(db_250k, array_id=array_id, \
					minPercUnCoveredByLerContig=minPercUnCoveredByLerContig, cnv_method_id=cnv_method_id, \
					replaceAmpWithMedianIntensity=False, deletedFractionType=deletedFractionType)
			
			problem = svm_problem(cnvFeatureData.class_label_ls, cnvFeatureData.feature_data)
			model = svm_model(problem, param)
			array_id2model[array_id] = model
		sys.stderr.write("%s models.\n"%(len(array_id2model)))
		return array_id2model
	
	def mapAnyArray2ModelArray(self, array_id2median_intensity, array_id2model, max_median_intensity_dist=150,\
							minNoOfModelArrays=5):
		"""
		2010-7-1
		"""
		sys.stderr.write("Mapping any array to a model array ...")
		array_id2model_array_id_ls = {}
		for array_id, median_intensity in array_id2median_intensity.iteritems():

			median_intensity_dist_array_id_ls = []
			for model_array_id in array_id2model:
				if model_array_id != array_id:	#not to predict itself
					model_array_median = array_id2median_intensity.get(model_array_id)
					median_intensity_dist = abs(model_array_median - median_intensity)
					median_intensity_dist_array_id_ls.append((median_intensity_dist, model_array_id))
			median_intensity_dist_array_id_ls.sort()
			no_of_model_arrays = 0
			for i in xrange(len(median_intensity_dist_array_id_ls)):
				if median_intensity_dist_array_id_ls[i][0]<=max_median_intensity_dist:
					no_of_model_arrays = i+1
			no_of_model_arrays = max(minNoOfModelArrays, no_of_model_arrays)	#minimum to ensure consensus
			
			if array_id not in array_id2model_array_id_ls:
				array_id2model_array_id_ls[array_id] = []
			for i in xrange(no_of_model_arrays):
				array_id2model_array_id_ls[array_id].append(median_intensity_dist_array_id_ls[i][1])
			
			"""
			# 2010-7-25
			# old way of doing it since only 4 Ler arrays qualified as model arrays
			if array_id in array_id2model:
				array_id2model_array_id_ls[array_id] = array_id
			else:
				median_intensity_dist_array_id_ls = []
				for model_array_id in array_id2model:
					model_array_median = array_id2median_intensity.get(model_array_id)
					median_intensity_dist = abs(model_array_median - median_intensity)
					median_intensity_dist_array_id_ls.append((median_intensity_dist, model_array_id))
				median_intensity_dist_array_id_ls.sort()
				model_array_id = median_intensity_dist_array_id_ls[0][1]	#the closest array
				array_id2model_array_id_ls[array_id] = model_array_id
			"""
		sys.stderr.write("%s arrays found their respective model array. Done.\n"%(len(array_id2model_array_id_ls)))
		return array_id2model_array_id_ls
	
	
	
	@classmethod
	def getCNVFeatureData(cls,  db_250k, array_id=None, \
					minPercUnCoveredByLerContig=0.6, cnv_method_id=6, \
					replaceAmpWithMedianIntensity=False, deletedFractionType=1):
		"""
		2010-7-25
			add argument deletedFractionType
				1: CNVCall.percUnCoveredByLerContig
				2: CNVCall.fractionDeletedInPECoverageData
		2010-7-1
			moved from CNV.CNVPredictionBySVM in misc.py
		"""
		sys.stderr.write("Getting CNV feature data (amplitude, #probes, probe density,) array %s, cnv_method %s, minPercUnCoveredByLerContig %s ... \n"%\
						(array_id, cnv_method_id, minPercUnCoveredByLerContig))
		i = 0
		block_size = 5000
		real_counter = 0
		TableClass = Stock_250kDB.CNVCall
		query = TableClass.query.filter_by(array_id=array_id).filter_by(cnv_method_id=cnv_method_id)
		rows = query.offset(i).limit(block_size)
		session = db_250k.session
		
		ecotype_id = None
		percUnCoveredByLerContig_ls = []
		feature_data = []
		class_label_ls = []
		c_ls = []
		while rows.count()!=0:
			for row in rows:
				ecotype_id = row.array.maternal_ecotype_id
				if deletedFractionType==1:
					deletedFraction = row.percUnCoveredByLerContig
				else:
					deletedFraction = row.fractionDeletedInPECoverageData
				if deletedFraction is not None:
					#x_ls.append(row.amplitude)
					no_of_probes = math.log10(row.no_of_probes_covered)
					probeDensity = row.no_of_probes_covered*1000.0/(row.stop-row.start+1.0)
					if deletedFraction>=minPercUnCoveredByLerContig:
						class_label = -1
						real_counter += 1
					else:
						class_label = 1
					class_label_ls.append(class_label)
					if replaceAmpWithMedianIntensity:
						amp = row.median_intensity
					else:
						amp = row.amplitude
					feature_data.append([amp, no_of_probes, probeDensity ])
					percUnCoveredByLerContig_ls.append(deletedFraction)
				
				i += 1
			if i%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, i, real_counter))
			rows = query.offset(i).limit(block_size)
		sys.stderr.write("%s%s\t%s\n"%('\x08'*80, i, real_counter))
		return PassingData(feature_data=feature_data, class_label_ls=class_label_ls, \
						percUnCoveredByLerContig_ls=percUnCoveredByLerContig_ls, ecotype_id=ecotype_id)
	
	def predictOneRow(self, cnv_segment_obj, model):
		"""
		2010-7-1
		"""
		probeDensity = cnv_segment_obj.no_of_probes*1000.0/cnv_segment_obj.segment_length
		one_feature_data = [cnv_segment_obj.amplitude, math.log10(cnv_segment_obj.no_of_probes), probeDensity]
		prediction, prediction2probability = model.predict_probability(one_feature_data)
		return prediction, prediction2probability
	
	@classmethod
	def saveSegmentObj(cls, param_obj, cnv_segment_obj):
		"""
		2010-7-1
			moved from CNV.putOneCNVSegmentIntoDBFunctor() from misc.py
		"""
		if not hasattr(param_obj, 'no_of_total'):
			setattr(param_obj, 'no_of_total', 0)
		if not hasattr(param_obj, 'no_of_into_db'):
			setattr(param_obj, 'no_of_into_db', 0)
		session = getattr(param_obj, 'session', None)
		if session is None:
			sys.stderr.write("Error: db session is not available.\n")
			return
		param_obj.no_of_total += 1
		
		cnv_method_id = getattr(param_obj, "cnv_method_id", None)
		cnv_type_id = getattr(param_obj, "cnv_type_id", None)
		array_id = cnv_segment_obj.array_id
		no_of_probes = cnv_segment_obj.no_of_probes
		rows = Stock_250kDB.CNVCall.query.filter_by(array_id=array_id).filter_by(start_probe_id=cnv_segment_obj.start_probe_id).\
				filter_by(stop_probe_id = cnv_segment_obj.stop_probe_id).filter_by(cnv_type_id = cnv_type_id).\
				filter_by(cnv_method_id=cnv_method_id)
		if rows.count()==0:	#make sure it's not in db yet.
			probability = getattr(cnv_segment_obj, 'probability', None)
			cnv_call = Stock_250kDB.CNVCall(array_id=array_id, chromosome = cnv_segment_obj.segment_chromosome,\
										start = cnv_segment_obj.segment_start_pos, stop = cnv_segment_obj.segment_stop_pos, \
										start_probe_id = cnv_segment_obj.start_probe_id, \
										stop_probe_id = cnv_segment_obj.stop_probe_id, \
										no_of_probes_covered = cnv_segment_obj.no_of_probes, \
										size_affected = cnv_segment_obj.segment_length,\
										amplitude = cnv_segment_obj.amplitude,\
										median_intensity = cnv_segment_obj.median_intensity,\
										probability = probability)
			cnv_call.cnv_method_id = cnv_method_id
			cnv_call.cnv_type_id = cnv_type_id
			cnv_call.comment = getattr(cnv_segment_obj, 'comment', None)	#2010-7-25
			session.add(cnv_call)
			param_obj.no_of_into_db += 1
		if param_obj.no_of_into_db%10000==0:
			session.flush()
			session.expunge_all()
		
		report = getattr(param_obj, 'report', False)
		if report and param_obj.no_of_total%10000==0:
			sys.stderr.write('%s%s into db out of %s\n'%('\x08'*100, param_obj.no_of_into_db, param_obj.no_of_total))
	
	def predictOneSegmentByMultipleModels(self, cnv_segment_obj, model_array_id_ls, array_id2model):
		"""
		2010-7-25
		"""
		label_predicted2count = {}
		label_probability_ls = []
		for model_array_id in model_array_id_ls:
			model = array_id2model.get(model_array_id)
			label_predicted, label_predicted2probability = self.predictOneRow(cnv_segment_obj, model)
			if label_predicted not in label_predicted2count:
				label_predicted2count[label_predicted] = 0
			label_predicted2count[label_predicted] += 1
			try:
				label_probability_ls.append(label_predicted2probability[-1])	#only care about deletion's probability
			except:
				print cnv_segment_obj
				print "model_array_id, label_predicted, label_predicted2probability:"
				print model_array_id, label_predicted, label_predicted2probability
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
		
		
		count_label_predicted_ls = []
		for label_predicted, count in label_predicted2count.iteritems():
			count_label_predicted_ls.append((count, label_predicted))
		count_label_predicted_ls.sort()
		count_label_predicted_ls.reverse()	#the label with bigger counts is in the front.
		
		if len(count_label_predicted_ls)>=2 and count_label_predicted_ls[0][0] == count_label_predicted_ls[1][0]:
			label_predicted = -1	#if two labels have same number of arrays supported, choose deletion.
		else:
			label_predicted = count_label_predicted_ls[0][1]
		
		# only care about deletion probability
		label_predicted2probability = {}
		label_predicted2probability[-1] = numpy.median(label_probability_ls)
		return label_predicted, label_predicted2probability
	
	def predictALLSegments(self, input_fname, array_id2model_array_id_ls, array_id2model,\
						max_amplitude=-0.1, param_obj=None):
		"""
		2010-7-25
			handle the situation that any arrays has >=3 model-arrays
		2010-7-1
		"""
		sys.stderr.write('Predicting for all segments from %s ... \n'%(input_fname))
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		median_col_index = col_name2index.get('median')
		ecotype_id_idx = col_name2index.get('ecotype_id', col_name2index.get('array_id'))
		counter = 0
		no_of_segments_in_model = 0
		no_of_predicted_deletions = 0
		for row in reader:
			counter += 1
			amplitude = float(row[col_name2index['amplitude']])
			if amplitude>max_amplitude:
				continue
			cnv_ecotype_id = int(row[ecotype_id_idx])
			array_id = int(row[col_name2index.get('array_id')])
			if array_id not in array_id2model_array_id_ls:
				continue
			no_of_probes = int(row[col_name2index['length']])
			
			start_probe = row[col_name2index['start_probe']].split('_')	# split chr_pos
			start_probe = map(int, start_probe)
			start_probe_id = row[col_name2index['start_probe_id']]
			stop_probe = row[col_name2index['end_probe']].split('_')
			stop_probe = map(int, stop_probe)
			stop_probe_id = row[col_name2index['end_probe_id']]
			
			segment_chromosome = start_probe[0]
			if start_probe[0]!=stop_probe[0]:	#spurious. on different chromosomes.
				continue
			segment_start_pos = start_probe[1]-12
			segment_stop_pos = stop_probe[1]+12
			segment_length = abs(segment_stop_pos-segment_start_pos+1)
			
			if median_col_index is not None:
				median_intensity = float(row[median_col_index])
			else:
				median_intensity = None
			cnv_segment_obj = PassingData(ecotype_id=cnv_ecotype_id, start_probe=start_probe, stop_probe=stop_probe,\
												no_of_probes=no_of_probes, amplitude=amplitude, segment_length=segment_length,\
												segment_chromosome=segment_chromosome, array_id=array_id,\
												start_probe_id=start_probe_id, stop_probe_id=stop_probe_id,\
												segment_start_pos=segment_start_pos, segment_stop_pos=segment_stop_pos,\
												median_intensity=median_intensity)
			model_array_id_ls = array_id2model_array_id_ls.get(array_id)
			no_of_segments_in_model += 1
			label_predicted, label_predicted2probability = self.predictOneSegmentByMultipleModels(cnv_segment_obj, \
																	model_array_id_ls, array_id2model)
			if label_predicted==-1:	# predicted to be deletion.
				cnv_segment_obj.probability = label_predicted2probability[-1]
				cnv_segment_obj.comment = 'model arrays: %s'%(repr(model_array_id_ls)[1:-1])
				self.saveSegmentObj(param_obj, cnv_segment_obj)
				no_of_predicted_deletions += 1
			if no_of_predicted_deletions%5000==0:
				sys.stderr.write('%s%s\t%s\t%s'%('\x08'*100, counter, no_of_segments_in_model, no_of_predicted_deletions))
		sys.stderr.write('%s%s\t%s\t%s\n'%('\x08'*100, counter, no_of_segments_in_model, no_of_predicted_deletions))
		sys.stderr.write('%s out of %s segments were used in prediction. %s predicted deletions.\n'%\
						(no_of_segments_in_model, counter, no_of_predicted_deletions))
		
	def getModelArrays(self, db_250k, training_cnv_method_id, array_id2median_intensity):
		"""
		2010-7-29
			bug-fix: require fractionDeletedInPECoverageData to be not null.
		2010-7-25
			
		"""
		sys.stderr.write("Getting model arrays ...")
		rows = db_250k.metadata.bind.execute("select v.* from view_array v, (select distinct array_id from %s \
				where cnv_method_id=%s and fractionDeletedInPECoverageData is not null) nt where nt.array_id=v.array_id order by median_intensity"%\
				(Stock_250kDB.CNVCall.table.name, training_cnv_method_id))
		arrays_to_form_model = []
		for row in rows:
			if row.array_id in array_id2median_intensity:
				arrays_to_form_model.append(row.array_id)
		sys.stderr.write(" %s arrays.\n"%(len(arrays_to_form_model)))
		return arrays_to_form_model
		
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				  			password=self.db_passwd, hostname=self.hostname, database=self.dbname, 
				   			schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		array_id2median_intensity = self.get_array_id2median_intensity(min_array_median_intensity=self.min_array_median_intensity)
		arrays_to_form_model = self.getModelArrays(db, self.training_cnv_method_id, array_id2median_intensity)
		if self.debug:	# 2010-7-25 for debug, temporary
			arrays_to_form_model = arrays_to_form_model[:4]
		
		array_id2model = self.constructSVMModels(db, arrays_to_form_model, array_id2median_intensity,\
						minPercUnCoveredByLerContig=self.minPercUnCoveredByLerContig, cnv_method_id=self.training_cnv_method_id,\
						C=self.SVM_C, gamma=self.SVM_gamma, eps=self.SVM_eps, deletedFractionType=self.deletedFractionType)
		
		array_id2model_array_id_ls = self.mapAnyArray2ModelArray(array_id2median_intensity, array_id2model, \
															max_median_intensity_dist=self.max_median_intensity_dist,\
															minNoOfModelArrays=self.minNoOfModelArrays)
		param_obj = PassingData(session=session, no_of_total=0, no_of_into_db=0, report=self.report,\
							cnv_method_id=self.cnv_method_id, cnv_type_id=self.cnv_type_id)
		
		self.predictALLSegments(self.input_fname, array_id2model_array_id_ls, array_id2model,\
						max_amplitude=self.max_amplitude, param_obj=param_obj)
		session.flush()
		session.expunge_all()
		session.commit()
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CNVPredictDeletionBySVM
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()