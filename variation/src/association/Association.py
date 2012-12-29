#!/usr/bin/env python
"""

Examples:
	#try kruskal wallis
	%s -i /tmp/250K_method_5_after_imputation_noRedundant_051908.tsv -P /Network/Data/250k/finalData_051808/phenotypes.tsv -E -r -o /tmp/250K_method_5_after_imputation_noRedundant_051908.LD.pvalue
	
	#try linear model
	%s -i /Network/Data/250k/tmp-yh/call_method_17_test.tsv -P ./banyan_fs/tmp/phenotype.tsv -o /tmp/call_method_17_lm.tsv -y2
	
	#try emma (linear mixture model) on 1st 7 phenotypes
	%s -i ./mnt2/panfs/250k/call_method_17.tsv -P ./banyan_fs/tmp/phenotype.tsv -o /tmp/call_method_17_y3.tsv  -y3 -w 1-7

	#linear model with principal components 0 to 9, phenotype from 1 to 7
	%s -i /Network/Data/250k/tmp-yh/call_method_17.tsv -P /Network/Data/250k/tmp-yh/phenotype.tsv -y4 -o /Network/Data/250k/tmp-yh/eigenstrat//call_method_17_lm_with_pc0_9 -W 0-9 -f /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.pca.evec -r -w 1-7
	
	#linear model with PCs 0 to 1, phenotype from 1 to 5. the PCs are calculated on the fly according to the snp input file.
	%s -i /Network/Data/250k/tmp-yh/250k_data/call_method_17_chr4_100000_700000.tsv -P /Network/Data/250k/tmp-yh/phenotype.tsv -o /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_chr4_100000_700000_y4_pc0_1 -W 0-1 -y 4 -w 1-5 -r
	
	#y = SNP + environment + noise, for phenotype 1 & 2
	%s -i /Network/Data/250k/tmp-yh/call_method_17.tsv -P /Network/Data/250k/tmp-yh/phenotype.tsv -o /Network/Data/250k/tmp-yh/association_results/lm/call_method_17_y5.tsv -y5 -w 1,2 -r
	
	#y = SNP + environment + PC1 + PC2 + noise, for phenotype 1 & 2
	%s -i /Network/Data/250k/tmp-yh/call_method_17.tsv -P /Network/Data/250k/tmp-yh/phenotype.tsv -o /Network/Data/250k/tmp-yh/association_results/lm_with_PC12/call_method_17_y5.tsv -W 0-1 -f /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.pca.evec -y5 -w 1,2 -r
	
	#y = SNP + environment + SNP X environ + noise, for phenotype 1 & 2
	%s -i /Network/Data/250k/tmp-yh/call_method_17.tsv -P /Network/Data/250k/tmp-yh/phenotype.tsv -o /Network/Data/250k/tmp-yh/association_results/lm/call_method_17_y6.tsv -y6 -w 1,2 -r
	
	#y = SNP + environment + SNP X environ + PC1 + PC2 + noise, for phenotype 1 & 2
	%s -i /Network/Data/250k/tmp-yh/call_method_17.tsv -P /Network/Data/250k/tmp-yh/phenotype.tsv -o /Network/Data/250k/tmp-yh/association_results/lm_with_PC12/call_method_17_y6.tsv -W 0-1 -f /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.pca.evec -y6 -w 1,2 -r
	
	# 2010-2-1 EMMAX
	%s -i /Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv -P /Network/Data/250k/tmp-yh//phenotype.tsv -o /tmp/call_method_17_y8.tsv  -y8 -w 1
	
	#2010-8-7 Run KW on binary (0-1) CNV deletion data. -n is added to turn off snpAlleleOrdinalConversion.
	%s -i ~/panfs/250k/CNV/NonOverlapCNVAsSNP_cnvMethod20.tsv -P ~/panfs/250k/phenotype/phenotype.tsv
		-o ~/panfs/250k/association_results/cnvMethod20/cnvMethod20_y1_pheno.tsv -y1 -w 1-7 -n
	
	# 2011-4-27 run EMMAX on cnv-turned-into-SNP dataset with precomputed kinship matrix (same format as snp data in inputFname).
	# "-n" is essential because in this CNV-SNP dataset, 0 is normal, not NA; 1 is deletion.
	%s -i /Network/Data/250k/db/dataset/call_method_57.tsv -K ~/script/variation/data/JBLabSeasonFlowering/data/K.tsv
		-P /Network/Data/250k/tmp-yh//phenotype/phenotype20100419.tsv -o /tmp/call_method_57_y8.tsv  -y8 -w 1 -n
	
Description:
	class to do association test on SNP data. option 'test_type' decides which test to run.
	
	Input genotype file format is Strain X SNP format (Yu's format, Output by DB_250k2data.py Or Output250KSNPs.py + ConvertBjarniSNPFormat2Yu.py).
		Each allele is either in atcgATCG... or 0(NA)1234(ACGT) digital format. It converts everything into integer matrix.
		If the input is in integer but has different meaning, like 0 is not NA, toggle noSNPAlleleOrdinalConversion.
	Input phenotype file format is Strain X phenotype format (Output by OutputPhenotype.py). 
	
	It requires a minimum number of ecotypes for either alleles of a single SNP to be eligible for kruskal wallis or linear model test.
	
	For all methods, it will automatically match strains in two files.
		NO worry for missing/extra data in either input file.
	
	All methods iterate through phenotypes given by '-w' except that Method "5: LM two phenotypes with PCs" takes two phenotypes from '-w'.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], \
				sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])
#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv, numpy, traceback
from pymodule import read_data, ProcessOptions, PassingData, SNPData, getListOutOfStr
from numpy import linalg

from Kruskal_Wallis import Kruskal_Wallis
from variation.src.db.output.OutputPhenotype import OutputPhenotype

if __name__ == '__main__':
	import rpy
from sets import Set
from pymodule import pca_module
#from DrawEcotypeOnMap import DrawEcotypeOnMap


class Association(Kruskal_Wallis):
	debug = 0
	report = 0
	__doc__ = __doc__
	option_default_dict = {}
	common_option_dict = Kruskal_Wallis.common_option_dict.copy()
	common_option_dict.update({
							('phenotype_method_id_ls', 0, ): ['1', 'w', 1, 'which phenotypes to work on. \
								a comma-dash-separated list phenotype_method ids in the phenotype file. \
								Check db Table phenotype_method. \
								if not available, take all phenotypes in the phenotype_fname.',],
							('eigen_vector_fname', 0, ): [None, 'f', 1, 'eigen vector file with PCs outputted by smartpca.perl from EIGENSOFT', ],\
							('kinship_fname', 0, ): [None, 'K', 1, 'file which contains the kinship matrix', ],\
							('genotype_fname_to_generate_kinship', 0, ): [None, 'G', 1, 'genotype file which is used to generate kinship, \
									if not given, kinship will be generated from inputFname', ],\
							('which_PC_index_ls', 0, ): [None, 'W', 1, 'list of indices indicating which PC(s) from eigen_vector_fname should be used. format: 0,1-3', ],\
							('noSNPAlleleOrdinalConversion', 0, ): [0, 'n', 0, 'by default (exept test-type 4), \
								it converts everything other than 0(=NA) into binary. toggle this for no such conversion.', ],\
							})
	#common_option_dict will be inherited by AssociationWorkflow.py
	option_default_dict.update(common_option_dict)
	option_default_dict.update({
						('output_fname', 1, ): ['', 'o', 1, 'file to store the pvalue. If multiple phenotypes are give, \
							phenotype id will be attached.', ],\
						('test_type', 1, int): [1, 'y', 1, 'Which type of test to do. \n\
	1: Kruskal_Wallis, 2:linear model(y=xb+e), \n\
	3: Emma, \n\
	4: LM with PCs, \n\
	5: LM two phenotypes with PCs, \n\
	6: LM two phenotypes with PCs, GeneXEnvironment Interaction, \n\
	7: Emma for genotype matrix without NA (no MAC and MAF output),\n\
	8: EMMAX (variance matrix is estimated once)'],\
						})
	def __init__(self, **keywords):
		"""
		2008-11-10
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		if self.phenotype_method_id_ls:
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
		
		self.which_PC_index_ls = getListOutOfStr(self.which_PC_index_ls, data_type=int)
		
		self.run_whole_matrix = {1:self._kruskal_wallis_whole_matrix,
								2:self.LM_whole_matrix,
								3:self.Emma_whole_matrix,
								4:self.LM_with_PCs_whole_matrix,
								5:self.LM_with_PCs_whole_matrix,
								6:self.LM_with_PCs_whole_matrix,
								7:self.Emma_whole_matrixForNoNAGenotypeMatrix,
								8:self.EMMAX}
		
		self.output_results = {1:self.output_kw_results,
							2:self.output_lm_results,
							3:self.output_lm_results,
							7:self.output_emma_results}
	
	@classmethod
	def getEigenValueFromFile(cls, eigen_value_fname):
		"""
		2008-12-03
		"""
		sys.stderr.write("Getting eigen values from %s ..."%eigen_value_fname)
		inf = open(eigen_value_fname)
		eigen_value_ls = []
		for line in inf:
			eigen_value = float(line.strip())
			eigen_value_ls.append(eigen_value)
		sys.stderr.write("Done.\n")
		return eigen_value_ls
	
	
	@classmethod
	def getPCFromFile(cls, eigen_vector_fname):
		"""
		2009-2-2
			add id_ls in the returning data
			wrap the data under PC_data
		2008-12-03
			
		"""
		sys.stderr.write("Getting principal components from %s ..."%eigen_vector_fname)
		inf = open(eigen_vector_fname)
		inf.next()	#skip first row (corresponding eigen values)
		PC_matrix = []
		id_ls = []
		for line in inf:
			row = line.split()
			PCs_for_one_entity = map(float, row[1:-1])	#1st entry in row is individual label. last entry is Case/Control.
			PC_matrix.append(PCs_for_one_entity)
			id_ls.append(row[0])
		PC_matrix = numpy.array(PC_matrix)
		PC_data = PassingData(PC_matrix=PC_matrix, id_ls=id_ls)
		sys.stderr.write("Done.\n")
		return PC_data
	
	@classmethod
	def multi_linear_model(cls, genotype_matrix, phenotype_ls, min_data_point=3, snp_index=None, kinship_matrix=None, eig_L=None, run_type=1):
		"""
		
		2008-11-13
			run_type
				1: pure_linear_model
				2: emma
			genotype_ls is already in binary (or integer starting from 0)
		2008-11-10
			similar to _kruskal_wallis() of Kruskal_Wallis
		"""
		non_NA_genotype_matrix = []
		non_NA_phenotype_ls = []
		non_NA_genotype2count = {}
		#non_NA_genotype2phenotype_ls = {}	#2008-08-06 try wilcox
		non_NA_index_ls = []
		for i in range(len(genotype_matrix)):
			if not numpy.isnan(phenotype_ls[i]):
				non_NA_index_ls.append(i)
				non_NA_phenotype = phenotype_ls[i]
				non_NA_phenotype_ls.append(non_NA_phenotype)
				
				non_NA_genotype_matrix.append(genotype_matrix[i,:])
				#non_NA_genotype2phenotype_ls[non_NA_genotype].append(phenotype_ls[i])	#2008-08-06 try wilcox
		"""
		#2008-08-06 try wilcox
		new_snp_allele2index = returnTop2Allele(non_NA_genotype2count)
		top_2_allele_ls = new_snp_allele2index.keys()
		non_NA_genotype2count = {top_2_allele_ls[0]: non_NA_genotype2count[top_2_allele_ls[0]],
								top_2_allele_ls[1]: non_NA_genotype2count[top_2_allele_ls[1]]}
		"""
		count_ls = non_NA_genotype2count.values()
		n = len(non_NA_phenotype_ls)
		if len(count_ls)>=2 and min(count_ls)>=min_data_point:	#require all alleles meet the min data point requirement
			if run_type==1:
				pdata = cls.pure_linear_model(non_NA_genotype_ls, non_NA_phenotype_ls)
			elif run_type==2:
				if kinship_matrix.shape[0]!=n:	#there is NA and need slicing
					new_kinship_matrix = kinship_matrix[non_NA_index_ls, non_NA_index_ls]
				else:
					new_kinship_matrix = kinship_matrix
				pdata = cls.emma(non_NA_genotype_ls, non_NA_phenotype_ls, new_kinship_matrix, eig_L)
			else:
				sys.stderr.write("run_type=%s not supported.\n"%run_type)
				return None
			pdata.snp_index = snp_index
			pdata.count_ls = count_ls
		else:
			pdata = None
		return pdata
	
	@classmethod
	def pure_linear_model(cls, non_NA_genotype_ls, non_NA_phenotype_ls, non_NA_phenotype2count=None,\
						lower_triangular_cholesky_inverse=None, beta_interval_percentage=0.95):
		"""
		2010-4-18
			now var_perc =  percentage of variance explained by the whole genotype matrix, rather than the first column
		2010-3-23
			do beta confidence interval estimation
			add argument beta_interval_percentage, default=0.95
			return S_square, d_jj_ls & beta_interval_delta_ls
		2010-2-28
			add argument lower_triangular_cholesky_inverse:
				run generalized least squares via ols (ordinary least squares after transforming the y and x
					using the lower-triangular matrix got from Cholesky decomposition of correlated variance_matrix).
				http://en.wikipedia.org/wiki/Generalized_least_squares or page 66 or Seber2003 LinearRegressionAnalysis
				
				1. non_NA_phenotype_ls is already transformed with the lower-triangular matrix but non_NA_genotype_ls is not yet.
				
		2009-12-18
			non_NA_phenotype2count is not used.
		2008-11-10
			split out of linear_model()
		"""
		genotype_matrix = cls.createDesignMatrix(non_NA_genotype_ls, add_intercept=True)	# need to add a constant vector as intercept.
		if lower_triangular_cholesky_inverse is not None:	# 2010-2-28
			genotype_matrix = numpy.dot(lower_triangular_cholesky_inverse, genotype_matrix)
		
		no_of_rows, no_of_cols = genotype_matrix.shape
		
		genotype_matrix2 = numpy.transpose(genotype_matrix)
		D = linalg.inv(numpy.inner(genotype_matrix2, genotype_matrix2))	#2nd matrix's shape is opposite to real matrix algebra.

		(p, residuals, rank, s) = linalg.lstsq(genotype_matrix, non_NA_phenotype_ls)	#all 4 returned elements are arrays
		# p is a list of coeffcients.
		#coeff_list = [p[0]]	#put the intercept beta there first
		#coeff_p_value_list = [-1]
		coeff_list = []
		coeff_p_value_list = []
		d_jj_ls = []	# 2010-3-23
		degrees_of_freedom = no_of_rows-len(p)
		S_square = residuals[0]/degrees_of_freedom	# 2010-3-23 estimate for the RSS page 106 of Seber2003LRA
		beta_interval_delta_ls = []	# 2010-3-23
		# 2010-3-23 the t-stat corresponding to the two-sided interval percentage
		interval_t_stat = rpy.r.qt((1-beta_interval_percentage)/2, degrees_of_freedom, \
								lower_tail=rpy.r.FALSE)
		
		var_perc_ls = []	# 2010-3-26
		phenotype_variance = numpy.var(non_NA_phenotype_ls)
		for i in range(len(p)):
			coeff_list.append(p[i])
			F = p[i]*p[i]/(S_square*D[i,i])	#page 106 of Seber2003LRA
				#page 108 of Seber2003LRA. F=r^2(n-2)/(1-r^2).
				#T=sqrt(F) is t-stat with n-2 df. it's used in r.glm(). F-test gives slightly bigger p-value than t-test.
			pvalue = rpy.r.pf(F,1, degrees_of_freedom, lower_tail=rpy.r.FALSE)
			coeff_p_value_list.append(pvalue)
			d_jj_ls.append(D[i,i])	# 2010-3-23 page 106 of Seber2003LRA
			# 2010-3-23 get beta interval delta
			beta_interval_delta = interval_t_stat*math.sqrt(S_square*D[i,i])
			beta_interval_delta_ls.append(beta_interval_delta)
			
			# 2010-3-26 calculate variance explained for every covariate
			genotype_var = numpy.var(genotype_matrix[:,i])
			geno_effect_var = genotype_var*p[i]*p[i]
			var_perc = geno_effect_var/phenotype_variance
			var_perc_ls.append(var_perc)
		
		pvalue = coeff_p_value_list[1]	#this is the pvalue to return, corresponding to the genotype vector
		# 2008-11-10 variance of only the first covariate (after intercept)
		#genotype_var = numpy.var(genotype_matrix[:,1]) 	#2008-11-10 var=\sum(x_i-\bar{x})^2/(n-1)
		#geno_effect_var = genotype_var*p[1]*p[1]*(no_of_rows-1)
		#var_perc = geno_effect_var/(residuals[0]+geno_effect_var)	# 2010-3-27 only good when there's only one non-intercept covariate.
		#geno_effect_var = genotype_var*p[1]*p[1]	# varia
		
		# 2010-4-18	variance of the whole genotype matrix.
		geno_effect_var = numpy.var(numpy.dot(genotype_matrix, coeff_list))
		var_perc = geno_effect_var/phenotype_variance
		
		
		#pvalue = rpy.r.kruskal_test(x=non_NA_phenotype_ls, g=rpy.r.as_factor(non_NA_genotype_ls))['p.value']
		#2008-08-06 try wilcox
		#pvalue = rpy.r.wilcox_test(non_NA_genotype2phenotype_ls[top_2_allele_ls[0]], non_NA_genotype2phenotype_ls[top_2_allele_ls[1]], conf_int=rpy.r.TRUE)['p.value']
		pdata = PassingData(pvalue=pvalue, var_perc=var_perc, coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list,
						S_square = S_square, d_jj_ls = d_jj_ls, beta_interval_delta_ls=beta_interval_delta_ls,\
						var_perc_ls=var_perc_ls)
		return pdata
	
	@classmethod
	def pure_linear_model_via_R(cls, non_NA_genotype_ls, non_NA_phenotype_ls, non_NA_phenotype2count=None):
		"""
		2010-2-25
			use createDesignMatrix() to generate a design matrix
		2009-8-28
			split out of pure_linear_model(). same functionality as pure_linear_model(), but invoke R to run regression.
		"""
		
		genotype_matrix = cls.createDesignMatrix(non_NA_genotype_ls)  # no need to add a constant vector.
		no_of_rows, no_of_cols = genotype_matrix.shape
		#2008-11-10 do linear regression by R
		genotype_var = numpy.var(genotype_matrix[:,0]) 	#2008-11-10 var=\sum(x_i-\bar{x})^2/(n-1)
		rpy.set_default_mode(rpy.NO_CONVERSION) #04-07-05
		#data_frame = rpy.r.as_data_frame({"phenotype":non_NA_phenotype_ls, "genotype":rpy.r.as_factor(genotype_matrix[:,1])})
		formula_list = []
		data_frame_dict = {"phenotype":non_NA_phenotype_ls}
		for i in range(genotype_matrix.shape[1]):
			var_name = 'genotype%s'%i
			formula_list.append(var_name)
			data_frame_dict.update({var_name: genotype_matrix[:,i]})
		data_frame = rpy.r.as_data_frame(data_frame_dict)
		formula = 'phenotype~%s'%'+'.join(formula_list)
		
		if non_NA_phenotype2count and len(non_NA_phenotype2count)==2:	#binary phenotype, use logistic regression
			lm_result = rpy.r.glm(rpy.r(formula), data=data_frame, family=rpy.r("binomial"))
		else:
			lm_result = rpy.r.glm(rpy.r(formula), data=data_frame)
		rpy.set_default_mode(rpy.BASIC_CONVERSION)
		#04-07-05 r.summary() requires lm_result in NO_CONVERSION state
		summary_stat = rpy.r.summary(lm_result)
		
		#06-30-05	index 0 in summary_stat['coefficients'] is intercept
		coeff_list = []
		coeff_p_value_list = []
		for i in range(len(summary_stat['coefficients'])):
			coeff_list.append(summary_stat['coefficients'][i][0])	#0 is the coefficient
			coeff_p_value_list.append(summary_stat['coefficients'][i][-1])	#-1 is the corresponding p-value
		#06-30-05	fill in other efficients based on bit_string, NOTE i+1
		pvalue = coeff_p_value_list[1]
		residuals = summary_stat['deviance']
		geno_effect_var = genotype_var*coeff_list[1]*coeff_list[1]*(no_of_rows-1)
		var_perc = geno_effect_var/(residuals+geno_effect_var)
		
		pdata = PassingData(pvalue=pvalue, var_perc=var_perc, coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list)
		return pdata
	
	@classmethod
	def createDesignMatrix(cls, genotype_ls, add_intercept=False,n=None):
		"""
		2009-12-23
			argument add_intercept determines whether a constant vector would be added to the design matrix or not
		"""
		if genotype_ls is None:
			if add_intercept:
				if n is None:
					sys.stderr.write("n (number of accessions) must be specified to create a design matrix with genotype_ls=None.\n")
					sys.exit(3)
				design_matrix = (numpy.ones([n ,1], numpy.int))
		else:
			#rpy.r.as_factor.local_mode(rpy.NO_CONVERSION)
			if not isinstance(genotype_ls, numpy.ndarray):
				design_matrix = numpy.array(genotype_ls)
			else:
				design_matrix = genotype_ls
			if len(design_matrix.shape)==1:	#transform into 2D array
				design_matrix = numpy.resize(design_matrix, [len(design_matrix), 1])				
			no_of_rows, no_of_cols = design_matrix.shape
			
			
			if add_intercept:
				#2008-11-10 do linear regression by numpy
				design_matrix = numpy.hstack((numpy.ones([len(design_matrix),1], numpy.int), design_matrix))	#the design variable for intercept has to be included.
			
		return design_matrix
	
	@classmethod
	def gls_with_transformed_y(cls, non_NA_genotype_ls, non_NA_phenotype_ls, non_NA_phenotype2count=None, \
							lower_triangular_cholesky_inverse=None):
		"""
		2010-2-28
			run generalized least squares via ols (ordinary least squares after transforming the y and x
				using the lower-triangular matrix got from Cholesky decomposition of correlated variance_matrix).
			http://en.wikipedia.org/wiki/Generalized_least_squares or page 66 or Seber2003 LinearRegressionAnalysis
			
			1. non_NA_phenotype_ls is already transformed with the lower-triangular matrix but non_NA_genotype_ls is not yet.
		"""
		return cls.pure_linear_model(non_NA_genotype_ls, non_NA_phenotype_ls, \
									lower_triangular_cholesky_inverse=lower_triangular_cholesky_inverse)
	
	@classmethod
	def gls_via_R(cls, non_NA_genotype_ls, non_NA_phenotype_ls, non_NA_phenotype2count=None, variance_matrix=None):
		"""
		2010-3-1
			deprecated. gls_with_transformed_y() replaces this.
		2009-12-23
			generalized least squares (not general least square) model via calling equivalent function in R.
			
		"""
		genotype_matrix = cls.createDesignMatrix(non_NA_genotype_ls)  # no need to add a constant vector.
		
		no_of_rows = len(non_NA_genotype_ls)
		
		#2008-11-10 do linear regression by R
		genotype_var = numpy.var(genotype_matrix[:,0]) 	#2008-11-10 var=\sum(x_i-\bar{x})^2/(n-1)
		rpy.set_default_mode(rpy.NO_CONVERSION) #04-07-05
		rpy.r.library("nlme")
		
		#data_frame = rpy.r.as_data_frame({"phenotype":non_NA_phenotype_ls, "genotype":rpy.r.as_factor(genotype_matrix[:,1])})
		formula_list = []
		data_frame_dict = {"phenotype":non_NA_phenotype_ls}
		for i in range(genotype_matrix.shape[1]):
			var_name = 'genotype%s'%i
			formula_list.append(var_name)
			data_frame_dict.update({var_name: genotype_matrix[:,i]})
		data_frame = rpy.r.as_data_frame(data_frame_dict)
		formula = 'phenotype~%s'%'+'.join(formula_list)
		
		if hasattr(cls, 'corStruct'):
			corStruct = cls.corStruct
		else:
			if variance_matrix is not None:
				corStruct = cls.generateCorStructForGLSFromVarianceMatrix(variance_matrix)
				setattr(cls, "corStruct", corStruct)
			else:
				corStruct = None
		
		lm_result = rpy.r.gls(rpy.r(formula), data=data_frame, correlation=corStruct)
		#2010-2-28 this error occurred while invoking gls(): function evaluation limit reached without convergence (9)
		
		#04-07-05 r.summary() requires lm_result in NO_CONVERSION state
		rpy.set_default_mode(rpy.BASIC_CONVERSION)
		summary_stat = rpy.r.summary(lm_result)
		
		#06-30-05	index 0 in summary_stat['coefficients'] is intercept
		coeff_list = []
		coeff_p_value_list = []
		# 2010-2-28 this summary_stat['coefficients'] has a different structure from the one got from glm().
		#	It also doesn't have pvalues attached. code below doesn't work in reality.
		for i in range(len(summary_stat['coefficients'])):
			coeff_list.append(summary_stat['coefficients'][i][0])	#0 is the coefficient
			coeff_p_value_list.append(summary_stat['coefficients'][i][-1])	#-1 is the corresponding p-value
		#06-30-05	fill in other efficients based on bit_string, NOTE i+1
		pvalue = coeff_p_value_list[1]
		residuals = summary_stat['deviance']
		geno_effect_var = genotype_var*coeff_list[1]*coeff_list[1]*(no_of_rows-1)
		var_perc = geno_effect_var/(residuals+geno_effect_var)
		
		pdata = PassingData(pvalue=pvalue, var_perc=var_perc, coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list)
		return pdata
	
	
	@classmethod
	def linear_model(cls, genotype_ls, phenotype_ls, min_data_point=3, snp_index=None, kinship_matrix=None, eig_L=None, \
					run_type=1, counting_and_NA_checking=True, variance_matrix=None, \
					lower_triangular_cholesky_inverse=None):
		"""
		2010-2-28
			add argument lower_triangular_cholesky_inverse for run_type 4 gls_with_transformed_y().
		2009-12-23
			add argument variance_matrix for gls()
		2009-12-18
			add argument counting_and_NA_checking, which if enabled, this function counts the number of accessions with either allele
				(, which allows calculation of MAF and makes sure that MAF wasn't too small)  
				and exclude the NA genotypes or phenotypes.
			If it's disabled, pure_linear_model_via_R() would always call glm() with no family specification. Others are same.
		2009-2-28
			run_type
				1: pure_linear_model
				2: emma
				3: pure_linear_model via R
				4: generalized least square with specified variance matrix
		2009-5-15
			report full error information if cls.debug is set.
		2009-2-11
			simply report the snp_index when numpy.linalg.linalg.LinAlgError is caught
		2008-12-04
			genotype_ls could be either 1D or 2D matrix
			if it's 2D, only the 1st column is genotype. 2nd or further are PCs (so far).
		2008-11-13
			run_type
				1: pure_linear_model
				2: emma
			genotype_ls is already in binary (or integer starting from 0)
		2008-11-10
			similar to _kruskal_wallis() of Kruskal_Wallis
		"""
		if counting_and_NA_checking:
			non_NA_genotype_ls = []
			non_NA_phenotype_ls = []
			non_NA_genotype2count = {}
			#non_NA_genotype2allele = {}
			non_NA_phenotype2count = {}
			#non_NA_genotype2phenotype_ls = {}	#2008-08-06 try wilcox
			non_NA_index_ls = []
			for i in range(len(genotype_ls)):
				if isinstance(genotype_ls[i], numpy.ndarray):
					genotype = genotype_ls[i][0]
				else:
					genotype = genotype_ls[i]
				if genotype!=-2 and not numpy.isnan(genotype) and not numpy.isnan(phenotype_ls[i]):
					non_NA_index_ls.append(i)
					non_NA_genotype = genotype
					non_NA_phenotype = phenotype_ls[i]
					if non_NA_phenotype not in non_NA_phenotype2count:
						non_NA_phenotype2count[non_NA_phenotype] = 0
					non_NA_phenotype2count[non_NA_phenotype] += 1
					non_NA_phenotype_ls.append(non_NA_phenotype)
					
					if non_NA_genotype not in non_NA_genotype2count:
						non_NA_genotype2count[non_NA_genotype] = 0
						#non_NA_genotype2allele[non_NA_genotype] = len(non_NA_genotype2allele)
						#non_NA_genotype2phenotype_ls[non_NA_genotype] = []	#2008-08-06 try wilcox
					#allele = non_NA_genotype2allele[non_NA_genotype]
					allele = genotype_ls[i]
					if not isinstance(allele, numpy.ndarray):	#make non_NA_genotype_ls 2D
						allele = [allele]
					non_NA_genotype_ls.append(allele)
					non_NA_genotype2count[non_NA_genotype] += 1
					#non_NA_genotype2phenotype_ls[non_NA_genotype].append(phenotype_ls[i])	#2008-08-06 try wilcox
		else:
			non_NA_genotype_ls = genotype_ls
			non_NA_phenotype_ls = phenotype_ls
			non_NA_genotype2count = {0:min_data_point, 1:min_data_point}	# fake one and allow it to pass the condition below
			non_NA_phenotype2count = {}	# based on its length, pure_linear_model_via_R would decide whether to call glm() with 'binomial'
										# if phenotypes are binary or not. So if this is empty, envoke glm() with no family specification.
		"""
		#2008-08-06 try wilcox
		new_snp_allele2index = returnTop2Allele(non_NA_genotype2count)
		top_2_allele_ls = new_snp_allele2index.keys()
		non_NA_genotype2count = {top_2_allele_ls[0]: non_NA_genotype2count[top_2_allele_ls[0]],
								top_2_allele_ls[1]: non_NA_genotype2count[top_2_allele_ls[1]]}
		"""
		count_ls = non_NA_genotype2count.values()
		n = len(non_NA_phenotype_ls)
		if len(count_ls)>=2 and min(count_ls)>=min_data_point:	#require all alleles meet the min data point requirement
			try:
				if run_type==1:
					pdata = cls.pure_linear_model(non_NA_genotype_ls, non_NA_phenotype_ls)
				elif run_type==3:
					pdata = cls.pure_linear_model_via_R(non_NA_genotype_ls, non_NA_phenotype_ls, non_NA_phenotype2count)
				elif run_type==2:
					if kinship_matrix.shape[0]!=n:	#there is NA and need slicing
						new_kinship_matrix = kinship_matrix[non_NA_index_ls, non_NA_index_ls]
						eig_L = rpy.r.emma_eigen_L(None, new_kinship_matrix)	#2008-1-5 generate new eig_L
					else:
						new_kinship_matrix = kinship_matrix
					pdata = cls.emma(non_NA_genotype_ls, non_NA_phenotype_ls, new_kinship_matrix, eig_L)
				elif run_type==4:
					pdata = cls.gls_with_transformed_y(non_NA_genotype_ls, non_NA_phenotype_ls, \
								lower_triangular_cholesky_inverse=lower_triangular_cholesky_inverse)
					#pdata = cls.gls_via_R(non_NA_genotype_ls, non_NA_phenotype_ls, variance_matrix=variance_matrix)
				else:
					sys.stderr.write("run_type=%s not supported.\n"%run_type)
					return None
				pdata.snp_index = snp_index
				pdata.count_ls = count_ls
			
			except numpy.linalg.linalg.LinAlgError:
				sys.stderr.write("Except while running pure_linear_model on snp_index=%s:\n"%(snp_index))
				sys.stderr.write('\t%s.\n'%repr(sys.exc_info()))
				pdata = None
			except:
				if cls.debug:	# 2009-5-15
					sys.stderr.write("Except while running pure_linear_model on snp_index=%s with non_NA_genotype_ls=%s, \
					non_NA_phenotype_ls=%s, non_NA_phenotype2count=%s.\n"%(snp_index, repr(non_NA_genotype_ls), \
																			repr(non_NA_phenotype_ls), repr(non_NA_phenotype2count)))
					traceback.print_exc()
					sys.stderr.write('%s.\n'%repr(sys.exc_info()))
				else:	#2009-4-5 simpilify output if not debug
					sys.stderr.write("Except (%s) while running pure_linear_model on snp_index=%s.\n"%(repr(sys.exc_info()), snp_index))
				pdata = None
		else:
			pdata = None
		return pdata
	
	
	def LM_whole_matrix(self, snpData, phenotype_ls, min_data_point=3, **keywords):
		"""
		2008-11-10
			adapted from _kruskal_wallis_whole_matrix() of Kruskal_Wallis
		"""
		sys.stderr.write("Association by pure linear model ...\n")
		no_of_rows, no_of_cols = snpData.data_matrix.shape
		results = []
		counter = 0
		real_counter = 0
		for j in range(no_of_cols):
			genotype_ls = snpData.data_matrix[:,j]
			pdata = self.linear_model(genotype_ls, phenotype_ls, min_data_point, snp_index=j)
			if pdata is not None:
				results.append(pdata)
				real_counter += 1
			counter += 1
			if self.report and counter%2000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		if self.report:
			sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		sys.stderr.write("Done.\n")
		return results
	
	def LM_with_PCs_whole_matrix(self, snpData, phenotype_ls, min_data_point=3, **keywords):
		"""
		2008-01-04
			add code to deal with environment_matrix
		2008-12-04
		"""
		sys.stderr.write("Association by linear model with principle components ...\n")
		no_of_rows, no_of_cols = snpData.data_matrix.shape
		which_PC_index_ls = keywords.get('which_PC_index_ls') or [0]
		PC_matrix = keywords.get('PC_matrix')
		environment_matrix = keywords.get('environment_matrix')
		gene_environ_interaction = keywords.get('gene_environ_interaction')
		
		results = []
		counter = 0
		real_counter = 0
		if PC_matrix is not None:
			sub_PC_matrix = PC_matrix[:, which_PC_index_ls]
		else:
			sub_PC_matrix = None
		for j in range(no_of_cols):
			genotype_ls = snpData.data_matrix[:,j]
			genotype_ls = numpy.resize(genotype_ls, [len(genotype_ls),1])	#make it 2D , so able to hstack with sub_PC_matrix
			
			if gene_environ_interaction:
				gene_environ_matrix = genotype_ls*environment_matrix
			
			if environment_matrix is not None:
				genotype_ls = numpy.hstack((genotype_ls, environment_matrix))
			
			if gene_environ_interaction:
				genotype_ls = numpy.hstack((genotype_ls, gene_environ_matrix))
			if sub_PC_matrix is not None:
				genotype_ls = numpy.hstack((genotype_ls, sub_PC_matrix))
			
			pdata = self.linear_model(genotype_ls, phenotype_ls, min_data_point, snp_index=j)
			if pdata is not None:
				results.append(pdata)
				real_counter += 1
			counter += 1
			if self.report and counter%2000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		if self.report:
			sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		sys.stderr.write("Done.\n")
		return results
	
	@classmethod
	def removeRowsWithNAPhenotypeFromKinshipAndSNPData(cls, phenotype_ls, snpData, kinshipData=None):
		"""
		2011-4-27
			a preprocessing step for EMMAX(), Emma_whole_matrixForNoNAGenotypeMatrix(), Emma_whole_matrix()
			1. generate a list of phenotypes in which each value is non-NA.
			2. generate a new SNPData in which only non-NA-phenotype rows of snpData are retained
			3. generate a new kinshipData so that each row corresponds to each row in new SNPData. 
		"""
		#remove NA phenotype and genotypes from the corresponding accessions
		non_NA_phenotype_row_index_ls = []
		non_NA_phenotype_ls = []
		for i in range(len(phenotype_ls)):
			if not numpy.isnan(phenotype_ls[i]):
				non_NA_phenotype_row_index_ls.append(i)
				non_NA_phenotype_ls.append(phenotype_ls[i])
		
		newSNPData = SNPData.keepRowsByRowIndex(snpData, non_NA_phenotype_row_index_ls)
		new_data_matrix = newSNPData.data_matrix
		if not kinshipData:	#if it's None, generate it.
			kinshipData = SNPData(row_id_ls=newSNPData.row_id_ls, col_id_ls=newSNPData.row_id_ls, data_matrix=newSNPData.get_kinship_matrix())
		kinshipData = cls.adjustKinshipBySNPData(kinshipData, newSNPData)[0]
		return PassingData(kinshipData=kinshipData, snpData=newSNPData, non_NA_phenotype_row_index_ls=non_NA_phenotype_row_index_ls,\
						non_NA_phenotype_ls=non_NA_phenotype_ls)
	
	@classmethod
	def adjustKinshipBySNPData(cls, kinshipData, snpData):
		"""
		2011-4-27
			In case that kinshipData.row_id_ls and snpData.row_id_ls are not same, 
				generate a kinshipData that has the same row_id_ls as snpData.
		"""
		if snpData.row_id_ls!=kinshipData.row_id_ls:
			Z = cls.createIndividualToLineIncidenceMatrix(snpData.row_id_ls, kinshipData.row_id_ls)
			new_K = numpy.dot(numpy.dot(Z, kinshipData.data_matrix), numpy.transpose(Z))
			kinshipData = SNPData(row_id_ls=snpData.row_id_ls, col_id_ls=snpData.row_id_ls, \
								data_matrix=new_K)
		else:
			Z = None
		return kinshipData,Z
	
	@classmethod
	def getExpandedKinship(cls, kinship_fname=None, genotype_fname_to_generate_kinship=None, independentSNPData=None):
		"""
		2011-4-25
			moved from class JBDataGWA of variation/src/misc.py
		2010-8-22
			return independentSNPData as well
		2010-4-23
			split out of getCholeskyInverseData()
		"""
		if kinship_fname and os.path.isfile(kinship_fname):
			kinshipData = SNPData(input_fname=kinship_fname, ignore_2nd_column=1, \
								matrix_data_type=float, turn_into_array=1)
		else:
			if genotype_fname_to_generate_kinship is not None:
				snpData = SNPData(input_fname=genotype_fname_to_generate_kinship, turn_into_array=1, ignore_2nd_column=1)
				newSnpData, allele_index2allele_ls = snpData.convert2Binary()
				snpData = newSnpData
			else:
				snpData = independentSNPData
			kinship_matrix = snpData.get_kinship_matrix()
			kinshipData = SNPData(row_id_ls=snpData.row_id_ls, col_id_ls=snpData.row_id_ls, data_matrix=kinship_matrix)
			del snpData
			if kinship_fname is not None:
				kinshipData.tofile(kinship_fname)
		
		#2010-8-22 a new independentSNPData is generated
		independentSNPData = independentSNPData.removeRowsNotInTargetSNPData(kinshipData)
		kinshipData, Z = cls.adjustKinshipBySNPData(kinshipData, independentSNPData)
		return kinshipData, independentSNPData, Z	#2010-8-22 return independentSNPData as well
	
	@classmethod
	def createIndividualToLineIncidenceMatrix(cls, individualID_ls, lineID_ls):
		"""
		2011-4-25
			moved from class JBDataGWA of variation/src/misc.py
		2010-3-22
			multiple individuals in individualID_ls share the same ID (lineID).
			The purpose of this function is to create a matrix mapping each individual to a line, 
				which would be used to create a new kinship matrix.
		"""
		sys.stderr.write("Creating individual-to-line incidence matrix ...")
		no_of_replicates = len(individualID_ls)
		no_of_lines = len(lineID_ls)
		Z = numpy.zeros([no_of_replicates, no_of_lines], numpy.int8)	# incidence matrix
		no_of_replicates_with_line = 0
		for i in range(no_of_replicates):
			individualID = individualID_ls[i]
			for j in range(no_of_lines):
				if individualID==lineID_ls[j]:
					Z[i,j] = 1
					no_of_replicates_with_line += 1
					break
		sys.stderr.write("%s out of %s in list 1 found corresponding entry in the 2nd %s-entry list.\n"%\
						(no_of_replicates_with_line, no_of_replicates, len(lineID_ls)))
		return Z
	
	
	emma_R_sourced = False	# 2010-4-23
	@classmethod
	def emma(cls, non_NA_genotype_ls=None, non_NA_phenotype_ls=None, kinship_matrix=None, eig_L = None, Z=None,):
		"""
		2010-9-8
			add argument Z, an incidence matrix mapping accessions from non_NA_phenotype_ls to kinship_matrix
			return LL, dLL, optLL as well
		2010-9-6
			return delta
		2010-9-6
			return heritability
		2010-4-23
			source the emma.R first if it's not sourced yet.
		2009-12-19
			non_NA_genotype_ls could be None/empty list, which means it won't be used.
		2009-8-26
			update coeff_p_value_list with pvalues returned from EMMA
		2008-11-13
			call emma.REMLE()
		"""
		if not cls.emma_R_sourced:	#2010-4-23 not sourced yet.
			rpy.r.source(os.path.expanduser('~/script/variation/src/gwa/emma/R/emma.R'))
			cls.emma_R_sourced = True
		len_functor = getattr(non_NA_genotype_ls, '__len__', None)
		if len_functor:
			n = len_functor()
			if n==0:	# if non_NA_genotype_ls is an empty list, same effect as None.
				non_NA_genotype_ls = None	#2010-8-22 empty list could mislead createDesignMatrix().
				n = None
		else:
			n = None
		
		if n is None:
			if non_NA_phenotype_ls is not None:
				n = len(non_NA_phenotype_ls)
			elif kinship_matrix is not None:
				n = kinship_matrix.shape[0]
			else:
				n = None	# this would cause error
		genotype_matrix = cls.createDesignMatrix(non_NA_genotype_ls, add_intercept=True, n=n)
		one_marker_rs = rpy.r.emma_REMLE(non_NA_phenotype_ls, genotype_matrix, kinship_matrix, Z=Z, eig_L=eig_L, cal_pvalue=rpy.r.TRUE)
		coeff_list = [beta[0] for beta in one_marker_rs['beta']]
		if non_NA_genotype_ls is None:    # no covariate used, so no pvalues
			coeff_p_value_list=[0.]
		else:
			coeff_p_value_list = [pvalue for pvalue in one_marker_rs['pvalues']]
		
		"""# 2009-3-24 temporary, get random_effect+residual
		coeff_ar = numpy.array(coeff_list)
		random_effect_and_residual_ls = list(numpy.array(non_NA_phenotype_ls)-numpy.inner(genotype_matrix, coeff_ar))
		coeff_list.extend(random_effect_and_residual_ls)
		"""
		# vg is the variance for the random effect (multi-small-effect-gene effect), ve is the residual variance. 
		ve=one_marker_rs['ve']
		vg=one_marker_rs['vg']
		heritability = vg/(vg+ve)
		pdata = PassingData(pvalue=one_marker_rs['pvalue'], var_perc=one_marker_rs['genotype_var_perc'][0][0], \
						coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list, ve=ve, \
						vg=vg, heritability=heritability, delta=one_marker_rs['delta'], logdelta=one_marker_rs['logdelta'],\
						optLL=one_marker_rs['optLL'], dLL=one_marker_rs['dLL'], LL=one_marker_rs['LL'])
		return pdata
	
	def Emma_whole_matrix(self, snpData, phenotype_ls, min_data_point=3, **keywords):
		"""
		2011-4-27
			argument data_matrix renamed to snpData.
			preprocessing steps are merged into self.removeRowsWithNAPhenotypeFromKinshipAndSNPData().
			If keywords contains kinshipData, row_id and col_id should be same as snpData's row_id.
			If keywords doesn't contain kinshipData, kinshipData will be generated based on snpData.
		2008-11-13
			iteratively call rpy.r.emma_REMLE() through self.linear_model() in order to get MAF, MAC etc
			stop calling rpy.r.emma_REML_t
		2008-11-11
			assume:
				1. no NA in data_matrix (imputed) but allow NA in phenotype_ls
				2. data_matrix is binary, used in get_kinship_matrix()
			procedure:
			
			1. remove rows that have NA in phenotype_ls
			2. calculate kinship
			3. run emma
		"""
		sys.stderr.write("Association by emma (linear mixture model) ...\n")
		returnData = self.removeRowsWithNAPhenotypeFromKinshipAndSNPData(phenotype_ls, snpData, keywords.get('kinshipData'))
		new_data_matrix = returnData.snpData.data_matrix
		kinship_matrix = returnData.kinshipData.data_matrix
		
		no_of_rows, no_of_cols = new_data_matrix.shape
		results = []
		counter = 0
		real_counter = 0
		rpy.r.source(os.path.expanduser('~/script/variation/src/gwa/emma/R/emma.R'))
		
		eig_L = rpy.r.emma_eigen_L(None, kinship_matrix)	#to avoid repeating the computation of eig_L inside emma.REMLE
		for j in range(no_of_cols):
			genotype_ls = new_data_matrix[:,j]
			pdata = self.linear_model(genotype_ls, returnData.non_NA_phenotype_ls, min_data_point, snp_index=j, \
									kinship_matrix=kinship_matrix, eig_L=eig_L, run_type=2)
			if pdata is not None:
				results.append(pdata)
				real_counter += 1
			counter += 1
			if self.report and counter%2000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		if self.report:
			sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		
		sys.stderr.write("Done.\n")
		return results
	
	def Emma_whole_matrixForNoNAGenotypeMatrix(self, snpData, phenotype_ls, min_data_point=3, **keywords):
		"""
		2011-4-27
			argument data_matrix renamed to snpData.
			preprocessing steps are merged into self.removeRowsWithNAPhenotypeFromKinshipAndSNPData().
			If keywords contains kinshipData, row_id and col_id should be same as snpData's row_id.
			If keywords doesn't contain kinshipData, kinshipData will be generated based on snpData.
		2009-3-18
			carved out of old Emma_whole_matrix() for data_matrix that has no NA in it.
			no MAC and MAF information.
		"""
		sys.stderr.write("Association by monolithic-emma (linear mixture model) ...\n")
		
		returnData = self.removeRowsWithNAPhenotypeFromKinshipAndSNPData(phenotype_ls, snpData, keywords.get('kinshipData'))
		new_data_matrix = returnData.snpData.data_matrix
		kinship_matrix = returnData.kinshipData.data_matrix
		
		non_NA_phenotype_ar = numpy.array(returnData.non_NA_phenotype_ls)
		
		if not cls.emma_R_sourced:	#2010-4-23 not sourced yet.
			rpy.r.source(os.path.expanduser('~/script/variation/src/gwa/emma/R/emma.R'))
			cls.emma_R_sourced = True
		#rpy.set_default_mode(rpy.NO_CONVERSION)
		non_NA_phenotype_ar.resize([len(non_NA_phenotype_ar),1])	#transform it into 2-D for emma
		results = rpy.r.emma_REML_t(non_NA_phenotype_ar.transpose(), new_data_matrix.transpose(), kinship_matrix)	#in emma, row is marker. col is strain.
		
		sys.stderr.write("Done.\n")
		return results
	
	@classmethod
	def generateCorStructForGLSFromVarianceMatrix(cls, variance_matrix, data_frame=None):
		"""
		2010-3-1
			deprecated. since gls_via_R() is not used anymore
		2010-2-26
			set the mode back to the original mode
		2009-12-23
			generate the corStruct for gls_via_R()
		"""
		sys.stderr.write("Generating corStruct for gls() from variance_matrix ...")
		
		rpy_default_mode = rpy.get_default_mode()
		rpy.set_default_mode(rpy.NO_CONVERSION) #04-07-05
		rpy.r.library("nlme")
		
		# bring the lower-triangle of variance_matrix into a list, row by row
		no_of_rows, no_of_cols = variance_matrix.shape
		lower_triangle_cor_vector = []
		for i in range(0, no_of_rows-1):
			for j in range(i+1, no_of_rows):
				lower_triangle_cor_vector.append(variance_matrix[i][j]/math.sqrt(variance_matrix[i][i]*variance_matrix[j][j]))
		
		csSymm = rpy.r.corSymm(value=lower_triangle_cor_vector)
		if data_frame is None:
			data_frame = rpy.r.as_data_frame({"fakedata":[1]*no_of_rows})
		csSymm = rpy.r.Initialize(csSymm, data=data_frame)
		rpy.set_default_mode(rpy_default_mode)	# 2010-2-26 set the mode back to the original mode
		sys.stderr.write("Done.\n")
		return csSymm
	
	@classmethod
	def preEMMAX(cls, non_NA_phenotype_ls, kinship_matrix, non_NA_genotype_ls=None, debug=False, Z=None):
		"""
		2010-9-8
			add argument Z, an incidence matrix mapping accessions from non_NA_phenotype_ls to kinship_matrix
			return one_emma_rs from cls.emma()
		2010-4-16
			add argument non_NA_genotype_ls, to allow co-factors
		2010-3-16
			pre-EMMAX procedures:
				1. run EMMA once to get variance matrix,
				2. cholesky decomposition
				3. transform the phenotype
		"""
		sys.stderr.write("Running pre-EMMAX procedures ...")
		if not cls.emma_R_sourced:	#2010-4-23 not sourced yet.
			rpy.r.source(os.path.expanduser('~/script/variation/src/gwa/emma/R/emma.R'))
			cls.emma_R_sourced = True
		# run EMMA here to get vg & ve, scalars for the two variance matrices (random effect, residual) 
		one_emma_rs = cls.emma(non_NA_genotype_ls=non_NA_genotype_ls, \
			non_NA_phenotype_ls=non_NA_phenotype_ls, kinship_matrix=kinship_matrix, eig_L=None, Z=Z)
		vg = one_emma_rs.vg
		ve = one_emma_rs.ve
		
		sys.stderr.write("vg=%s, ve=%s.\n"%(vg, ve))
		no_of_rows, no_of_cols = kinship_matrix.shape
		
		variance_matrix = vg*kinship_matrix + ve*numpy.identity(no_of_rows)
		
		# 2010-2-28
		if debug:
			sys.stderr.write("\t cholesky decomposition on variance_matrix ...")
		L = numpy.linalg.cholesky(variance_matrix)
		if debug:
			sys.stderr.write("Done.\n")
		if debug:
			sys.stderr.write("\t getting inverse ...")
		L_inverse = numpy.linalg.inv(L)
		if debug:
			sys.stderr.write("Done.\n")
		non_NA_phenotype_ar = numpy.array(non_NA_phenotype_ls)		
		non_NA_phenotype_ar = numpy.dot(L_inverse, non_NA_phenotype_ar)	# numpy.dot and numpy.inner has subtle difference.
		sys.stderr.write("Done.\n")
		return PassingData(variance_matrix=variance_matrix, non_NA_phenotype_ar=non_NA_phenotype_ar, L_inverse=L_inverse, 
					one_emma_rs = one_emma_rs)
		
	def EMMAX(self, snpData, phenotype_ls, min_data_point=3, **keywords):
		"""
		2011-4-27
			argument data_matrix renamed to snpData.
			preprocessing steps are merged into self.removeRowsWithNAPhenotypeFromKinshipAndSNPData().
			If keywords contains kinshipData, row_id and col_id should be same as snpData's row_id.
			If keywords doesn't contain kinshipData, kinshipData will be generated based on snpData.
		2010-4-23
			the bug fixed on 2010-4-22 (below) is not fixed.
				no_of_cols shouldn't be derived from kinship_matrix, but from new_data_matrix.
		2010-4-22
			fix a bug. no_of_cols is not defined as it's moved into preEMMAX()
		2009-12-18
			fast-version EMMA
				estimate variance first by calling EMMA with only intercept. then call regular glm() with the new variance matrix.
		"""
		sys.stderr.write("Association by EMMAX (fast-version EMMA) ...\n")
		
		returnData = self.removeRowsWithNAPhenotypeFromKinshipAndSNPData(phenotype_ls, snpData, keywords.get('kinshipData'))
		new_data_matrix = returnData.snpData.data_matrix
		kinship_matrix = returnData.kinshipData.data_matrix
		
		preEMMAX_data = self.preEMMAX(returnData.non_NA_phenotype_ls, kinship_matrix, debug=self.debug)
		variance_matrix = preEMMAX_data.variance_matrix
		non_NA_phenotype_ar = preEMMAX_data.non_NA_phenotype_ar
		L_inverse = preEMMAX_data.L_inverse
		
		results = []
		counter = 0
		real_counter = 0
		
		no_of_rows, no_of_cols = new_data_matrix.shape	#2010-4-23
		for j in range(no_of_cols):
			genotype_ls = new_data_matrix[:,j]
			pdata = self.linear_model(genotype_ls, non_NA_phenotype_ar, min_data_point, snp_index=j, \
									kinship_matrix=None, eig_L=None, run_type=4, counting_and_NA_checking=True,\
									variance_matrix=None, lower_triangular_cholesky_inverse=L_inverse)
			# counting_and_NA_checking=True to get MAF and MAC. doesn't matter if it's False.
			if pdata is not None:
				results.append(pdata)
				real_counter += 1
			counter += 1
			if self.report and counter%2000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		if self.report:
			sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		
		sys.stderr.write("Done.\n")
		return results
	
	@classmethod
	def output_emma_results(cls, results, SNP_header, output_fname, log_pvalue=0):
		"""
		2011-2-17
			call returnSNPIdLs() to format the SNP id output
		2008-11-12
			this is written when rpy.r.emma_REML_t() is used in Emma_whole_matrix().
		"""
		#collect results
		sys.stderr.write("Outputting pvalue results ...")
		writer = csv.writer(open(output_fname,'w'), delimiter='\t')
		
		no_of_snps = len(results['ps'])
		counter = 0
		real_counter = 0
		for i in range(no_of_snps):
			SNPIdLs = cls.returnSNPIdLs(SNP_header, snp_index=i)	#2011-2-17
			
			pvalue = results['ps'][i][0]
			var_perc = results['genotype_var_perc'][i][0]
			
			coeff_list = [results['beta0_est'][i][0], results['beta1_est'][i][0]]
			
			#pdata = PassingData(snp_index=i, pvalue=pvalue, count_ls=[0.5,20], var_perc=var_perc, coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list)
			MAF = 0.5
			MAC = 20
			writer.writerow(SNPIdLs + [pvalue, MAF, MAC, var_perc] + coeff_list)	#2011-2-17
			
			real_counter += 1
			counter += 1
			if cls.report and counter%2000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		if cls.report:
			sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		
		del writer
		sys.stderr.write("Done.\n")
	
	@classmethod
	def output_lm_results(cls, results, SNP_header, output_fname, log_pvalue=0):
		"""
		2011-2-17
			call returnSNPIdLs() to format the SNP id output
		2010-3-22
			pdata.snp_index could be a tuple or list, which means a couple of SNPs.
			right now only handles first 2 SNPs.
		2009-9-19 fix a bug
			if coeff_pvalue is not None:	#2009-9-19 different from "if coeff_pvalue:", which would skip 0.0
		2008-11-11
			adapted from output_kw_results() of Kruskal_Wallis.py
		2008-05-27
			log10
		2008-05-21
			more stuff in results
			each kw result is wrapped in PassingData
		2008-02-14
		"""
		sys.stderr.write("Outputting pvalue results ...")
		writer = csv.writer(open(output_fname,'w'), delimiter='\t')
		for i in range(len(results)):
			pdata = results[i]
			pvalue = pdata.pvalue
			snp_index = pdata.snp_index
			SNPIdLs = cls.returnSNPIdLs(SNP_header, snp_index)	#2010-3-22, 2011-2-17
			if log_pvalue:
				if pvalue>0:
					pvalue = -math.log10(pvalue)
				else:
					pvalue = 'NA'
			MAC = min(pdata.count_ls)
			MAF = float(MAC)/sum(pdata.count_ls)
			coeff_and_pvalue_ls = []
			for j in range(len(pdata.coeff_list)):
				coeff = pdata.coeff_list[j]
				coeff_pvalue = pdata.coeff_p_value_list[j]
				coeff_and_pvalue = '%s'%coeff
				if coeff_pvalue is not None:	#2009-9-19 different from "if coeff_pvalue:", which would skip 0.0
					coeff_and_pvalue += ':%s'%coeff_pvalue
				coeff_and_pvalue_ls.append(coeff_and_pvalue)
			writer.writerow(SNPIdLs + [pvalue, MAF, MAC, pdata.var_perc] + coeff_and_pvalue_ls)	#2011-2-17
		del writer
		sys.stderr.write("Done.\n")
	
	@classmethod
	def output_multi_variate_lm_results(cls, results, writer, log_pvalue=0, run_genome_scan=True, \
									extraVariateNameLs = ['S_square', 'var_perc', ]):
		"""
		2010-4-23
			fix a bug that arises when beta_interval_delta_ls, coeff_p_value_list, coeff_list are of different length
		2010-4-18
			add argument extraVariateNameLs, a list of stats to be outputted right after the test-variate
		2010-4-5
			add argument run_genome_scan
				if True, the values for the last covaraite (test_variate_name) are put out in first 4 columns.
				otherwise, output in the original order.
		2010-3-23
			mainly for GWA.checkEpistasisInJBLabData() from misc.py
		"""
		#sys.stderr.write("Outputting multi-variate pvalue results ...")
		for i in range(len(results)):
			pdata = results[i]
			pvalue = pdata.pvalue
			snp_index = pdata.snp_index
			test_variate_name = snp_index
			S_square = getattr(pdata, 'S_square', None)
			beta_interval_delta_ls = getattr(pdata, 'beta_interval_delta_ls', None)
			if beta_interval_delta_ls is None:
				beta_interval_delta_ls = [-1]
			coeff_p_value_list = getattr(pdata, 'coeff_p_value_list', None)
			if coeff_p_value_list is None:
				coeff_p_value_list = [-1]
			coeff_list = getattr(pdata, 'coeff_list', None)
			if coeff_list is None:
				coeff_list = [-1]
			
			if run_genome_scan:
				row = [test_variate_name, coeff_p_value_list[-1], coeff_list[-1], beta_interval_delta_ls[-1]]
				no_of_coeffs = len(coeff_list)-1	# excluding the last, which is the test statistic
			else:
				row = []
				no_of_coeffs = len(coeff_list)
			for extra_variate_name in extraVariateNameLs:	#2010-4-18
				row.append(getattr(pdata, extra_variate_name, None))
			for j in range(no_of_coeffs):
				if len(coeff_p_value_list)>j:
					coeff_pvalue = coeff_p_value_list[j]
				else:
					coeff_pvalue = None
				if log_pvalue:
					if coeff_pvalue>0:
						coeff_pvalue = -math.log10(coeff_pvalue)
					else:
						coeff_pvalue = 'NA'
				if len(coeff_list)>j:
					coeff = coeff_list[j]
				else:
					coeff = None
				if len(beta_interval_delta_ls)>j:	#2010-4-23
					beta_interval_delta = beta_interval_delta_ls[j]
				else:
					beta_interval_delta = None
				row.extend([coeff_pvalue, coeff, beta_interval_delta])
			writer.writerow(row)
		#sys.stderr.write("Done.\n")
	
	@classmethod
	def removeUnPhenotypedSNPData(clf, snpData, header_phen, strain_acc_list_phen, data_matrix_phen, phenotype_method_id_ls):
		"""
		2011-4-27
			check if each row_id of snpData.row_id_ls is of type list/tuple or not.
		2010-4-21
			replace PlotGroupOfSNPs.findOutWhichPhenotypeColumn() with phenData.getColIndexLsGivenQuerySet()
		2010-2-25
			remove un-phenotyped ecotypes from the SNP data in order to keep the snp dataset small 
		"""
		sys.stderr.write("Removing un-phenotyped ecotypes from the SNP data ...")
		phenData = SNPData(header=header_phen, strain_acc_list=strain_acc_list_phen, data_matrix=data_matrix_phen)
		if phenotype_method_id_ls:
			which_phenotype_ls = phenData.getColIndexLsGivenQuerySet(Set(phenotype_method_id_ls), \
												colIDHashFunction=OutputPhenotype.extractPhenotypeIDFromMethodIDName)
		else:	#if not available, take all phenotypes
			which_phenotype_ls = range(len(phenData.col_id_ls))
		
		phenotyped_ecotype_id_set = set()
		for i in range(len(phenData.row_id_ls)):
			ecotype_id = phenData.row_id_ls[i]
			keep_this_ecotype = False
			for col_index in which_phenotype_ls:
				if phenData.data_matrix[i][col_index]!='NA':	# 2010-2-25 phenotype values are in raw string.
					keep_this_ecotype = True
					break
			if keep_this_ecotype:
				phenotyped_ecotype_id_set.add(ecotype_id)
		
		row_ids_to_be_kept = set()	# 2010-2-21
		no_of_ecotypes_in_total = len(snpData.row_id_ls)
		for row_id in snpData.row_id_ls:
			if type(row_id)==list or type(row_id)==tuple:	#2011-4-27
				ecotype_id = row_id[0]	#1st column is ecotype_id, 2nd is array id
			else:	#single-column
				ecotype_id = row_id
			if ecotype_id in phenotyped_ecotype_id_set:
				row_ids_to_be_kept.add(row_id)
		snpData = SNPData.keepRowsByRowID(snpData, row_ids_to_be_kept)
		no_of_removed = no_of_ecotypes_in_total - len(row_ids_to_be_kept)
		no_of_kept = len(row_ids_to_be_kept)
		sys.stderr.write("%s removed. %s kept. Done.\n"%(no_of_removed, no_of_kept))
		return snpData
	
	@classmethod
	def readInData(cls, phenotype_fname, input_fname, eigen_vector_fname=None, phenotype_method_id_ls=[], test_type=1, report=0,\
				snpAlleleOrdinalConversion=True, removeUnPhenotypedSNPData=True, ignore_2nd_column=False, kinship_fname=None,\
				genotype_fname_to_generate_kinship=None, needKinship=False):
		"""
		2011-4-27
			add argument kinship_fname and genotype_fname_to_generate_kinship
		2010-4-21
			add argument removeUnPhenotypedSNPData to offer choice whether to remove un-phenotyped ecotypes from snpData.
			add argument ignore_2nd_column, default is False. Ignore 2nd column in reading phenotype_fname and input_fname.
			The two arguments are used in PlotGroupOfSNPs.py
			replace PlotGroupOfSNPs.findOutWhichPhenotypeColumn() with phenData.getColIndexLsGivenQuerySet()
		2010-2-28
			set default of eigen_vector_fname to None, phenotype_method_id_ls to []
			add argument snpAlleleOrdinalConversion (encode SNP alleles into 0, 1, 2, ...). default=True.
				If a SNP is bi-allelic, it's binary after conversion.
		2010-2-25
			call removeUnPhenotypedSNPData() to shrink the snp dataset by removing un-phenotyped ecotypes
			pass snpData.strain_acc_list to get_phenotype_matrix_in_data_matrix_order()
		2009-3-20
			refactored out of run(), easy for MpiAssociation.py to call
		"""
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname)
		snpData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
						data_matrix=data_matrix, turn_into_array=1, ignore_2nd_column=ignore_2nd_column)
		
		if phenotype_fname and os.path.isfile(phenotype_fname):	# 2010-2-28 make sure it exists
			header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(phenotype_fname, turn_into_integer=0)
			if phenotype_method_id_ls and removeUnPhenotypedSNPData:	#2010-4-21
				snpData = cls.removeUnPhenotypedSNPData(snpData, header_phen, strain_acc_list_phen, data_matrix_phen, phenotype_method_id_ls)
			data_matrix_phen = cls.get_phenotype_matrix_in_data_matrix_order(snpData.strain_acc_list, strain_acc_list_phen, data_matrix_phen)
			# 2010-2-25 snpData.strain_acc_list is different from the old standalone strain_acc_list.
			phenData = SNPData(header=header_phen, strain_acc_list=snpData.strain_acc_list, data_matrix=data_matrix_phen,\
							ignore_2nd_column=ignore_2nd_column)
		else:
			phenData = None
		
		if snpAlleleOrdinalConversion or test_type==4:	# test_type 4 involves PCA_svd, which requires bi-allelic SNPs.
			newSnpData, allele2index_ls = snpData.convertSNPAllele2Index(report)	#0 (NA) or -2 (untouched) is all converted to -2 as 0 is used to denote allele
			newSnpData.header = snpData.header
			snpData = newSnpData	
		if eigen_vector_fname:
			PC_data = cls.getPCFromFile(eigen_vector_fname)
			PC_matrix = PC_data.PC_matrix
		else:
			if test_type==4:	#eigen_vector_fname not given for this test_type. calcualte PCs.
				T, P, explained_var = pca_module.PCA_svd(snpData.data_matrix, standardize=False)
				PC_matrix = T
			else:
				PC_matrix = None
		
		if phenData:	# 2010-2-28 make sure it exists
			if phenotype_method_id_ls:
				which_phenotype_ls = phenData.getColIndexLsGivenQuerySet(Set(phenotype_method_id_ls), \
																	colIDHashFunction=OutputPhenotype.extractPhenotypeIDFromMethodIDName)
			else:	#if not available, take all phenotypes
				which_phenotype_ls = range(len(phenData.col_id_ls))
		else:
			which_phenotype_ls = []
		
		if needKinship:
			#rows of snpData that are not in kinshipData will be removed. so a new snpData is returned.
			kinshipData, snpData, Z = cls.getExpandedKinship(kinship_fname=kinship_fname, \
										genotype_fname_to_generate_kinship=genotype_fname_to_generate_kinship, \
										independentSNPData=snpData)
			if snpData.row_id_ls!=phenData.row_id_ls:	#2011-4-27 snpData might have been changed by getExpandedKinship()
				#2011-4-27 the order of rows in snpData and phenData have to be same. get_phenotype_matrix_in_data_matrix_order() was run already.
				phenData = phenData.removeRowsNotInTargetSNPData(snpData)
		else:
			kinshipData = None
		pdata = PassingData(snpData=snpData, phenData=phenData, PC_matrix=PC_matrix, which_phenotype_ls=which_phenotype_ls, \
						phenotype_method_id_ls=phenotype_method_id_ls, kinshipData=kinshipData)
		sys.stderr.write("%s phenotypes, %s accessions, %s SNPs.\n"%(len(which_phenotype_ls), len(snpData.row_id_ls), \
															len(snpData.col_id_ls)))
		return pdata
	
	def preprocessForTwoPhenotypeAsso(self, snpData, which_phenotype_ls, data_matrix_phen):
		"""
		2009-3-20
			refactored out of run(), easy for MpiAssociation.py to call
		"""
		if len(which_phenotype_ls)<2:
			sys.stderr.write("Error: Require to specify 2 phenotypes in order to carry out this test type.\n")
			sys.exit(3)
		snpData.data_matrix = numpy.vstack((snpData.data_matrix, snpData.data_matrix))
		#stack the two phenotypes, fake a data_matrix_phen, which_phenotype_ls
		no_of_strains = len(strain_acc_list)
		which_phenotype1, which_phenotype2 = which_phenotype_ls[:2]
		phenotype1 = data_matrix_phen[:, which_phenotype1]
		phenotype1 = numpy.resize(phenotype1, [no_of_strains, 1])	#phenotype1.resize([10,1]) doesn't work. ValueError: 'resize only works on single-segment arrays'
		phenotype2 = data_matrix_phen[:, which_phenotype2]
		phenotype2 = numpy.resize(phenotype2, [no_of_strains, 1])
		data_matrix_phen = numpy.vstack((phenotype1, phenotype2))
		which_phenotype_index_ls = [0]
		#stack the PC_matrix as well
		if PC_matrix is not None:
			PC_matrix = numpy.vstack((PC_matrix, PC_matrix))
		
		#create an environment variable
		a = numpy.zeros(no_of_strains)
		a.resize([no_of_strains, 1])
		b = numpy.ones(no_of_strains)
		b.resize([no_of_strains, 1])
		environment_matrix = numpy.vstack((a, b))
		return which_phenotype_index_ls, environment_matrix, data_matrix_phen
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		if self.test_type in [3,7,8]:
			needKinship = True
		else:
			needKinship = False
		
		#2011-4-27 ignore_2nd_column has to be turned on all the time now.
		initData = self.readInData(self.phenotype_fname, self.input_fname, self.eigen_vector_fname, \
								phenotype_method_id_ls=self.phenotype_method_id_ls, \
								test_type=self.test_type, report=self.report, \
								snpAlleleOrdinalConversion=1-self.noSNPAlleleOrdinalConversion,\
								ignore_2nd_column=True, \
								kinship_fname=self.kinship_fname, \
								genotype_fname_to_generate_kinship=self.genotype_fname_to_generate_kinship,\
								needKinship=needKinship)
		
		if self.test_type==5 or self.test_type==6:
			which_phenotype_index_ls, environment_matrix, initData.phenData.data_matrix = \
					self.preprocessForTwoPhenotypeAsso(initData.snpData, \
													initData.which_phenotype_ls,\
													initData.phenData.data_matrix)
		else:
			which_phenotype_index_ls = initData.which_phenotype_ls
			environment_matrix = None
		
		if self.test_type==6:
			gene_environ_interaction = True
		else:
			gene_environ_interaction = False
		
		for which_phenotype in which_phenotype_index_ls:
			if self.test_type==5 or self.test_type==6:
				which_phenotype1, which_phenotype2 = initData.which_phenotype_ls[:2]
				phenotype1_name = initData.phenData.col_id_ls[which_phenotype1]
				phenotype2_name = initData.phenData.col_id_ls[which_phenotype2]
				phenotype_name = phenotype1_name+'_'+phenotype2_name
			else:
				phenotype_name = initData.phenData.col_id_ls[which_phenotype]
			phenotype_name = phenotype_name.replace('/', '_')	#'/' will be recognized as directory in output_fname
			if len(which_phenotype_index_ls)==1:	#2012.6.5 only one phenotype, don't bother change output_fname
				output_fname = self.output_fname
			else:
				output_fname='%s_pheno_%s.tsv'%(os.path.splitext(self.output_fname)[0], phenotype_name)	#make up a new name corresponding to this phenotype
			results = self.run_whole_matrix[self.test_type](initData.snpData, initData.phenData.data_matrix[:, which_phenotype], \
											self.min_data_point, PC_matrix=initData.PC_matrix, \
											which_PC_index_ls=self.which_PC_index_ls, environment_matrix=environment_matrix,\
											gene_environ_interaction=gene_environ_interaction, kinshipData=initData.kinshipData)
			output_results_func = self.output_results.get(self.test_type)
			if output_results_func is None:
				output_results_func = self.output_lm_results
			output_results_func(results, initData.snpData.col_id_ls, output_fname, self.minus_log_pvalue)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Association
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
