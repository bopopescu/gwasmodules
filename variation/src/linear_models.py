"""
Contains functions to perform various linear regression schemes, such as simple, and mixed models.
"""
try:
	import scipy as sp
	sp.alterdot()
	from scipy import linalg
	from scipy import stats
	from scipy import optimize
	#import scipy as sp
	anti_decoder = {1:0, 0:1}
	get_anti_snp = sp.vectorize(lambda x: anti_decoder[x])  #Creating a vectorized function for anti-snp
except Exception, err_str:
	print 'scipy is missing:', err_str
import snpsdata
import warnings
import pdb
import cPickle
import os
import sys
import time
from pdb import Pdb
import itertools as it
import math
import platform
import os
import analyze_gwas_results as agr
import ipdb

os.putenv('OMP_NUM_THREADS', '4')




def cholesky(V):
	try:
		H_sqrt = sp.mat(linalg.cholesky(V))
	except Exception, err_str:
		import rpy2
		from rpy2 import robjects
		robjects.r.library('Matrix')
		robjects.r("""
		f <- function(v,len){
			m <- matrix(v,len,len);
		  	return(as.matrix(nearPD(m)['mat'][[1]]));
		}
		""")
		import rpy2.robjects.numpy2ri
		near_V = sp.array(robjects.r.f(sp.array(V), len(V)))
		rel_change = (near_V - V) / V
		avg_rel_change = sp.average(sp.absolute(rel_change))
		max_rel_change = rel_change.max()
		print 'Average absolute relative change in matrix: %0.6f' % avg_rel_change
		print 'Maximum relative change in matrix: %0.6f' % max_rel_change
		H_sqrt = sp.mat(linalg.cholesky(near_V))
	return H_sqrt


def qr_decomp(X):
	if sp.__version__ >= '0.9':
		return linalg.qr(X, mode='economic')  #Do the QR-decomposition for the Gram-Schmidt process.			
	else:
		return linalg.qr(X, econ=True)  #Do the QR-decomposition for the Gram-Schmidt process.




class LinearModel(object):
	"""
	A simple linear model
	"""
	def __init__(self, Y=None):
		"""
		The fixed effects should be a list of fixed effect lists (SNPs)
		"""
		self.n = len(Y)
		self.Y = sp.matrix(Y).T
		self.X = sp.mat(sp.repeat(1, self.n)).T #The intercept
		self.p = 1
		self.beta_est = None
		self.cofactors = []



	def add_factor(self, x, lin_depend_thres=1e-8):
		"""
		Adds an explanatory variable to the X matrix.
		"""
		#Checking whether this new cofactor in linearly independent.
		new_x = sp.matrix(x).T
		(beta, rss, rank, sigma) = linalg.lstsq(self.X, new_x)
		if float(rss) < lin_depend_thres:
			warnings.warn('A factor was found to be linearly dependent on the factors already in the X matrix.  Hence skipping it!')
			return False
		self.X = sp.hstack([self.X, new_x])
		self.cofactors.append(x)
		self.p += 1
		return True


	def set_factors(self, factors, include_intercept=True):
		"""
		Set the cofactors.
		"""
		self.p = 0
		if include_intercept:
			self.X = sp.mat(sp.repeat(1, self.n), dtype='single').T #The intercept
			self.p = 1
			if len(factors) > 0:
				self.X = sp.hstack([self.X, sp.matrix(factors, dtype='single').T])
			self.p += len(factors)

		else:
			self.X = sp.matrix(factors, dtype='single').T
			self.p = len(factors)





	def get_hat_matrix(self):
		self.X_squared_inverse = (self.X.T * self.X).I
		self.hat_matrix = self.X * self.X_squared_inverse * self.X.T
		return self.hat_matrix



	def least_square_estimate(self):
		"""
		Via Normal equations, get LSEs
		"""
		self.X_squared_inverse = (self.X.T * self.X).I
		self.beta_est = self.X_squared_inverse * self.X.T * self.Y
		return self.beta_est

	def get_residuals(self):
		"""
		Calculates and returns the residuals as a column vector.
		"""
		#if self.beta_est==None:
		self.least_square_estimate()
		self.residuals = self.Y - self.X * self.beta_est
		return self.residuals


	def general_f_test(self, A, c):
		"""
		A general F-test implementation.
		Where the hypothesis is A*beta=c, a constraint.
		
		Here A is a matrix, and c a column vector		
		"""
		#if not self.residuals:
		self.get_residuals()
		q, p = shape(A)
		assert p == self.p, 'Shape of A is wrong!'
		B = (A * self.beta_est - c)
		M = A * (self.X_squared_inverse) * A.T
		f_stat = (B.T * M.I * B) / ((self.residuals.T * self.residuals) / (self.n - self.p))
		p_value = 1 - stats.f.cdf(f_stat, q, self.n - self.p)
		return p_value, f_stat



	def get_rss(self, dtype='single'):
		"""
		Returns the RSS
		"""
		h0_X = sp.mat(self.X, dtype=dtype)
		(betas, rss, r, s) = linalg.lstsq(h0_X, self.Y)
		return rss


	def get_ll(self, rss=None, dtype='single'):
		"""
		Returns the log-likelihood
		"""
		if not rss:
			rss = self.get_rss(dtype)
		return (-self.n / 2) * (1 + sp.log(2 * sp.pi) + sp.log(rss / self.n))


	def fast_f_test(self, snps, verbose=True, Z=None, with_betas=False):
		"""
		LM implementation 
		Single SNPs
					
		"""
		dtype = 'single'
		q = 1  # Single SNP is being tested
		p = len(self.X.T) + q
		n = self.n
		n_p = n - p
		num_snps = len(snps)

		h0_X = sp.mat(self.X, dtype=dtype)
		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, self.Y)
		Y = sp.mat(self.Y - h0_X * h0_betas, dtype=dtype)
		h0_betas = map(float, list(h0_betas))

		if not with_betas:
			(Q, R) = qr_decomp(h0_X)  #Do the QR-decomposition for the Gram-Schmidt process.
			Q = sp.mat(Q, dtype=dtype)
			Q2 = Q * Q.T
			M = sp.mat((sp.eye(n) - Q2), dtype=dtype)
		else:
			betas_list = [h0_betas] * num_snps

		rss_list = sp.repeat(h0_rss, num_snps)
		chunk_size = len(Y)
		for i in range(0, len(snps), chunk_size): #Do the dot-product in chuncks!
			snps_chunk = sp.matrix(snps[i:i + chunk_size])
			if with_betas:
				Xs = snps_chunk
			else:
				Xs = sp.mat(snps_chunk, dtype=dtype) * M
			for j in range(len(Xs)):
				if with_betas:
					(betas, rss, p, sigma) = linalg.lstsq(sp.hstack([h0_X, Xs[j].T]), Y, \
									overwrite_a=True)
					if not rss:
						if verbose: print 'No predictability in the marker, moving on...'
						continue
					betas_list[i + j] = map(float, list(betas))
				else:
					(betas, rss, p, sigma) = linalg.lstsq(Xs[j].T, Y, overwrite_a=True)
				rss_list[i + j] = rss[0]

				if verbose and num_snps >= 10 and (i + j + 1) % (num_snps / 10) == 0: #Print dots
					sys.stdout.write('.')
					sys.stdout.flush()

		if verbose and num_snps >= 10:
			sys.stdout.write('\n')
		rss_ratio = h0_rss / rss_list
		var_perc = 1 - 1 / rss_ratio
		f_stats = (rss_ratio - 1) * n_p / float(q)
		p_vals = stats.f.sf(f_stats, q, n_p)

		res_d = {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'var_perc':var_perc,
			'h0_rss':h0_rss, 'h0_betas':h0_betas}
		if with_betas:
			res_d['betas'] = betas_list
		return res_d


	def two_snps_ftest(self, snps, verbose=True):
		"""
		Every pair of SNPs, Vincent's results.
		"""
		import util
		num_snps = len(snps)

		ftest_res = self.fast_f_test(snps)
		full_rss = ftest_res['h0_rss']
		h0_X = self.X
		Y = self.Y	#The transformed outputs.

		#Contructing result containers		
		p3_f_stat_array = sp.zeros((num_snps, num_snps))
		p3_p_val_array = sp.ones((num_snps, num_snps))
		p4_f_stat_array = sp.zeros((num_snps, num_snps))
		p4_p_val_array = sp.ones((num_snps, num_snps))
		f_stat_array = sp.zeros((num_snps, num_snps))
		p_val_array = sp.ones((num_snps, num_snps))
		rss_array = sp.repeat(full_rss, num_snps * num_snps)
		rss_array.shape = (num_snps, num_snps)
		var_perc_array = sp.zeros((num_snps, num_snps))
		haplotype_counts = [[{} for j in range(i + 1)] for i in range(num_snps)]

		#Fill the diagonals with the single SNP emmax
		for i, snp in enumerate(snps):
			hap_set, hap_counts, haplotypes = snpsdata.get_haplotypes([snp], self.n,
										count_haplotypes=True)
			d = {'num_haps':hap_counts}
			for hap, hap_c in zip(hap_set, hap_counts):
				d[hap] = hap_c
			haplotype_counts[i][i] = d
			p_val_array[i, i] = ftest_res['ps'][i]
			p3_p_val_array[i, i] = p_val_array[i, i]
			p4_p_val_array[i, i] = p_val_array[i, i]
			f_stat_array[i, i] = ftest_res['f_stats'][i]
			p3_f_stat_array[i, i] = f_stat_array[i, i]
			p4_f_stat_array[i, i] = f_stat_array[i, i]
			rss_array[i, i] = ftest_res['rss'][i]
			var_perc_array[i, i] = ftest_res['var_perc'][i]


		identical_snp_count = 0
		no_interaction_count = 0
		for i, snp1 in enumerate(snps):
			snp1 = snps[i]
			for j in range(i):
				snp2 = snps[j]
				if i == j: continue #Skip diagonals..

				#Haplotype counts 
				hap_set, hap_counts, haplotypes = snpsdata.get_haplotypes([snp1, snp2], self.n,
											count_haplotypes=True)
				groups = set(haplotypes)
				d = {'num_haps':len(hap_counts)}
				for hap, hap_c in zip(hap_set, hap_counts):
					d[hap] = hap_c
				haplotype_counts[i][j] = d

				#Fill the upper right part with more powerful of two tests.

				if ftest_res['ps'][i] < ftest_res['ps'][j]:
					rss_array[j, i] = ftest_res['rss'][i]
					max_i = i
				else:
					rss_array[j, i] = ftest_res['rss'][j]
					max_i = j

				if d['num_haps'] == 2:
					identical_snp_count += 1
					continue
				elif d['num_haps'] == 3:
					n_p = self.n - 3
					no_interaction_count += 1
					#Do ANOVA
					l = []
					for g in groups:
						l.append(sp.int8(haplotypes == g))
					X = sp.mat(l)
					(betas, rss, p, sigma) = linalg.lstsq(X.T, Y)
					rss_array[i, j] = rss[0]
					var_perc_array[i, j] = 1 - rss / full_rss
					f_stat = (rss_array[j, i] / rss - 1) * n_p #FINISH
					p_val = stats.f.sf([f_stat], 1, n_p)[0]
					p3_f_stat_array[j, i] = f_stat
					p3_f_stat_array[i, j] = f_stat
					p3_p_val_array[j, i] = p_val
					p3_p_val_array[i, j] = p_val

					f_stat = ((full_rss - rss) / 2) / (rss / n_p)
					f_stat_array[j, i] = f_stat
					p_val_array[j, i] = stats.f.sf([f_stat], 2, n_p)[0]


				elif d['num_haps'] == 4: #With interaction
					n_p = self.n - 3
					#Do additive model
					snp_mat = sp.mat([snp1, snp2]).T #Transformed inputs
					X = sp.hstack([h0_X, snp_mat])
					(betas, rss, rank, s) = linalg.lstsq(X, Y)
					f_stat = (rss_array[j, i] / rss - 1) * n_p #Compared to only one SNP
					p_val = stats.f.sf([f_stat], 1, n_p)[0]
					rss_array[i, j] = rss
					p3_f_stat_array[j, i] = f_stat
					p3_p_val_array[j, i] = p_val

#					v_f_stat_array[j, i] = f_stat
#					v_p_val_array[j, i] = stats.f.sf([f_stat], 1, n_p)[0]

					f_stat = ((full_rss - rss) / 2) / (rss / n_p) #Compared to only the intercept
					f_stat_array[j, i] = f_stat
					p_val_array[j, i] = stats.f.sf([f_stat], 2, n_p)[0]

					#Generate the interaction, and test it.
					i_snp = snp1 & snp2
					snp_mat = sp.mat([i_snp]).T
					X = sp.hstack([h0_X, sp.mat([snps[max_i]]).T, snp_mat])
					(betas, rss, rank, s) = linalg.lstsq(X, Y)
					f_stat = (rss_array[j, i] / rss - 1) * n_p #Compared to only one SNP
					p_val = stats.f.sf([f_stat], 1, n_p)[0]
					p3_f_stat_array[i, j] = f_stat
					p3_p_val_array[i, j] = p_val


					#full model p-value
					n_p = self.n - 4
					l = []
					for g in groups:
						l.append(sp.int8(haplotypes == g))
					X = sp.mat(l)
					(betas, rss, p, sigma) = linalg.lstsq(X.T, Y)

					f_stat = ((rss_array[j, i] - rss) / 2) / (rss / n_p) #Compared to only one SNP
					p_val = stats.f.sf([f_stat], 2, n_p)[0]
					p4_f_stat_array[j, i] = f_stat
					p4_p_val_array[j, i] = p_val

					f_stat = (rss_array[i, j] / rss - 1) * n_p #Compared to two SNPs
					p4_f_stat_array[i, j] = f_stat
					p4_p_val_array[i, j] = stats.f.sf([f_stat], 1, n_p)[0]

					f_stat = ((full_rss - rss) / 3) / (rss / n_p) #Compared to only intercept
					f_stat_array[i, j] = f_stat
					p_val_array[i, j] = stats.f.sf([f_stat], 3, n_p)[0]
					rss_array[j, i] = rss

			if num_snps >= 10 and (i + 1) % (num_snps / 10) == 0: #Print dots
				sys.stdout.write('.')
				sys.stdout.flush()

		print no_interaction_count, identical_snp_count

		#FINISH res dict!!!
		res_dict = {'p3_ps':p3_p_val_array, 'p3_f_stats':p3_f_stat_array, 'p4_ps':p4_p_val_array,
			'p4_f_stats':p4_f_stat_array, 'rss':rss_array, 'var_perc':var_perc_array,
			'haplotype_counts':haplotype_counts,
			'f_stats':f_stat_array, 'ps':p_val_array}
		return res_dict



	def anova_f_test(self, snps, dtype='int8'):
		"""
		A standard ANOVA, using a F-test
		"""
		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(self.X, self.Y)
		num_snps = len(snps)
		rss_list = sp.repeat(h0_rss, num_snps)
		h0_betas = map(float, list(h0_betas)) + [0.0]
		betas_list = [h0_betas] * num_snps
		var_perc = sp.zeros(num_snps)
		f_stats = sp.zeros(num_snps)
		dfs = sp.zeros(num_snps)
		p_vals = sp.ones(num_snps)
		n = self.n
		p_0 = len(self.X.T)

		for i, snp in enumerate(snps):
			groups = sp.unique(snp)
			q = len(groups) - 1  # Null model has 1 df.
			p = p_0 + q
			n_p = n - p
			x = sp.zeros((len(groups), n), dtype=dtype)
			for g_i, g in enumerate(groups):
				x[g_i] = sp.int32(snp == g)
			(betas, rss, p, sigma) = linalg.lstsq(sp.mat(x).T, self.Y)

  			if not rss:
				print 'No predictability in the marker, moving on...'
				continue
			rss_list[i] = rss[0]
			betas_list[i] = map(float, list(betas))
			rss_ratio = h0_rss / rss
			var_perc[i] = 1 - 1 / rss_ratio
			f_stat = (rss_ratio - 1) * n_p / float(q)
			p_vals[i] = stats.f.sf([f_stat], q, n_p)[0]
			f_stats[i] = f_stat
			dfs[i] = n_p
			if num_snps >= 10 and (i + 1) % (num_snps / 10) == 0: #Print dots
				sys.stdout.write('.')
				sys.stdout.flush()

		if num_snps >= 10:
			sys.stdout.write('\n')

		return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list,
			'var_perc':var_perc, 'dfs':dfs}



	def anova_f_test_w_missing(self, snps):
		"""
		A standard ANOVA, using a F-test
		
		Handles SNPs w. missing data...
		"""
		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(self.X, self.Y)
		num_snps = len(snps)
		rss_list = sp.repeat(h0_rss, num_snps)
		h0_betas = map(float, list(h0_betas)) + [0.0]
		betas_list = [h0_betas] * num_snps
		var_perc = sp.zeros(num_snps)
		f_stats = sp.zeros(num_snps)
		dfs = sp.zeros(num_snps)
		p_vals = sp.ones(num_snps)
		n = self.n
		p_0 = len(self.X.T)

		for i, snp in enumerate(snps):
			groups = sp.unique(snp)
			q = len(groups) - 1  # Null model has 1 df.
			p = p_0 + q
			n_p = n - p
			x = []
			for g in groups:
				x.append(sp.int8(snp == g))
			(betas, rss, p, sigma) = linalg.lstsq(sp.mat(x).T, self.Y)

  			if not rss:
				print 'No predictability in the marker, moving on...'
				continue
			rss_list[i] = rss[0]
			betas_list[i] = map(float, list(betas))
			rss_ratio = h0_rss / rss
			var_perc[i] = 1 - 1 / rss_ratio
			f_stat = (rss_ratio - 1) * n_p / float(q)
			p_vals[i] = stats.f.sf([f_stat], q, n_p)[0]
			f_stats[i] = f_stat
			dfs[i] = n_p
			if num_snps >= 10 and (i + 1) % (num_snps / 10) == 0: #Print dots
				sys.stdout.write('.')
				sys.stdout.flush()

		if num_snps >= 10:
			sys.stdout.write('\n')

		return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list,
			'var_perc':var_perc, 'dfs':dfs}





	def test_explanatory_variable(self, x):
		"""
		Returns a p-value for whether adding this variable to the model explains the model better.
		
		Hopefully a sped-up version of a specific F-test.
		"""

		#THIS CAN BE SPED-UP MORE IF WE CHECK WHETHER self.X IS A VECTOR.  
		#AND USE t-TEST. 
		res_1 = self.get_residuals()

		X_mat = sp.hstack([self.X, sp.matrix([[v] for v in x])])
		X_sq = X_mat.T * X_mat
		try:
			X_sq_inv = X_sq.I
		except Exception, err_str:
			print err_str
			raise Exception('Andskotinn!!')

		res_2 = self.Y - X_mat * X_sq_inv * X_mat.T * self.Y
		rss_1 = res_1.T * res_1
		rss_2 = res_2.T * res_2
		f_stat = (rss_1 - rss_2) / (rss_2 / (self.n - self.p + 1))
		p_value = stats.f.sf(f_stat, 1, self.n - self.p + 1)
		return p_value, f_stat






class LinearMixedModel(LinearModel):
	"""
	A class for linear mixed models
	"""
	def __init__(self, Y=None, dtype='single'):
		"""
		The fixed effects should be a list of fixed effect lists (SNPs)
		"""
		self.n = len(Y)
		self.y_var = sp.var(Y, ddof=1, dtype=dtype)
		self.Y = sp.matrix([[y] for y in Y], dtype=dtype)
		self.X = sp.matrix([[1] for y in Y], dtype=dtype) #The intercept
		self.p = 1
		self.beta_est = None
		self.cofactors = []

		#A list of random effect type, and the cov matrix.
		self.random_effects = [('normal', sp.matrix(sp.identity(self.n)))] #The first random effect is the IID error.


	def add_random_effect(self, cov_matrix=None, effect_type='normal'):
		if effect_type != 'normal':
			raise Exception('Currently, only Normal random effects are allowed.')
		self.random_effects.append((effect_type, cov_matrix))

	def _get_eigen_L_(self, K=None, dtype='single'):
		if K == None:
			K = self.random_effects[1][1]
		if sp.__version__ < '0.8':
			K = sp.mat(K, dtype=dtype)
		evals, evecs = linalg.eigh(K)
		evals = sp.array(evals, dtype=dtype)
		return {'values':evals, 'vectors':sp.mat(evecs, dtype=dtype).T}



	def _get_eigen_R_(self, X=None, K=None, hat_matrix=None, dtype='single'):
		if X == None:
			X = self.X
		q = X.shape[1]
		if not hat_matrix:
			X_squared_inverse = linalg.pinv(X.T * X) #(X.T*X).I
			hat_matrix = X * X_squared_inverse * X.T
		if K == None:
			K = self.random_effects[1][1]
		S = sp.mat(sp.identity(self.n)) - hat_matrix	#S=I-X(X'X)^{-1}X'
		M = sp.mat(S * (K + self.random_effects[0][1]) * S, dtype='double')
		if sp.__version__ < '0.8':
			M = sp.mat(M, dtype=dtype)
		evals, evecs = linalg.eigh(M, overwrite_a=True) #eigen of S(K+I)S
		eig_values = sp.array(evals[q:], dtype=dtype) - 1  #Because of S(K+I)S?
		return {'values':eig_values, 'vectors':(sp.mat(evecs, dtype=dtype).T[q:])}


	def _rell_(self, delta, eig_vals, sq_etas):
		num_eig_vals = len(eig_vals)
		c_1 = 0.5 * num_eig_vals * (sp.log(num_eig_vals / (2.0 * sp.pi)) - 1)
		v = eig_vals + delta
		res = c_1 - 0.5 * (num_eig_vals * sp.log(sp.sum(sq_etas.flatten() / v)) + sp.sum(sp.log(v)))
		return res  #log-likelihoods (eq. 7 from paper)


	def _redll_(self, delta, eig_vals, sq_etas):
		num_eig_vals = len(eig_vals)
		v1 = eig_vals + delta
		v2 = sq_etas.flatten() / v1
		res = (num_eig_vals * sp.sum(v2 / v1) / sp.sum(v2) - sp.sum(1.0 / v1))
		return res  #diffrentiated log-likelihoods (*2) (eq. 9 from paper)


	def _ll_(self, delta, eig_vals, eig_vals_L, sq_etas):
		n = self.n
		c_1 = 0.5 * n * (sp.log(n / (2.0 * sp.pi)) - 1)
		v1 = eig_vals + delta
		v2 = eig_vals_L + delta
		res = c_1 - 0.5 * (n * sp.log(sp.sum(sq_etas.flatten() / v1)) + sp.sum(sp.log(v2)))
		return res  #log-likelihoods (eq. 6 from paper)


	def _dll_(self, delta, eig_vals, eig_vals_L, sq_etas):
		num_eig_vals = len(eig_vals)
		v1 = eig_vals + delta
		v2 = sq_etas.flatten() / v1
		v3 = eig_vals_L + delta
		res = (self.n * sp.sum(v2 / v1) / sp.sum(v2) - sp.sum(1.0 / v3))
		return res  #diffrentiated log-likelihoods (*2) (eq. 8 from paper)



	def get_REML(self, ngrids=100, llim= -10, ulim=10, esp=1e-6, eig_L=None, eig_R=None):
		"""
		Get REML estimates for the effect sizes, as well as the random effect contributions.
		
		This is EMMA
		"""
		if not eig_L:
			K = self.random_effects[1][1]
			eig_L = self._get_eigen_L_(K)

		#Get the variance estimates..
		res = self.get_estimates(eig_L=eig_L, ngrids=ngrids, llim=llim, ulim=ulim, esp=esp, method='REML',
					eig_R=eig_R)
		#res = self.get_estimates(ngrids=ngrids, llim=llim, ulim=ulim, esp=esp, method='REML', K=K)
		res['eig_L'] = eig_L
		return res



	def get_ML(self, ngrids=100, llim= -10, ulim=10, esp=1e-6, eig_L=None, eig_R=None):
		"""
		Get REML estimates for the effect sizes, as well as the random effect contributions.
		
		This is EMMA
		"""
		if not eig_L:
			K = self.random_effects[1][1]
			eig_L = self._get_eigen_L_(K)
		#Get the variance estimates..
		return self.get_estimates(eig_L=eig_L, ngrids=ngrids, llim=llim, ulim=ulim, esp=esp, method='ML',
					eig_R=eig_R)


	def get_estimates_3(self, xs=None, ngrids=[5, 5, 5, 5, 5, 5], llim= -5, ulim=5, method='REML',
				verbose=False, dtype='single', debug=False):
		"""
		Handles two K matrices, and one I matrix.
		
		Methods available are 'REML', and 'ML'		
		"""

		if verbose:
			print 'Retrieving %s variance estimates' % method
		if xs != None:
			X = sp.hstack([self.X, xs])
		else:
			X = self.X
#		xs1 = [[] for i in range(len(ngrids))]
#		ys1 = [[] for i in range(len(ngrids))]
#		xs2 = []
#		ys2 = []
		for it_i in range(len(ngrids)):
			delta = float(ulim - llim) / ngrids[it_i]
			#print llim, ulim
			#print delta
			log_k_ratio = llim
			lls = []
			res_list = []
			for i in range(ngrids[it_i] + 1):
				#xs1[it_i].append(log_k_ratio)
				#xs2.append(log_k_ratio)
				k_ratio = sp.exp(log_k_ratio)
				a = k_ratio / (k_ratio + 1.0)
				K = a * self.random_effects[1][1] + (1 - a) * self.random_effects[2][1]
				eig_L = self._get_eigen_L_(K=K)
				#Now perform EMMA
				res_dict = self.get_estimates(eig_L, K=K, xs=xs, ngrids=10, method=method, llim= -8, ulim=8)
				res_list.append(res_dict)
				lls.append(res_dict['max_ll'])
				#ys1[it_i].append(res_dict['max_ll'])
				#ys2.append(res_dict['pseudo_heritability'])
				log_k_ratio += delta
			max_ll_i = sp.argmax(lls)
			#print 'max_ll_i:', max_ll_i
			#print 'delta:', delta
			#Update the ulim and llim
			ulim = llim + delta * (max_ll_i + 1)
			llim = llim + delta * (max_ll_i - 1)
		opt_log_k_ratio = llim + delta * max_ll_i
		res_dict = res_list[max_ll_i]
		opt_k_ratio = sp.exp(log_k_ratio)
		a = opt_k_ratio / (opt_k_ratio + 1)
		res_dict['opt_k'] = scale_k(a * self.random_effects[1][1] + (1 - a) * self.random_effects[2][1])
		res_dict['opt_k_ratio'] = opt_k_ratio
		res_dict['perc_var1'] = a * res_dict['pseudo_heritability']
		res_dict['perc_var2'] = (1 - a) * res_dict['pseudo_heritability']
#		if debug:
#			import pylab
#			pylab.figure()
#			for i in range(3):
#				pylab.plot(xs1[i], ys1[i], marker='o', ls='None', alpha=0.4)
#			pylab.xlabel('log(K ratio)')
#			pylab.ylabel('log-likelihood')
#			pylab.savefig('/Users/bjarni.vilhjalmsson/tmp/%s.png' % str(debug))
#			pylab.figure()
#			pylab.plot(xs2, ys2, marker='o', ls='None', alpha=0.5)
#			pylab.xlabel('log(K ratio)')
#			pylab.ylabel('pseudo_heritability')
#			pylab.savefig('/Users/bjarni.vilhjalmsson/tmp/%s.png' % str(debug))
		return res_dict



	def get_estimates(self, eig_L, K=None, xs=None, ngrids=50, llim= -10, ulim=10, esp=1e-6,
				return_pvalue=False, return_f_stat=False, method='REML', verbose=False,
				dtype='single', eig_R=None):
		"""
		Get ML/REML estimates for the effect sizes, as well as the random effect contributions.
		Using the EMMA algorithm (Kang et al., Genetics, 2008).
		
		Methods available are 'REML', and 'ML'		
		"""
		if verbose:
			print 'Retrieving %s variance estimates' % method
		if xs != None:
			X = sp.hstack([self.X, xs])
		else:
			X = self.X

		if not (eig_R and xs != None):
			eig_R = self._get_eigen_R_(X=X, K=K)
		q = X.shape[1] #number of columns
		n = self.n
		p = n - q
		m = ngrids + 1

		etas = sp.array(eig_R['vectors'] * self.Y, dtype=dtype)
		sq_etas = etas * etas
		log_deltas = (sp.arange(m, dtype=dtype) / ngrids) * (ulim - llim) + llim  #a list of deltas to search
		deltas = sp.exp(log_deltas)
		assert len(deltas) == m, 'Number of deltas is incorrect.'
		eig_vals = sp.array(eig_R['values'], dtype=dtype)
		assert len(eig_vals) == p, 'Number of eigenvalues is incorrect.'

  		lambdas = sp.reshape(sp.repeat(eig_vals, m), (p, m)) + \
  			sp.reshape(sp.repeat(deltas, p), (m, p)).T
	  	s1 = sp.sum(sq_etas / lambdas, axis=0)
  		if method == 'REML':
	  		if verbose: print 'Calculating grid REML values'
	  		s2 = sp.sum(sp.log(lambdas), axis=0)
	  		lls = 0.5 * (p * (sp.log((p) / (2.0 * sp.pi)) - 1 - sp.log(s1)) - s2)
	 		s3 = sp.sum(sq_etas / (lambdas * lambdas), axis=0)
	 		s4 = sp.sum(1 / lambdas, axis=0)
	  		dlls = 0.5 * (p * s3 / s1 - s4)
  		elif method == 'ML':
	  		if verbose: print 'Calculating grid ML values'
  			#Xis < -matrix(eig.L$values, n, m) + matrix(delta, n, m, byrow=TRUE)
  			eig_vals_L = sp.array(eig_L['values'], dtype=dtype)
  			xis = sp.reshape(sp.repeat(eig_vals_L, m), (n, m)) + \
  				sp.reshape(sp.repeat(deltas, n), (m, n)).T
  			#LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Xis)))	
  			#dLL <- 0.5*delta*(n*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Xis))	

	  		s2 = sp.sum(sp.log(xis), axis=0)
	  		lls = 0.5 * (n * (sp.log((n) / (2.0 * sp.pi)) - 1 - sp.log(s1)) - s2)
	 		s3 = sp.sum(sq_etas / (lambdas * lambdas), axis=0)
	 		s4 = sp.sum(1 / xis, axis=0)
	  		dlls = 0.5 * (n * s3 / s1 - s4)

		max_ll_i = sp.argmax(lls)
		max_ll = lls[max_ll_i]

		last_dll = dlls[0]
		last_ll = lls[0]
		zero_intervals = []
		for i in range(1, len(dlls)):
			if dlls[i] < 0 and last_dll > 0:
				zero_intervals.append(((lls[i] + last_ll) * 0.5, i))
			last_ll = lls[i]
			last_dll = dlls[i]

		if len(zero_intervals) > 0:
			if verbose: print 'Found a zero interval... now performing Newton-Rhapson alg.'
			opt_ll, opt_i = max(zero_intervals)
			opt_delta = 0.5 * (deltas[opt_i - 1] + deltas[opt_i])
			#Newton-Raphson
			try:
				with warnings.catch_warnings():
	    				warnings.simplefilter("ignore")
	    	   			if method == 'REML':
						new_opt_delta = optimize.newton(self._redll_, opt_delta,
										args=(eig_vals, sq_etas),
										tol=esp, maxiter=100)
					elif method == 'ML':
						new_opt_delta = optimize.newton(self._dll_, opt_delta,
										args=(eig_vals, eig_vals_L, sq_etas),
										tol=esp, maxiter=100)
			except Exception, err_str:
				if verbose:
					print 'Problems with Newton-Raphson method:', err_str
					print "Using the maximum grid value instead."
					print 'opt_i:', opt_i
					print 'Grid optimal delta:', opt_delta
				new_opt_delta = opt_delta
			#Validating the delta
			if opt_i > 1 and deltas[opt_i - 1] - esp < new_opt_delta < deltas[opt_i] + esp:
				opt_delta = new_opt_delta
				opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
			#Cheking lower boundary
			elif opt_i == 1 and 0.0 < new_opt_delta < deltas[opt_i] + esp:
				opt_delta = new_opt_delta
				opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
			#Cheking upper boundary
			elif opt_i == len(deltas) - 1 and new_opt_delta > deltas[opt_i - 1] - esp \
						and not sp.isinf(new_opt_delta):
				opt_delta = new_opt_delta
				opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
			else:
				if verbose:
					print 'Local maximum outside of suggested possible areas?'
					print 'opt_i:', opt_i
					print 'Grid optimal delta:', opt_delta
					print "Newton's optimal delta:", new_opt_delta
					print 'Using the maximum grid value instead.'

			if verbose: print 'Done with Newton-Rahpson'
			if method == 'REML':
				opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
			elif method == 'ML':
				opt_ll = self._ll_(opt_delta, eig_vals, eig_vals_L, sq_etas)

			if opt_ll < max_ll:
				opt_delta = deltas[max_ll_i]
		else:
			if verbose: print 'No zero-interval was found, taking the maximum grid value.'
			opt_delta = deltas[max_ll_i]
			opt_ll = max_ll

		if verbose: print 'Finishing up.. calculating H_sqrt_inv.'
		l = sq_etas / (eig_vals + opt_delta)
		opt_vg = sp.sum(l) / p  #vg   
		opt_ve = opt_vg * opt_delta  #ve

		H_sqrt_inv = sp.mat(sp.diag(1.0 / sp.sqrt(eig_L['values'] + opt_delta)), dtype=dtype) * eig_L['vectors']
		#V = opt_vg * K + opt_ve * sp.eye(len(K))
		#H_sqrt = cholesky(V).T
		#H_sqrt_inv = H_sqrt.I
		X_t = H_sqrt_inv * X
		Y_t = H_sqrt_inv * self.Y
		(beta_est, mahalanobis_rss, rank, sigma) = linalg.lstsq(X_t, Y_t)
		x_beta = X * beta_est
		residuals = self.Y - x_beta
		rss = residuals.T * residuals
		x_beta_var = sp.var(x_beta, ddof=1)
		var_perc = x_beta_var / self.y_var
		res_dict = {'max_ll':opt_ll, 'delta':opt_delta, 'beta':beta_est, 've':opt_ve, 'vg':opt_vg,
			    'var_perc':var_perc, 'rss':rss, 'mahalanobis_rss':mahalanobis_rss, 'H_sqrt_inv':H_sqrt_inv,
			    'pseudo_heritability':1.0 / (1 + opt_delta)}

		if xs != None and return_f_stat:
#			rss_ratio = h0_rss / rss_list
#			var_perc = 1 - 1 / rss_ratio
#			f_stats = (rss_ratio - 1) * n_p / float(q)

			h0_X = H_sqrt_inv * self.X
			(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y_t)
			f_stat = (h0_rss / mahalanobis_rss - 1) * p / xs.shape[1]
			#f_stat = ((h0_rss - mahalanobis_rss) / (xs.shape[1])) / (h0_rss / p)

			res_dict['f_stat'] = float(f_stat)
		if return_pvalue:
			p_val = stats.f.sf(f_stat, (xs.shape[1]), p)
			res_dict['p_val'] = float(p_val)
		return res_dict #, lls, dlls, sp.log(deltas)



	def expedited_REML_t_test(self, snps, ngrids=50, llim= -4, ulim=10, esp=1e-6, verbose=True, eig_L=None):
		"""
		Single SNP analysis, i.e. EMMA, but faster than R EMMA.
		"""
		assert len(self.random_effects) == 2, "Expedited REMLE only works when we have exactly two random effects."
		K = self.random_effects[1][1]
		if eig_L == None:
			eig_L = self._get_eigen_L_(K)
		num_snps = len(snps)
		f_stats = sp.empty(num_snps)
		vgs = sp.empty(num_snps)
		ves = sp.empty(num_snps)
		max_lls = sp.empty(num_snps)
		var_perc = sp.empty(num_snps)
		rss_list = sp.empty(num_snps)
		betas = []
		p_vals = sp.empty(num_snps)

		for i, snp in enumerate(snps):
			res = self.get_estimates(eig_L=eig_L, xs=sp.matrix(snp).T, ngrids=ngrids, llim=llim, ulim=ulim,
						       esp=esp, return_pvalue=True, return_f_stat=True)
			f_stats[i] = res['f_stat']
			vgs[i] = res['vg']
			ves[i] = res['ve']
			max_lls[i] = res['max_ll']
			var_perc[i] = res['var_perc']
			betas.append(map(float, list(res['beta'])))
			p_vals[i] = res['p_val']
			rss_list[i] = res['rss']
			if verbose and num_snps >= 1000 and (i + 1) % (num_snps / 1000) == 0: #Print dots
				sys.stdout.write('.')
				sys.stdout.flush()


		return {'ps':p_vals, 'f_stats':f_stats, 'vgs':vgs, 'ves':ves, 'var_perc':var_perc,
			'max_lls':max_lls, 'betas':betas, 'rss':rss_list}



	def emmax_f_test_w_interactions(self, snps, int_af_threshold=15):
		"""
		EMMAX implementation (in python)
		Single SNPs
		
		With interactions between SNP and possible cofactors.
		"""
		assert len(self.random_effects) == 2, "Expedited REMLE only works when we have exactly two random effects."
		p_0 = len(self.X.T)
		n = self.n


		K = self.random_effects[1][1]
		eig_L = self._get_eigen_L_(K)
		res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
		delta = res['delta']
		print 'pseudo_heritability:', 1.0 / (1 + delta)
		H_sqr = res['H_sqrt']
		H_sqrt_inv = H_sqr.I
		Y = H_sqrt_inv * self.Y	#The transformed outputs.
		h0_X = H_sqrt_inv * self.X
		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
		h0_betas = map(float, list(h0_betas))
		f_stats = []
		rss_list = []
		betas_list = []
		p_vals = []
		var_perc = []
		cofactors = sp.array(self.cofactors)
		num_cof = len(cofactors)
		low_int_af_count = 0
		for snp in snps:
			snp_a = sp.array(snp)
			if sum(snp_a) > int_af_threshold:
				interactions = cofactors * snp_a
				interactions = interactions[sum(interactions, axis=1) > int_af_threshold]
				low_int_af_count += num_cof - len(interactions)
				X = sp.hstack([h0_X, H_sqrt_inv * sp.matrix(snp_a).T, H_sqrt_inv * sp.matrix(interactions).T])
			else:
				low_int_af_count += num_cof
				X = sp.hstack([h0_X, H_sqrt_inv * sp.matrix(snp_a).T])
			(betas, rss, p, sigma) = linalg.lstsq(X, Y)
			q = p - p_0
			n_p = n - p
			if not rss:
				if q == 0:
					print 'No predictability in the marker, moving on...'
					p_vals.append(1)
					f_stats.append(0)
					rss_list.append(h0_rss)
					betas_list.append(h0_betas)
					var_perc.append(0)
					continue
				else:
					res = (Y - X * betas)
					rss = res.T * res
			f_stat = ((h0_rss - rss) / q) / (rss / n_p)
			p_val = stats.f.sf(f_stat, q, n_p)
			p_vals.append(p_val[0])
			f_stats.append(f_stat[0])
			rss_list.append(rss[0])
			betas_list.append(map(float, list(betas)))
			var_perc.append(float(1 - rss / h0_rss))

		print 'low_int_af_count:', low_int_af_count

		return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list,
			'delta':delta, 'pseudo_heritability': 1.0 / (1 + delta), 'var_perc':var_perc}



	def emmax_anova_f_test(self, snps):
		"""
		EMMAX implementation (in python)
		Single SNPs
		
		With interactions between SNP and possible cofactors.
		"""
		K = self.random_effects[1][1]
		eig_L = self._get_eigen_L_(K)
		res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
		print 'pseudo_heritability:', res['pseudo_heritability']

		r = self._emmax_anova_f_test_(snps, res['H_sqrt'])
		r['pseudo_heritability'] = res['pseudo_heritability']
		r['max_ll'] = res['max_ll']
		return r



	def _emmax_anova_f_test_(self, snps, H_sqrt, verbose=True):
		"""
		EMMAX implementation (in python)
		Single SNPs
		
		With interactions between SNP and possible cofactors.
		"""
		n = self.n
		p_0 = len(self.X.T)

		H_sqrt_inv = H_sqrt.I
		Y = H_sqrt_inv * self.Y	#The transformed outputs.
		h0_X = H_sqrt_inv * self.X
		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
		h0_betas = map(float, list(h0_betas))
		num_snps = len(snps)
		rss_list = sp.repeat(h0_rss, num_snps)
		betas_list = [h0_betas] * num_snps
		var_perc = sp.zeros(num_snps)
		f_stats = sp.zeros(num_snps)
		dfs = sp.zeros(num_snps)
		p_vals = sp.ones(num_snps)

		for i, snp in enumerate(snps):
			groups = sp.unique(snp)
			q = len(groups) - 1  # Null model has 1 df.
			p = p_0 + q
			n_p = n - p
			l = []
			for g in groups:
				l.append(sp.int8(snp == g))
			X = sp.mat(l) * (H_sqrt_inv.T)
			(betas, rss, p, sigma) = linalg.lstsq(X.T, Y, overwrite_a=True)
			if not rss:
				if verbose: print 'No predictability in the marker, moving on...'
				continue
			rss_list[i] = rss[0]
			betas_list[i] = map(float, list(betas))
			rss_ratio = h0_rss / rss
			var_perc[i] = 1 - 1 / rss_ratio
			f_stat = (rss_ratio - 1) * n_p / float(q)
			p_vals[i] = stats.f.sf([f_stat], q, n_p)[0]
			f_stats[i] = f_stat
			dfs[i] = n_p
			if num_snps >= 10 and (i + 1) % (num_snps / 10) == 0: #Print dots
				sys.stdout.write('.')
				sys.stdout.flush()

		if num_snps >= 10:
			sys.stdout.write('\n')
		#rss_ratio = h0_rss / rss_list
		#var_perc = 1 - 1 / rss_ratio
		#f_stats = (rss_ratio - 1) * n_p / float(q)
		#p_vals = stats.f.sf(f_stats, q, n_p)


		return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list,
			'var_perc':var_perc, 'dfs':dfs}



	def emmax_permutations(self, snps, num_perm, method='REML'):
		"""
		EMMAX permutation test
		Single SNPs
		
		Returns the list of max_pvals and max_fstats 
		"""
		K = self.random_effects[1][1]
		eig_L = self._get_eigen_L_(K)
		#s = time.time()
		res = self.get_estimates(eig_L=eig_L, method=method) #Get the variance estimates..
		#print 'Took % .6f secs.' % (time.time() - s)
		#print 'pseudo_heritability:', res['pseudo_heritability']
		q = 1  # Single SNP is being tested
		p = len(self.X.T) + q
		n = self.n
		n_p = n - p
		H_sqrt_inv = res['H_sqrt_inv']
		Y = H_sqrt_inv * self.Y	#The transformed outputs.
		h0_X = H_sqrt_inv * self.X
		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
		Y = Y - h0_X * h0_betas
		num_snps = len(snps)
		max_fstat_list = []
		min_pval_list = []
		chunk_size = len(Y)
		Ys = sp.mat(sp.zeros((chunk_size, num_perm)))
		for perm_i in range(num_perm):
			#print 'Permutation nr. % d' % perm_i
			sp.random.shuffle(Y)
			Ys[:, perm_i] = Y

		min_rss_list = sp.repeat(h0_rss, num_perm)
		for i in range(0, num_snps, chunk_size): #Do the dot-product in chuncks!
			snps_chunk = sp.matrix(snps[i:i + chunk_size])
			Xs = snps_chunk * (H_sqrt_inv.T)
			Xs = Xs - sp.mat(sp.mean(Xs, axis=1))
			for j in range(len(Xs)):
				(betas, rss_list, p, sigma) = linalg.lstsq(Xs[j].T, Ys, overwrite_a=True)
				for k, rss in enumerate(rss_list):
					if not rss:
						print 'No predictability in the marker, moving on...'
						continue
					if min_rss_list[k] > rss:
						min_rss_list[k] = rss
				if num_snps >= 10 and (i + j + 1) % (num_snps / num_perm) == 0: #Print dots
					sys.stdout.write('.')
					sys.stdout.flush()

		if num_snps >= 10:
			sys.stdout.write('\n')
		min_rss = min(rss_list)
		max_f_stats = ((h0_rss / min_rss_list) - 1.0) * n_p / float(q)
		min_pvals = (stats.f.sf(max_f_stats, q, n_p))



		res_d = {'min_ps':min_pvals, 'max_f_stats':max_f_stats}
		return res_d


	def emmax_f_test(self, snps, snp_priors=None, Z=None, with_betas=False, method='REML', return_stepw_stats=False,
			eig_L=None, eig_R=None, emma_num=200, progress_file_writer=None):
		"""
		EMMAX implementation (in python)
		Single SNPs
		
		Returns various stats useful for stepwise regression if return_stepw_stats flag is set to True.
		"""
		if not eig_L:
			print 'Calculating the eigenvalues of K'
			eig_L = self._get_eigen_L_()
			print 'Done.'
		if not eig_R:
			print "Calculating the eigenvalues of S(K+I)S where S = I-X(X'X)^-1X'"
			eig_R = self._get_eigen_R_(X=self.X)
			print 'Done'

		print 'Getting variance estimates'
		res = self.get_estimates(eig_L, method=method, eig_R=eig_R) #Get the variance estimates..
		print 'Done.'
		print 'pseudo_heritability:', res['pseudo_heritability']

		r = self._emmax_f_test_(snps, res['H_sqrt_inv'], snp_priors=snp_priors, Z=Z, with_betas=with_betas,
				progress_file_writer=progress_file_writer, emma_num=emma_num, eig_L=eig_L)
		r['pseudo_heritability'] = res['pseudo_heritability']
		r['ve'] = res['ve']
		r['vg'] = res['vg']
		r['max_ll'] = res['max_ll']

		return r




	def _emmax_f_test_(self, snps, H_sqrt_inv, snp_priors=None, verbose=True, return_transformed_snps=False,
			Z=None, with_betas=False, progress_file_writer=None, emma_num=200, eig_L=None):
		"""
		EMMAX implementation (in python)
		Single SNPs
		
		Methods:
			normal - Normal regression
			qr - Uses QR decomposition to speed up regression with many co-factors.
			
		"""
		dtype = 'single'
		q = 1  # Single SNP is being tested
		p = len(self.X.T) + q
		n = self.n
		n_p = n - p
		num_snps = len(snps)

		h0_X = sp.mat(H_sqrt_inv * self.X, dtype=dtype)
		Y = H_sqrt_inv * self.Y	#The transformed outputs.
		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
		Y = sp.mat(Y - h0_X * h0_betas, dtype=dtype)
		h0_betas = map(float, list(h0_betas))

		if Z != None:
			H_sqrt_inv = H_sqrt_inv * Z

		if not with_betas:
			(Q, R) = qr_decomp(h0_X)  #Do the QR-decomposition for the Gram-Schmidt process.
			Q = sp.mat(Q)
			Q2 = Q * Q.T
			M = sp.mat(H_sqrt_inv.T * (sp.eye(n) - Q2), dtype=dtype)
		else:
			betas_list = [h0_betas] * num_snps
			M = H_sqrt_inv.T

		rss_list = sp.repeat(h0_rss, num_snps)
		if return_transformed_snps:
			t_snps = []
		if snp_priors != None:
			snp_priors = sp.array(snp_priors)
			log_h0_rss = sp.log(h0_rss)
			log_bfs = sp.zeros(num_snps) #Bayes factors
		chunk_size = len(Y)
		if not progress_file_writer == None:
			progress_file_writer.update_progress_bar(progress=0.45, task_status='Performing EMMAX')
		for i in range(0, num_snps, chunk_size): #Do the dot-product in chuncks!
			snps_chunk = sp.matrix(snps[i:i + chunk_size], dtype=dtype)
			Xs = snps_chunk * M
			for j, X_j in enumerate(Xs):
				if return_transformed_snps:
					t_snps.append(sp.array(X_j).flatten())
				if with_betas:
					(betas, rss, p, sigma) = linalg.lstsq(sp.hstack([h0_X, X_j.T]), Y, \
									overwrite_a=True)
					if rss:
						betas_list[i + j] = map(float, list(betas))
				else:
					(betas, rss, p, sigma) = linalg.lstsq(X_j.T, Y, overwrite_a=True)
				if rss:
					rss_list[i + j] = rss[0]
					if snp_priors != None:
						log_bfs[i + j] = log_h0_rss - sp.log(rss) #-(1/2)*log(n)

				if verbose and num_snps >= 10 and (i + j + 1) % (num_snps / 10) == 0: #Print dots
					if not progress_file_writer == None:
						progress_file_writer.update_progress_bar(task_status='Performing EMMAX (SNPs completed: %d %%)' % round((float(i + j + 1) / num_snps) * 100))
					sys.stdout.write('.')
					sys.stdout.flush()

		if verbose and num_snps >= 10:
			sys.stdout.write('\n')


		rss_ratio = h0_rss / rss_list
		var_perc = 1 - 1 / rss_ratio
		f_stats = (rss_ratio - 1) * n_p / float(q)
		p_vals = stats.f.sf(f_stats, q, n_p)

		res_d = {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'var_perc':var_perc,
			'h0_rss':h0_rss, 'h0_betas':h0_betas}
		if with_betas:
			res_d['betas'] = betas_list
		if return_transformed_snps:
			res_d['t_snps'] = t_snps
		if snp_priors != None:
			bfs = sp.exp((log_bfs * n - sp.log(n)) * 1 / 2) #Bayes factors
			res_d['bfs'] = bfs
			pos = bfs * snp_priors / (1 - snp_priors) #Posterior odds
			res_d['pos'] = pos
			ppas = pos / (1 + pos) #Posterior probablities of association
			res_d['ppas'] = ppas

		if emma_num > 0:
			print 'Updating p-values using EMMA for the smallest %d p-values.' % emma_num
			pval_indices = sorted(zip(res_d['ps'], range(len(snps))))[:emma_num]
			l = map(list, zip(*pval_indices))
			top_snps = [snps[pi] for pi in l[1]]
			top_emma_res = self.expedited_REML_t_test(top_snps, eig_L=eig_L)
			for pi, p, f, r, v in it.izip(l[1], top_emma_res['ps'], top_emma_res['f_stats'], \
						top_emma_res['rss'], top_emma_res['var_perc']):
				res_d['ps'][pi] = p
				res_d['f_stats'][pi] = f
				res_d['rss'][pi] = r
				res_d['var_perc'][pi] = v

		return res_d



	def emmax_two_snps(self, snps, verbose=True):
		"""
		Every pair of SNPs, Vincent's results.
		"""
		import util
		K = self.random_effects[1][1]
		eig_L = self._get_eigen_L_(K)
		num_snps = len(snps)

		res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
		delta = res['delta']
		pseudo_her = 1.0 / (1 + delta)
		H_sqrt = res['H_sqrt']
		H_sqrt_inv = H_sqrt.I
		emmax_res = self._emmax_f_test_(snps, H_sqrt, verbose=False, return_transformed_snps=True)
		t_snps = emmax_res['t_snps']
		full_rss = emmax_res['h0_rss']
		h0_X = H_sqrt_inv * self.X
		Y = H_sqrt_inv * self.Y	#The transformed outputs.

		#Contructing result containers		
		p3_f_stat_array = sp.zeros((num_snps, num_snps))
		p3_p_val_array = sp.ones((num_snps, num_snps))
		p4_f_stat_array = sp.zeros((num_snps, num_snps))
		p4_p_val_array = sp.ones((num_snps, num_snps))
		f_stat_array = sp.zeros((num_snps, num_snps))
		p_val_array = sp.ones((num_snps, num_snps))
		rss_array = sp.repeat(full_rss, num_snps * num_snps)
		rss_array.shape = (num_snps, num_snps)
		var_perc_array = sp.zeros((num_snps, num_snps))
		haplotype_counts = [[{} for j in range(i + 1)] for i in range(num_snps)]

		#Fill the diagonals with the single SNP emmax
		for i, snp in enumerate(snps):
			hap_set, hap_counts, haplotypes = snpsdata.get_haplotypes([snp], self.n,
										count_haplotypes=True)
			d = {'num_haps':hap_counts}
			for hap, hap_c in zip(hap_set, hap_counts):
				d[hap] = hap_c
			haplotype_counts[i][i] = d
			p_val_array[i, i] = emmax_res['ps'][i]
			p3_p_val_array[i, i] = p_val_array[i, i]
			p4_p_val_array[i, i] = p_val_array[i, i]
			f_stat_array[i, i] = emmax_res['f_stats'][i]
			p3_f_stat_array[i, i] = f_stat_array[i, i]
			p4_f_stat_array[i, i] = f_stat_array[i, i]
			rss_array[i, i] = emmax_res['rss'][i]
			var_perc_array[i, i] = emmax_res['var_perc'][i]


		identical_snp_count = 0
		no_interaction_count = 0
#		p = len(self.X.T) + q + 1 #Adding one SNP as a cofactor.
#		n_p = self.n - p
		for i, snp1 in enumerate(snps):
			snp1 = snps[i]
			for j in range(i):
				snp2 = snps[j]
				if i == j: continue #Skip diagonals..

				#Haplotype counts 
				hap_set, hap_counts, haplotypes = snpsdata.get_haplotypes([snp1, snp2], self.n,
											count_haplotypes=True)
				groups = set(haplotypes)
				d = {'num_haps':len(hap_counts)}
				for hap, hap_c in zip(hap_set, hap_counts):
					d[hap] = hap_c
				haplotype_counts[i][j] = d

				#Fill the upper right part with more powerful of two tests.

				if emmax_res['ps'][i] < emmax_res['ps'][j]:
					rss_array[j, i] = emmax_res['rss'][i]
					max_i = i
				else:
					rss_array[j, i] = emmax_res['rss'][j]
					max_i = j

				if d['num_haps'] == 2:
					identical_snp_count += 1
					continue
				elif d['num_haps'] == 3:
					n_p = self.n - 3
					no_interaction_count += 1
					#Do ANOVA
					l = []
					for g in groups:
						l.append(sp.int8(haplotypes == g))
					X = sp.mat(l) * (H_sqrt_inv.T)
					(betas, rss, p, sigma) = linalg.lstsq(X.T, Y)
					rss_array[i, j] = rss[0]
					var_perc_array[i, j] = 1 - rss / full_rss
					f_stat = (rss_array[j, i] / rss - 1) * n_p #FINISH
					p_val = stats.f.sf([f_stat], 1, n_p)[0]
					p3_f_stat_array[j, i] = f_stat
					p3_f_stat_array[i, j] = f_stat
					p3_p_val_array[j, i] = p_val
					p3_p_val_array[i, j] = p_val

					f_stat = ((full_rss - rss) / 2) / (rss / n_p)
					f_stat_array[j, i] = f_stat
					p_val_array[j, i] = stats.f.sf([f_stat], 2, n_p)[0]


				elif d['num_haps'] == 4: #With interaction
					n_p = self.n - 3
					#Do additive model
					#snp_mat = H_sqrt_inv * sp.mat([snp1, snp2]).T #Transformed inputs
					snp_mat = sp.mat([t_snps[i], t_snps[j]]).T #Transformed inputs
					X = sp.hstack([h0_X, snp_mat])
					(betas, rss, rank, s) = linalg.lstsq(X, Y)
					f_stat = (rss_array[j, i] / rss - 1) * n_p #Compared to only one SNP
					p_val = stats.f.sf([f_stat], 1, n_p)[0]
					rss_array[i, j] = rss
					p3_f_stat_array[j, i] = f_stat
					p3_p_val_array[j, i] = p_val

#					v_f_stat_array[j, i] = f_stat
#					v_p_val_array[j, i] = stats.f.sf([f_stat], 1, n_p)[0]

					f_stat = ((full_rss - rss) / 2) / (rss / n_p) #Compared to only the intercept
					f_stat_array[j, i] = f_stat
					p_val_array[j, i] = stats.f.sf([f_stat], 2, n_p)[0]

					#Generate the interaction, and test it.
					i_snp = snp1 & snp2
					snp_mat = H_sqrt_inv * sp.mat([i_snp]).T
					X = sp.hstack([h0_X, sp.mat([t_snps[max_i]]).T, snp_mat])
					(betas, rss, rank, s) = linalg.lstsq(X, Y)
					f_stat = (rss_array[j, i] / rss - 1) * n_p #Compared to only one SNP
					p_val = stats.f.sf([f_stat], 1, n_p)[0]
					p3_f_stat_array[i, j] = f_stat
					p3_p_val_array[i, j] = p_val


					#full model p-value
					n_p = self.n - 4
					l = []
					for g in groups:
						l.append(sp.int8(haplotypes == g))
					X = sp.mat(l) * (H_sqrt_inv.T)
					(betas, rss, p, sigma) = linalg.lstsq(X.T, Y)

					f_stat = ((rss_array[j, i] - rss) / 2) / (rss / n_p) #Compared to only one SNP
					p_val = stats.f.sf([f_stat], 2, n_p)[0]
					p4_f_stat_array[j, i] = f_stat
					p4_p_val_array[j, i] = p_val

					f_stat = (rss_array[i, j] / rss - 1) * n_p #Compared to two SNPs
					p4_f_stat_array[i, j] = f_stat
					p4_p_val_array[i, j] = stats.f.sf([f_stat], 1, n_p)[0]

					f_stat = ((full_rss - rss) / 3) / (rss / n_p) #Compared to only intercept
					f_stat_array[i, j] = f_stat
					p_val_array[i, j] = stats.f.sf([f_stat], 3, n_p)[0]
					rss_array[j, i] = rss

			if num_snps >= 10 and (i + 1) % (num_snps / 10) == 0: #Print dots
				sys.stdout.write('.')
				sys.stdout.flush()

		print no_interaction_count, identical_snp_count

		#FINISH res dict!!!
		res_dict = {'p3_ps':p3_p_val_array, 'p3_f_stats':p3_f_stat_array, 'p4_ps':p4_p_val_array,
			'p4_f_stats':p4_f_stat_array, 'rss':rss_array, 'var_perc':var_perc_array,
			'pseudo_heritability': pseudo_her, 'haplotype_counts':haplotype_counts,
			'f_stats':f_stat_array, 'ps':p_val_array}
		return res_dict


	def calc_statistics(self):
		"""
		Returns all sorts of statistics used in stepwise regression
		
		Log Likelihood, BIC, modified BIC, extended BIC, RSS, mahalnobis RSS 
		"""
		pass


class MultipleTraitLMM(LinearMixedModel):
	"""
	A class to encompass multiple traits/environments mixed model.
	"""
	def __init__(self, phenotype_list, ecotype_list, K=None, K_ets=None, unique_ets=None):
		"""
		The fixed effects should be a list of fixed effect lists (SNPs)
		"""
		joint_phen_vals = []
		if phenotype_list:
			for phen in phenotype_list:
				joint_phen_vals.extend(phen)
		self.n = len(joint_phen_vals)
		self.Y = sp.matrix(joint_phen_vals).T

		self.joint_ets = []
		if ecotype_list:
			for ets in ecotype_list:
				self.joint_ets.append(ets)
		if not unique_ets:
			unique_ets = list(set(self.joint_ets))
		self.unique_ets = unique_ets

		self.ys = phenotype_list
		self.et_list = ecotype_list
		self.num_ets = map(len, ecotype_list)
		self.corr_mat = self.get_phen_corr(phenotype_list, ecotype_list)
		self.X = sp.mat(sp.repeat(1, self.n)).T #The intercept
		self.num_unique = len(unique_ets)

		print 'Creating ecotype maps'
		self.ecotype_maps = []
		joint_ets = []
		for ets, phens in zip(ecotype_list, phenotype_list):
			l = []
			for et in ets:
				l.append(unique_ets.index(et))
			self.ecotype_maps.append(l)
			joint_ets.extend(ets)

		self.Z = sp.int8(sp.mat(joint_ets).T == sp.mat(unique_ets))
		self.p = len(phenotype_list)

		if K != None and K_ets != None:
			K = prepare_k(K, K_ets, unique_ets)
			self.K = scale_k(K)


		#E in the model.
		E = []
		shift = self.num_ets[0]
		for num_ets in self.num_ets[1:]: #Skip the first one..
			E.append([0] * shift + [1] * num_ets + [0] * (self.n - num_ets - shift))
			shift += num_ets
		self.E = sp.mat(E).T





	def get_marginal_estimates(self):
		print 'Retrieving marginal estimates'
		r_list = []
		for y, ets_map in zip(self.ys, self.ecotype_maps):
			K = self.K[ets_map, :][:, ets_map] #The K matrix
			X = self.X[ets_map].T #The cofactors...
			res = get_emma_reml_estimates(y, K, X, include_intercept=False)
			r_list.append(res)

		return r_list


	def get_variance_matrix(self):
		"""
		Assumes kinship is there...
		"""
		gen_var_list = []
		err_var_list = []
		her_list = []
		r_list = self.get_marginal_estimates()
		for res in r_list:
			her_list.append(res['pseudo_heritability'])
			err_var_list.append(res['ve'])
			gen_var_list.append(res['vg'])
		print err_var_list, gen_var_list, her_list
		V = sp.zeros((self.n, self.n))
		m_i = 0
		for i, ets_map1 in enumerate(self.ecotype_maps):
			num_ets1 = len(ets_map1)
			m_j = 0
			for j, ets_map2 in enumerate(self.ecotype_maps):
				num_ets2 = len(ets_map2)
				if i == j:
					var_sclice = gen_var_list[i] * self.K[ets_map1, :][:, ets_map2] \
							+ err_var_list[i] * sp.eye(num_ets1)
				else:
					if  her_list[i] > 0 and her_list[j] > 0:
						rho = self.corr_mat[i, j] / sp.sqrt(her_list[i] * her_list[j])
						if rho > 1.0: rho = 1.0
						if rho < -1.0: rho = -1.0
					else:
						rho = 0
					print i, j, rho
					var_sclice = rho * math.sqrt(gen_var_list[i] * gen_var_list[j]) * \
							self.K[ets_map1, :][:, ets_map2]
				V[m_i:m_i + num_ets1, m_j:m_j + num_ets2] = var_sclice
				m_j += num_ets2
			m_i += num_ets1
		H_sqrt = cholesky(V)
		H_sqrt_inv = sp.mat(H_sqrt.T).I
		return {'V':V, 'H_sqrt':H_sqrt, 'H_sqrt_inv':H_sqrt_inv}




	def get_phen_corr(self, phenotype_list, ecotype_list):
		print 'Calculating phenotype correlation'
		num_traits = len(phenotype_list)
		corr_mat = sp.ones((num_traits, num_traits))
		for i, (pvs1, ets1) in enumerate(zip(phenotype_list, ecotype_list)):
			for j in range(i):
				ets2 = ecotype_list[j]
				pvs2 = phenotype_list[j]
				common_ets = set(ets1).intersection(set(ets2))
				ets_ix1 = map(ets1.index, common_ets)
				ets_ix2 = map(ets2.index, common_ets)
				vs1 = [pvs1[et_i] for et_i in ets_ix1]
				vs2 = [pvs2[et_i] for et_i in ets_ix2]
				corr_mat[i, j] = sp.corrcoef(vs1, vs2)[0, 1]
				corr_mat[j, i] = corr_mat[i, j]
		print corr_mat
		return corr_mat



	def emmax_f_test_fm(self, snps, with_betas=False, type='common'):
		"""
		EMMAX Full model f-test 
		returns three sets of p-values
		
		"""
		s1 = time.time()
		print 'Getting variance estimates'
		res = self.get_variance_matrix()
		print 'Done.'

		X = sp.hstack([self.X, self.E])
		print 'Perfoming GLS scan'
		r = self._emmax_f_test_(snps, X, res['H_sqrt_inv'], E=self.E, with_betas=with_betas)
		secs = time.time() - s1
		if secs > 60:
			mins = int(secs) / 60
			secs = secs - mins * 60
			print 'Took %d mins and %f seconds to run EMMAX on the traits...' % (mins, secs)
		else:
			print 'Took %f seconds to run EMMAX on the traits..' % (secs)
		return r



	def emmax_f_test(self, snps, with_betas=False, type='common'):
		"""
		EMMAX implementation
		Single SNPs, multiple traits.
		
		"""
		s1 = time.time()
		print 'Getting variance estimates'
		res = self.get_variance_matrix()
		print 'Done.'

		X = sp.hstack([self.X, self.E])
		print 'Perfoming GLS scan'
		r = self._emmax_f_test_(snps, X, res['H_sqrt_inv'], with_betas=with_betas)
		secs = time.time() - s1
		if secs > 60:
			mins = int(secs) / 60
			secs = secs - mins * 60
			print 'Took %d mins and %f seconds to run EMMAX on the traits...' % (mins, secs)
		else:
			print 'Took %f seconds to run EMMAX on the traits..' % (secs)
		return r



	def _emmax_f_test_(self, snps, X, H_sqrt_inv, E=None, verbose=True, with_betas=False, return_transformed_snps=False):
		"""
		EMMAX implementation (in python)
		Single SNPs, multiple traits
		
		Methods:
			normal - Normal regression
			qr - Uses QR decomposition to speed up regression with many co-factors.
			
		If E (environment) is provided then GxE factors are considered as well.
			
		"""
		dtype = 'single'
		q = 1  # Single SNP is being tested
		p = len(X.T) + q
		n = self.n
		n_p = n - p
		num_snps = len(snps)

		h0_X = sp.mat(H_sqrt_inv * X, dtype=dtype)
		Y = H_sqrt_inv * self.Y	#The transformed outputs.
		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
		Y = sp.mat(Y - h0_X * h0_betas, dtype=dtype)
		h0_betas = map(float, list(h0_betas))

		#H_sqrt_inv = H_sqrt_inv * self.Z
		(Q, R) = qr_decomp(h0_X)  #Do the QR-decomposition for the Gram-Schmidt process.
		Q = sp.mat(Q)
		Q2 = Q * Q.T
		M_0 = sp.mat(H_sqrt_inv.T * (sp.eye(n) - Q2), dtype=dtype)
		M = self.Z.T * M_0

		if E != None:
			q_e = q + E.shape[1]
			p_e = p + E.shape[1]
			n_p_e = n - p_e
			M_Es = []  #A list of matrices to genereate each GxE interaction
			for i in range(E.shape[1]):
				M_Es.append(sp.mat(self.Z * sp.array(E[:, i:i + 1]), dtype=dtype).T * M_0)
			rss_full_list = sp.repeat(h0_rss, num_snps)


		rss_list = sp.repeat(h0_rss, num_snps)
		chunk_size = len(Y)
		for i in range(0, num_snps, chunk_size): #Do the dot-product in chuncks!
			snps_chunk = sp.matrix(snps[i:i + chunk_size], dtype=dtype)
			Xs = snps_chunk * M
			if E != None:
				X_Es = []
				for M_E in M_Es:
					X_Es.append(snps_chunk * M_E)
			for j, X_j in enumerate(Xs):
				(betas, rss, p, sigma) = linalg.lstsq(X_j.T, Y)
				if rss:
					rss_list[i + j] = rss[0]
					if E != None:
						l = [X_E[j] for X_E in X_Es]
						l.append(X_j)
						X_mat = sp.mat(sp.vstack(l), dtype=dtype).T
						(betas, rss_full, p, sigma) = linalg.lstsq(X_mat, Y, overwrite_a=True)
						if rss_full and (rss - rss_full) / h0_rss > -0.0001:
							rss_full_list[i + j] = rss_full[0] if rss_full < rss else rss[0]
						elif rss_full and (rss - rss_full) / h0_rss <= -0.0001:
							print "Funny business"
							pdb.set_trace()

				if verbose and num_snps >= 10 and (i + j + 1) % (num_snps / 10) == 0: #Print dots
					sys.stdout.write('.')
					sys.stdout.flush()

		if verbose and num_snps >= 10:
			sys.stdout.write('\n')
		rss_ratio = h0_rss / rss_list
		var_perc = 1 - 1 / rss_ratio
		f_stats = (rss_ratio - 1) * n_p / float(q)
		p_vals = stats.f.sf(f_stats, q, n_p)

		res_d = {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'var_perc':var_perc,
			'h0_rss':h0_rss, 'h0_betas':h0_betas}

		if E != None:
			rss_ratio = h0_rss / rss_full_list
			res_d['var_perc_full'] = 1 - 1 / rss_ratio
			f_stats_full = (rss_ratio - 1) * n_p_e / float(q_e)
			res_d['ps_full'] = stats.f.sf(f_stats_full, q_e, n_p_e)
			res_d['f_stats_full'] = f_stats_full


			rss_ratio = rss_list / rss_full_list
			res_d['var_perc_ge'] = 1 - 1 / rss_ratio
			f_stats_ge = (rss_ratio - 1) * (n_p_e + q) / float(q_e - q)
			res_d['ps_ge'] = stats.f.sf(f_stats_ge, q_e - q, n_p_e + q)
			res_d['f_stats_ge'] = f_stats_ge
		return res_d



def _test_joint_analysis_(run_id='FT_suzi_log'):
	import dataParsers as dp
	import phenotypeData as pd
	import env
	#filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/seed_size.csv"
	filename = env.env['phen_dir'] + 'phen_raw_112210.csv'
	phed = pd.parse_phenotype_file(filename)
	phed.convert_to_averages()

	mac_threshold = 15
	#pids = [187, 188, 189, 190]
	pids = [5, 6, 7]
	#pids = [164, 165, 166]
	#pids = [1, 2, 3, 4]

	joint_phen_vals = []
	joint_ecotypes = []
	gen_var_list = []
	err_var_list = []
	her_list = []
	k_list = []
	for pid in pids:
		print phed.get_name(pid)
		sd = dp.load_250K_snps(debug_filter=0.05)
		sd.coordinate_w_phenotype_data(phed, pid)
		phed.normalize_values(pid)
		phenotypes = phed.get_values(pid)
		phed.log_transform(pid)
		#Create a new phenotype.
		joint_phen_vals.extend(phenotypes)
		ecotypes = phed.get_ecotypes(pid)
		#Figure out all ecotypes which are in 250K and for which we have phenotypic values.
		joint_ecotypes.extend(ecotypes)
		K = load_kinship_from_file(env.env['data_dir'] + 'kinship_matrix_cm72.pickled', ecotypes)
		#Get pseudoheritabilities for all phenotypes
		lmm = LinearMixedModel(phenotypes)
		lmm.add_random_effect(K)
		res = lmm.get_REML()
		print 'Pseudo-heritatbility:', res['pseudo_heritability']
		her_list.append(res['pseudo_heritability'])
		err_var_list.append(res['ve'])
		gen_var_list.append(res['vg'])


	#E in the model.
	num_vals = len(joint_phen_vals)
	joint_X = [[1] * num_vals]
	shift = phed.num_vals(pids[0])
	for pid in pids[1:]:
		n = phed.num_vals(pid)
		joint_X.append([0] * shift + [1] * n + [0] * (num_vals - n - shift))
		shift += n


	#Get correlations between traits..
	corr_mat, pids = phed.get_correlations(pids)
	#corr_mat = corr_mat * corr_mat
	print gen_var_list, err_var_list, corr_mat, pids

	unique_ecotypes = list(set(joint_ecotypes))
	sd = dp.load_250K_snps()
	sd.filter_accessions(unique_ecotypes)
	if mac_threshold:
		sd.filter_mac_snps(mac_threshold) #Filter MAF SNPs!
	#unique_ecotypes = sd.accessions
	K = load_kinship_from_file(env.env['data_dir'] + 'kinship_matrix_cm72.pickled', unique_ecotypes)

	#Create ecotype maps
	ecotype_maps = {}
	for pid in pids:
		l = []
		ets = phed.get_ecotypes(pid)
		for et in ets:
			l.append(unique_ecotypes.index(et))
		ecotype_maps[pid] = l

	#Construct new variance matrix
	V = sp.zeros((len(joint_ecotypes), len(joint_ecotypes)))
	m_i = 0
	for i, pid1 in enumerate(pids):
		ets_map1 = ecotype_maps[pid1]
		num_ets1 = len(ets_map1)
		m_j = 0
		for j, pid2 in enumerate(pids):
			ets_map2 = ecotype_maps[pid2]
			num_ets2 = len(ets_map2)
			if i == j:
				var_sclice = math.sqrt(gen_var_list[i] * gen_var_list[j]) * \
						K[ets_map1, :][:, ets_map2] + err_var_list[i] * sp.eye(num_ets1)
			else:
				if  her_list[i] > 0 and her_list[j] > 0:
					rho = corr_mat[i, j] / sp.sqrt(her_list[i] * her_list[j])
					if rho > 1.0: rho = 1.0
					if rho < -1.0: rho = 1.0
				else:
					rho = 0
				print i, j, rho
				var_sclice = rho * math.sqrt(gen_var_list[i] * gen_var_list[j]) * \
						K[ets_map1, :][:, ets_map2]
			V[m_i:m_i + num_ets1, m_j:m_j + num_ets2] = var_sclice
			m_j += num_ets2
		m_i += num_ets1

	print'Performing Cholesky decomposition'
	H_sqrt = cholesky(V)
	H_sqrt_inv = sp.mat(H_sqrt.T).I

	#pdb.set_trace()

	#5. Run analysis..
	lmm = LinearMixedModel(joint_phen_vals)
	lmm.set_factors(joint_X, False)
	#print lmm.X
	#Need to fix the SNPs!!!
	Z = sp.int16(sp.mat(joint_ecotypes).T == sp.mat(unique_ecotypes))
	snps = sd.getSnps()
 	positions = sd.getPositions()
 	chromosomes = sd.get_chr_list()
	#pdb.set_trace()
	res = lmm.emmax_full_model_gxe(snps, joint_X[1:], H_sqrt_inv, Z=Z)
	#pdb.set_trace()
	import gwaResults as gr
	res_h01 = gr.Result(scores=res['ps_h01'].tolist(), positions=positions, chromosomes=chromosomes)
	res_h02 = gr.Result(scores=res['ps_h02'].tolist(), positions=positions, chromosomes=chromosomes)
	res_h12 = gr.Result(scores=res['ps_h12'].tolist(), positions=positions, chromosomes=chromosomes)
	agr.qq_plot({'EMMAX_G':res_h01, 'EMMAX_FULL':res_h02, 'EMMAX_GxE':res_h12}, 1000, method_types=['emma', 'emma', 'emma'],
		mapping_labels=['EMMAX_G', 'EMMAX_FULL', 'EMMAX_GxE'], phenName='joint', pngFile=env.env['results_dir'] + 'qq_plot_%s.png' % run_id)
	agr.log_qq_plot({'EMMAX_G':res_h01, 'EMMAX_FULL':res_h02, 'EMMAX_GxE':res_h12}, 1000, 7, method_types=['emma', 'emma', 'emma'],
		mapping_labels=['EMMAX_G', 'EMMAX_FULL', 'EMMAX_GxE'], phenName='joint', pngFile=env.env['results_dir'] + 'log_qq_plot%s.png' % run_id)
	png_file_name = env.env['results_dir'] + 'manhattan_h01%s.png' % run_id
	res_h01.neg_log_trans()
	res_h01.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)
	png_file_name = env.env['results_dir'] + 'manhattan_h02%s.png' % run_id
	res_h02.neg_log_trans()
	res_h02.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)
	png_file_name = env.env['results_dir'] + 'manhattan_h12%s.png' % run_id
	res_h12.neg_log_trans()
	res_h12.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)




def get_emma_reml_estimates(y, K, cofactors=None, include_intercept=True):
	lmm = LinearMixedModel(y)
	lmm.add_random_effect(K)
	if cofactors != None:
		lmm.set_factors(cofactors, include_intercept=include_intercept)
	res = lmm.get_REML()
	res['Y_t'] = res['H_sqrt_inv'] * lmm.Y	#The transformed outputs.
	res['X_t'] = res['H_sqrt_inv'] * lmm.X
	res['lmm'] = lmm
	return res



def emma(snps, phenotypes, K, cofactors=None):
	"""
	Run EMMAX
	"""
	lmm = LinearMixedModel(phenotypes)
	lmm.add_random_effect(K)
	if cofactors:
		for cofactor in cofactors:
			lmm.add_factor(cofactor)

	print "Running EMMA (python)"
	s1 = time.time()
	res = lmm.expedited_REML_t_test(snps)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)
	return res



def emmax(snps, phenotypes, K, cofactors=None, Z=None, with_betas=False, progress_file_writer=None, emma_num=0):
	"""
	Run EMMAX
	"""
	lmm = LinearMixedModel(phenotypes)
	if Z != None:
		lmm.add_random_effect(Z * K * Z.T)
		if cofactors:
			for cofactor in cofactors:
				lmm.add_factor(Z * cofactor)
	else:
		lmm.add_random_effect(K)
		if cofactors:
			for cofactor in cofactors:
				lmm.add_factor(cofactor)

	if not progress_file_writer == None:
		progress_file_writer.update_progress_bar(0.4, 'Performing EMMAX')
	print "Running EMMAX"
	s1 = time.time()
	res = lmm.emmax_f_test(snps, Z=Z, with_betas=with_betas, progress_file_writer=progress_file_writer,
			emma_num=emma_num)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)
	return res


def emmax_perm_test(snps, phenotypes, K, num_perm=100):
	"""
	Run EMMAX
	"""
	lmm = LinearMixedModel(phenotypes)
	lmm.add_random_effect(K)
	print "Running %d EMMAX-permutation (writes %d dots)" % (num_perm, num_perm)
	s1 = time.time()
	res = lmm.emmax_permutations(snps, num_perm)
	p_f_list = zip(res['min_ps'], res['max_f_stats'])
	p_f_list.sort()
	print p_f_list[:10]
	threshold = p_f_list[len(p_f_list) / 20]
	res['threshold_05'] = threshold
	print 'Tresholds should be:', threshold
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)
	return res


def anova(snps, phenotypes):
	"""
	Run EMMAX
	"""
	lmm = LinearModel(phenotypes)

	print "Running ANOVA"
	s1 = time.time()
	res = lmm.anova_f_test(snps)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)
	return res



def emmax_anova(snps, phenotypes, K):
	"""
	Run EMMAX
	"""
	lmm = LinearMixedModel(phenotypes)
	lmm.add_random_effect(K)

	print "Running EMMAX-ANOVA"
	s1 = time.time()
	res = lmm.emmax_anova_f_test(snps)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)
	return res


#
#def _get_interactions_(isnp, snps, mac_threshold=15):
#	isnps = []
#	cofactor_indices = []
#	anti_isnp = sp.vectorize(lambda x: 1 if x == 0 else 0)(isnp)
#	for i, snp in enumerate(snps):
#		min_count = min(min(sp.bincount(isnp & snp)), min(sp.bincount(isnp | snp)), \
#			min(sp.bincount(snp & anti_isnp)), min(sp.bincount(snp | anti_isnp)))
#		if min_count > mac_threshold and min:
#			isnps.append(isnp & snp)
#			cofactor_indices.append(i)
#	return isnps, cofactor_indices



def _log_choose_(n, k):
 	#log_fact = lambda n: n * sp.log(n) - n + (sp.log(n * (1 + 4 * n * (1 + 2 * n))) / 6) + sp.log(sp.pi) / 2
 	#Srinivasa Ramanujan approximation of log(n!)
	if k == 0 or n == k:
		return 0
	if n < k:
		raise Exception('Out of range.')
	return sum(map(sp.log, range(n, n - k, -1))) - sum(map(sp.log, range(k, 0, -1)))


def _calc_bic_(ll, num_snps, num_par, n):
	bic = -2 * (ll) + num_par * sp.log(n)
	extended_bic = bic + \
		2 * _log_choose_(num_snps, num_par - 2)
	modified_bic = bic + \
		2 * (num_par) * sp.log(num_snps / 2.2 - 1)
	return (bic, extended_bic, modified_bic)


def _plot_manhattan_and_qq_(file_prefix, step_i, pvals, quantiles_dict=None, plot_bonferroni=True,
			highlight_markers=None, cand_genes=None, plot_xaxis=True, log_qq_max_val=5, with_qq_plots=True,
			num_dots=1000, simple_qq=False, highlight_loci=None, write_pvals=False,
			highlight_ppa_markers=None, ppas=None, **kwargs):
	import pylab
	import gwaResults as gr
	cm = pylab.get_cmap('hsv')

	png_file_name = '%s_step%d.png' % (file_prefix, step_i)
	res = gr.Result(scores=pvals, **kwargs)
	if write_pvals:
		pval_file_name = '%s_step%d.pvals' % (file_prefix, step_i)
		res.write_to_file(pval_file_name)


	res.neg_log_trans()
	res.plot_manhattan(png_file=png_file_name, plot_bonferroni=True, highlight_markers=highlight_markers,
				cand_genes=cand_genes, plot_xaxis=plot_xaxis, highlight_loci=highlight_loci)
	ret_dict = {'manhattan':png_file_name}


	if ppas != None:
		png_file_name = '%s_ppa_step%d.png' % (file_prefix, step_i)

		res = gr.Result(scores=ppas, **kwargs)
		res.plot_manhattan(png_file=png_file_name, highlight_markers=highlight_ppa_markers,
					cand_genes=cand_genes, plot_xaxis=plot_xaxis,
					highlight_loci=highlight_loci, ylab='Post. prob. of assoc.',
					min_score=0.0, max_score=1.0, markersize=3)
		ret_dict['ppa_manhattan'] = png_file_name
		if write_pvals:
			pval_file_name = '%s_step%d.ppas' % (file_prefix, step_i)
			res.write_to_file(pval_file_name)

	#Plot QQ plot
	if with_qq_plots and quantiles_dict != None:
		#calculate qq-plot line..
		log_quantiles = agr.get_log_quantiles(pvals, num_dots=num_dots, max_val=log_qq_max_val)
		quantiles = agr.get_quantiles(pvals, num_dots=num_dots)
		quantiles_dict['log'].append(log_quantiles)
		quantiles_dict['norm'].append(quantiles)
		quantiles_dict['labels'].append('Step %d' % (step_i))
		#plot all lines
		log_png_file_name = '%s_step%d_log_qqplot.png' % (file_prefix, step_i)
		png_file_name = '%s_step%d_qqplot.png' % (file_prefix, step_i)
		if step_i == 0:
			color_list = [cm(0.0)]
		else:
			color_list = [cm(i / float(step_i) * 0.7) for i in range(step_i + 1)]
		if simple_qq:
			log_qs = [quantiles_dict['log'][0], log_quantiles]
			qs = [quantiles_dict['norm'][0], quantiles]
			q_labs = [quantiles_dict['labels'][0], quantiles_dict['labels'][-1]]
			lcols = [cm(0), cm(0.7)]
			agr.simple_log_qqplot(log_qs, log_png_file_name, q_labs, line_colors=lcols,
					num_dots=num_dots, max_val=log_qq_max_val)
			agr.simple_qqplot(qs, png_file_name, q_labs, line_colors=color_list, num_dots=num_dots)

		elif step_i > 6:
			log_qs = quantiles_dict['log'][0:6] + [log_quantiles]
			qs = quantiles_dict['norm'][0:6] + [quantiles]
			q_labs = quantiles_dict['labels'][0:6] + [quantiles_dict['labels'][-1]]
			lcols = [cm(i / 6.0 * 0.7) for i in range(7)]
			agr.simple_log_qqplot(log_qs, log_png_file_name, q_labs, line_colors=lcols, num_dots=num_dots,
					max_val=log_qq_max_val)
			agr.simple_qqplot(qs, png_file_name, q_labs, line_colors=lcols, num_dots=num_dots)
		else:
			agr.simple_log_qqplot(quantiles_dict['log'], log_png_file_name, quantiles_dict['labels'],
					line_colors=color_list, num_dots=num_dots, max_val=log_qq_max_val)
			agr.simple_qqplot(quantiles_dict['norm'], png_file_name, quantiles_dict['labels'],
					line_colors=color_list, num_dots=num_dots)
		ret_dict['log_qq'] = log_png_file_name
		ret_dict['qq'] = png_file_name
	return ret_dict


def _analyze_opt_criterias_(criterias, sign_threshold, max_num_cofactors, file_prefix, with_qq_plots, lm, step_info_list,
			chr_pos_list, quantiles_dict, plot_bonferroni=True,
			cand_genes=None, plot_xaxis=True, log_qq_max_val=5, eig_L=None, type='emmax',
			highlight_loci=None, write_pvals=False, snp_priors=None, ppa_threshold=0.95, **kwargs):
	"""
	Copies or plots optimal criterias
	"""
	ret_dict = {}
	opt_indices = {}
	opt_file_dict = {}
	for c in criterias:
		print 'GWAs for optimal %s criteria:' % c
		if c == 'bonf':
			opt_list = sp.arange(max_num_cofactors + 1)
			for i, pval in enumerate(criterias['bonf']):
				if pval > sign_threshold:
					opt_list[i] = -1
			i_opt = opt_list.argmax()
		elif c == 'mbonf':
			fw_opt_list = sp.arange(max_num_cofactors + 1)
			for i in range(max_num_cofactors + 1):
				pval = criterias[c][i]
				if pval > sign_threshold:
					fw_opt_list[i] = -1
			fw_i_opt = fw_opt_list.argmax()
			fw_max = fw_opt_list[fw_i_opt]

			if max_num_cofactors > 1:
				shift = max_num_cofactors + 1
				bw_opt_list = sp.arange(max_num_cofactors - 1, 0, -1)
				for i in range(len(bw_opt_list)):
					pval = criterias[c][i + shift]
					if pval > sign_threshold:
						bw_opt_list[i] = -1
				bw_i_opt = bw_opt_list.argmax()
				bw_max = bw_opt_list[bw_i_opt]
				bw_i_opt = bw_opt_list.argmax() + shift
				if bw_max == fw_max:
					i_opt = bw_i_opt if criterias[c][fw_i_opt] > criterias[c][bw_i_opt] else fw_i_opt
				else:
					i_opt = bw_i_opt if bw_max > fw_max else fw_i_opt
			else:
				i_opt = fw_i_opt
		elif c == 'min_cof_ppa':
			fw_opt_list = sp.arange(max_num_cofactors + 1)
			for i in range(max_num_cofactors + 1):
				ppa = criterias[c][i]
				if ppa < ppa_threshold:
					fw_opt_list[i] = -1
			fw_i_opt = fw_opt_list.argmax()
			fw_max = fw_opt_list[fw_i_opt]

			if max_num_cofactors > 1:
				shift = max_num_cofactors + 1
				bw_opt_list = sp.arange(max_num_cofactors - 1, 0, -1)
				for i in range(len(bw_opt_list)):
					ppa = criterias[c][i + shift]
					if ppa < ppa_threshold:
						bw_opt_list[i] = -1
				bw_i_opt = bw_opt_list.argmax()
				bw_max = bw_opt_list[bw_i_opt]
				bw_i_opt = bw_opt_list.argmax() + shift
				if bw_max == fw_max:
					i_opt = bw_i_opt if criterias[c][fw_i_opt] > criterias[c][bw_i_opt] else fw_i_opt
				else:
					i_opt = bw_i_opt if bw_max > fw_max else fw_i_opt
			else:
				i_opt = fw_i_opt

		else:
			cur_min_val = criterias[c][0]
			min_indices = [0]
			for i in range(1, len(criterias[c])):
				v = criterias[c][i]
				if v < cur_min_val:
					cur_min_val = v
					min_indices = [i]
				if v == cur_min_val:
					min_indices.append(i)
			i_opt = min(min_indices)
			#i_opt = sp.argmin(criterias[c])
		print "    %d'th step was optimal." % i_opt
		ret_dict[c] = i_opt
		if i_opt <= max_num_cofactors:
			#Copy the pngs...
			if file_prefix:
				png_file_name = '%s_step%d.png' % (file_prefix, i_opt)
				opt_png_file_name = '%s_step%d_opt_%s.png' % (file_prefix, i_opt, c)
				if platform.system() == 'Linux' or platform.system() == 'Darwin':
					os.spawnlp(os.P_NOWAIT, 'cp', 'cp', png_file_name, opt_png_file_name)
					if snp_priors != None:
						png_file_name = '%s_ppa_step%d.png' % (file_prefix, i_opt)
						opt_png_file_name = '%s_ppa_step%d_opt_%s.png' % (file_prefix, i_opt, c)
						os.spawnlp(os.P_NOWAIT, 'cp', 'cp', png_file_name, opt_png_file_name)
					if with_qq_plots:
						qq_png_file_name = '%s_step%d_qqplot.png' % (file_prefix, i_opt)
						opt_qq_png_file_name = '%s_step%d_opt_%s_qqplot.png' % (file_prefix, i_opt, c)
						os.spawnlp(os.P_NOWAIT, 'cp', 'cp', qq_png_file_name, opt_qq_png_file_name)
						log_qq_png_file_name = '%s_step%d_log_qqplot.png' % (file_prefix, i_opt)
						opt_log_qq_png_file_name = '%s_step%d_opt_%s_log_qqplot.png' % (file_prefix, i_opt, c)
						os.spawnlp(os.P_NOWAIT, 'cp', 'cp', log_qq_png_file_name, opt_log_qq_png_file_name)
		elif i_opt in opt_file_dict:
			if file_prefix:
				png_file_name = opt_file_dict[i_opt]['manhattan']
				opt_png_file_name = '%s_step%d_opt_%s.png' % (file_prefix, i_opt, c)
				if platform.system() == 'Linux' or platform.system() == 'Darwin':
					os.spawnlp(os.P_NOWAIT, 'cp', 'cp', png_file_name, opt_png_file_name)
					if snp_priors != None:
						png_file_name = opt_file_dict[i_opt]['ppa_manhattan']
						opt_png_file_name = '%s_ppa_step%d_opt_%s.png' % (file_prefix, i_opt, c)
						os.spawnlp(os.P_NOWAIT, 'cp', 'cp', png_file_name, opt_png_file_name)

					if with_qq_plots:
						qq_png_file_name = opt_file_dict[i_opt]['qq']
						opt_qq_png_file_name = '%s_step%d_opt_%s_qqplot.png' % (file_prefix, i_opt, c)
						os.spawnlp(os.P_NOWAIT, 'cp', 'cp', qq_png_file_name, opt_qq_png_file_name)
						log_qq_png_file_name = opt_file_dict[i_opt]['log_qq']
						opt_log_qq_png_file_name = '%s_step%d_opt_%s_log_qqplot.png' % (file_prefix, i_opt, c)
						os.spawnlp(os.P_NOWAIT, 'cp', 'cp', log_qq_png_file_name, opt_log_qq_png_file_name)

		elif not i_opt in opt_indices:
			#Perfom GWAS witht he optimal cofactors
			cofactor_snps = step_info_list[i_opt]['cofactor_snps']
			cofactors = step_info_list[i_opt]['cofactors']
			print cofactors
			lm.set_factors(cofactor_snps)
			if type == 'emmax':
				eig_R = lm._get_eigen_R_(X=lm.X)
				reml_res = lm.get_REML(eig_L=eig_L, eig_R=eig_R)
				H_sqrt_inv = reml_res['H_sqrt_inv']
				l_res = lm._emmax_f_test_(kwargs['snps'], H_sqrt_inv, snp_priors=snp_priors)
				min_pval_i = l_res['ps'].argmin()
				mahalnobis_rss = l_res['rss'][min_pval_i]
				print 'Min Mahalanobis RSS:', mahalnobis_rss
			elif type == 'lm':
				l_res = lm.fast_f_test(kwargs['snps'])
				min_pval_i = l_res['ps'].argmin()

			min_pval = l_res['ps'][min_pval_i]
			min_pval_chr_pos = chr_pos_list[min_pval_i]
			print 'Min p-value:', min_pval
			l_pvals = l_res['ps'].tolist()
			if file_prefix:
				opt_file_prefix = '%s_opt_%s' % (file_prefix, c)
				if snp_priors:
					ppa_cofactors = step_info_list[i_opt]['ppa_cofactors']
					ppas = l_res['ppas'].tolist()
				else:
					ppas = None
					ppa_cofactors = None
				opt_file_dict[i_opt] = _plot_manhattan_and_qq_(opt_file_prefix, i_opt, l_pvals, quantiles_dict,
								plot_bonferroni=True, highlight_markers=cofactors,
								cand_genes=cand_genes, plot_xaxis=plot_xaxis,
								log_qq_max_val=log_qq_max_val, with_qq_plots=with_qq_plots,
								simple_qq=True, highlight_loci=highlight_loci,
								write_pvals=write_pvals, highlight_ppa_markers=ppa_cofactors,
								ppas=ppas, **kwargs)

			opt_indices[i_opt] = {'min_pval':min_pval, 'min_pval_chr_pos':min_pval_chr_pos,
						'kolmogorov_smirnov':agr.calc_ks_stats(l_pvals),
						'pval_median':agr.calc_median(l_pvals)}
			if type == 'emmax':
				opt_indices[i_opt]['mahalanobis_rss'] = mahalnobis_rss
	return ret_dict, opt_indices


def _cofactors_to_string_(cofactors):
	st = ''
	if len(cofactors) > 0:
		for tup in cofactors[:-1]:
			st += '%d_%d_%f,' % (tup[0], tup[1], tup[2])
		st += '%d_%d_%f' % (cofactors[-1][0], cofactors[-1][1], cofactors[-1][2])
	return st


def _plot_stepwise_stats_(file_prefix, step_info_list, sign_threshold, type='emmax'):
	import pylab
	rss_list = []
	ll_list = []
	bic_list = []
	e_bic_list = []
	m_bic_list = []
	mbonf_list = []
	min_pval_list = []
	ks_stat_list = []
	pval_median_list = []
	full_steps = []
	d_keys = [ 'rss', 'll', 'bic', 'e_bic', 'm_bic', 'min_pval', 'mbonf', 'pval_median']
	if type == 'emmax':
		p_her_list = []
		reml_mahalanobis_rss_list = []
		mahalanobis_rss_list = []
		d_keys += ['pseudo_heritability', 'reml_mahalanobis_rss', 'mahalanobis_rss']

	f = open(file_prefix + "_stats.csv", 'w')
	f.write(','.join(['step_nr'] + d_keys + ['kolmogorov_smirnov', 'min_pval_pos_chr', 'cofactors']) + '\n')
	for i, si in enumerate(step_info_list):
		st = ','.join(map(str, [i] + [si[k] for k in d_keys]))
		if si['kolmogorov_smirnov']:
			st += ',%d' % si['kolmogorov_smirnov']['D']
		else:
			st += ','
		if si['min_pval_chr_pos']:
			st += ',%d_%d,' % si['min_pval_chr_pos']
		else:
			st += ',,'
		st += _cofactors_to_string_(si['cofactors'])
		st += '\n'
		f.write(st)
		rss_list.append(float(si['rss']))
		if type == 'emmax':
			p_her_list.append(float(si['pseudo_heritability']))
			reml_mahalanobis_rss_list.append(float(si['reml_mahalanobis_rss']))
			if si['mahalanobis_rss']:
				mahalanobis_rss_list.append(float(si['mahalanobis_rss']))
		ll_list.append(si['ll'])
		bic_list.append(si['bic'])
		e_bic_list.append(si['e_bic'])
		m_bic_list.append(si['m_bic'])
		mbonf_list.append(si['mbonf'])
		if si['min_pval']:
			min_pval_list.append(float(si['min_pval']))
			ks_stat_list.append(si['kolmogorov_smirnov']['D'])
			pval_median_list.append(si['pval_median'])
			full_steps.append(i)
	f.close()
	pylab.figure(figsize=(6, 4))
	pylab.axes([0.15, 0.1, 0.83, 0.85])
	num_steps = len(step_info_list)
	rss_list.append(rss_list[0])
	y_max = max(rss_list)
	pylab.axis([-0.025 * num_steps, 1.025 * num_steps, 0, 1.025 * y_max])
	pylab.plot(range(len(rss_list)), rss_list, 'o-', alpha=0.7)
	pylab.ylabel('RSS')
	pylab.xlabel('Number of steps')
	pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
	pylab.savefig(file_prefix + '_stats_rss.png')
	pylab.clf()
	pylab.axes([0.15, 0.1, 0.83, 0.85])
	pylab.axis([-0.025 * num_steps / 2, 1.025 * num_steps / 2, 0, 1.025 * y_max])
	pylab.plot(range((num_steps / 2) + 1), rss_list[:(num_steps / 2) + 1], 'o-', alpha=0.7)
	pylab.ylabel('RSS')
	pylab.xlabel('Number of steps')
	pylab.savefig(file_prefix + '_stats_rss2.png')
	pylab.clf()
	ll_list.append(ll_list[0])
	pylab.plot(range(len(ll_list)), ll_list, 'o-', alpha=0.7)
	pylab.ylabel('Log likelihood')
	pylab.xlabel('Number of steps')
	pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
	pylab.savefig(file_prefix + '_stats_ll.png')
	pylab.clf()
	bic_list.append(bic_list[0])
	pylab.plot(range(len(bic_list)), bic_list, 'o-', alpha=0.7)
	pylab.ylabel('BIC')
	pylab.xlabel('Number of steps')
	pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
	pylab.savefig(file_prefix + '_stats_bic.png')
	pylab.clf()
	e_bic_list.append(e_bic_list[0])
	pylab.plot(range(len(e_bic_list)), e_bic_list, 'o-', alpha=0.7)
	pylab.ylabel('Extended BIC')
	pylab.xlabel('Number of steps')
	pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
	pylab.savefig(file_prefix + '_stats_ebic.png')
	pylab.clf()
	m_bic_list.append(m_bic_list[0])
	pylab.plot(range(len(m_bic_list)), m_bic_list, 'o-', alpha=0.7)
	pylab.ylabel('Modified BIC')
	pylab.xlabel('Number of steps')
	pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
	pylab.savefig(file_prefix + '_stats_mbic.png')
	pylab.clf()
	pylab.plot(full_steps, -sp.log10(min_pval_list), 'o-', alpha=0.7)
	pylab.axhline(y= -sp.log10(sign_threshold), c='k', linestyle='--')
	pylab.ylabel('Min. p-value')
	pylab.xlabel('Number of steps')
	pylab.savefig(file_prefix + '_stats_pval.png')
	pylab.clf()
	pylab.plot(full_steps, ks_stat_list, 'o-', alpha=0.7)
	pylab.ylabel('Kolmogorov-Smirnov statistic')
	pylab.xlabel('Number of steps')
	pylab.savefig(file_prefix + '_stats_ks.png')
	pylab.clf()
	pylab.plot(full_steps, pval_median_list, 'o-', alpha=0.7)
	pylab.axhline(y=0, c='k', linestyle='--')
	pylab.ylabel('0.5 - median p-value')
	pylab.xlabel('Number of steps')
	pylab.savefig(file_prefix + '_stats_mpval.png')
	pylab.clf()
	pylab.plot(range(1, len(mbonf_list)), -sp.log10(mbonf_list[1:]), 'o-', alpha=0.7)
	pylab.axhline(y= -sp.log10(sign_threshold), c='k', linestyle='--')
	pylab.ylabel('Max cofactor pvalue (-log10[pval])')
	pylab.xlabel('Number of steps')
	pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
	pylab.savefig(file_prefix + '_stats_mbonf.png')
	pylab.clf()

	#Plotting variance partition plots
	if type == 'emmax':
		p_her_list.append(p_her_list[0])
		pylab.plot(range(len(p_her_list)), p_her_list, 'o-', alpha=0.7)
		pylab.ylabel('Pseudo-heritability')
		pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
		pylab.savefig(file_prefix + '_stats_p_her.png')
		pylab.clf()
		reml_mahalanobis_rss_list.append(reml_mahalanobis_rss_list[0])
		pylab.plot(range(len(reml_mahalanobis_rss_list)), reml_mahalanobis_rss_list, 'o-', alpha=0.7)
		pylab.ylabel('REML Mahalanobis RSS')
		pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
		pylab.savefig(file_prefix + '_stats_reml_mahalanobis_rss.png')
		pylab.clf()
		pylab.plot(full_steps, mahalanobis_rss_list, 'o-', alpha=0.7)
		pylab.ylabel('Mahalanobis RSS')
		pylab.savefig(file_prefix + '_stats_mahalanobis_rss.png')
		pylab.clf()

		max_rss = max(rss_list)
		rss_array = sp.array(rss_list) / max_rss
		p_her_array = rss_array * sp.array(p_her_list)
		genetic_variance = p_her_array + (1 - rss_array)
		variance_explained = (1 - rss_array)
		pylab.figure(figsize=(12, 6))
		pylab.axes([0.05, 0.08, 0.94, 0.90])
		pylab.fill_between([0, num_steps], 0, 1, color='#DD3333', alpha=0.8, label='Variance explained')
		pylab.fill_between(sp.arange(num_steps + 1), 0, genetic_variance, color='#22CC44', alpha=0.8, label='Genetic variance')
		pylab.fill_between(sp.arange(num_steps + 1), 0, variance_explained, color='#2255AA', alpha=0.8, label='Variance explained')
		pylab.ylabel('Percentage of variance')
		pylab.xlabel('Step number')
		pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
		pylab.legend(loc=1, ncol=3, shadow=True)
		pylab.axis([0, num_steps, 0, 1])
		pylab.savefig(file_prefix + '_stats_variances.png', format='png')

		pylab.figure(figsize=(10, 6))
		pylab.axes([0.06, 0.08, 0.92, 0.90])
		num_steps = num_steps / 2
		pylab.fill_between([0, num_steps ], 0, 1, color='#DD3333', alpha=0.8, label='Variance explained')
		pylab.fill_between(sp.arange(num_steps + 1), 0, genetic_variance[:num_steps + 1], color='#22CC44', alpha=0.8, label='Genetic variance')
		pylab.fill_between(sp.arange(num_steps + 1), 0, variance_explained[:num_steps + 1], color='#2255AA', alpha=0.8, label='Variance explained')
		pylab.ylabel('Percentage of variance')
		pylab.xlabel('Step number')
		pylab.legend(loc=1, ncol=3, shadow=True)
		pylab.axis([0, num_steps, 0, 1])
		pylab.savefig(file_prefix + '_stats_variances_forward.png', format='png')



def emmax_step(phen_vals, sd, K, cof_chr_pos_list, eig_L=None, eig_R=None, progress_file_writer=None, plot_prefix=None):
	"""
	EMMAX single SNPs
	
	Returns various stats useful for stepwise regression if return_stepw_stats flag is set to True.
	"""
	import bisect
	import gwaResults as gr
	s1 = time.time()
       	lmm = LinearMixedModel(phen_vals)
	lmm.add_random_effect(K)
	h0_rss = lmm.get_rss()[0]
	chrom_pos_list = sd.getChrPosList()
	print 'Looking up cofactors'
	cof_indices = []
	for chrom_pos in cof_chr_pos_list:
		i = bisect.bisect(chrom_pos_list, chrom_pos) - 1
		assert chrom_pos_list[i] == chrom_pos, 'Cofactor missing??'
		cof_indices.append(i)
	snps = sd.getSnps()
	cof_snps = [snps[i] for i in cof_indices]
	lmm.set_factors(cof_snps)

	if not eig_L:
		print 'Calculating the eigenvalues of K'
		eig_L = lmm._get_eigen_L_()
		print 'Done.'
	if not eig_R:
		print "Calculating the eigenvalues of S(K+I)S where S = I-X(X'X)^-1X'"
		eig_R = lmm._get_eigen_R_()
		print 'Done'


	print 'Getting variance estimates'
	res_dict = {'REML':lmm.get_estimates(eig_L, method='REML', eig_R=eig_R),
			'ML':lmm.get_estimates(eig_L, method='ML', eig_R=eig_R)}
	print 'Done.'
	print 'pseudo_heritability:', res_dict['REML']['pseudo_heritability']

	H_sqrt_inv = res_dict['REML']['H_sqrt_inv']

	if not progress_file_writer == None:
		progress_file_writer.update_progress_bar(progress=0.45, task_status='Performing EMMAX')
		progress_file_writer.set_step(0.05)


	r = lmm._emmax_f_test_(snps, H_sqrt_inv, progress_file_writer=progress_file_writer)
	min_pval_i = sp.argmin(r['ps'])
	step_dict = {}
	step_dict['min_pval_i'] = min_pval_i
	step_dict['min_pval'] = r['ps'][min_pval_i]
	step_dict['mahalanobis_rss'] = r['rss'][min_pval_i]
	step_dict['min_pval_chr_pos'] = chrom_pos_list[min_pval_i]
	step_dict['pseudo_heritability'] = res_dict['REML']['pseudo_heritability']
	step_dict['ve'] = res_dict['REML']['ve']
	step_dict['vg'] = res_dict['REML']['vg']
	ll = res_dict['ML']['max_ll']
	step_dict['max_ll'] = ll
	step_dict['rss'] = float(res_dict['REML']['rss'])
	step_dict['h0_rss'] = h0_rss
	step_dict['perc_var_expl'] = min(0.0, 1.0 - step_dict['rss'] / h0_rss)
	step_dict['reml_mahalanobis_rss'] = float(res_dict['REML']['mahalanobis_rss'])
	num_snps = len(snps)
	num_par = lmm.X.shape[1] + 1
	step_dict['num_snps'] = num_snps
	step_dict['num_par'] = num_par
	(bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lmm.n) #Calculate the BICs
	step_dict['ebic'] = extended_bic
	step_dict['mbic'] = modified_bic
	step_dict['bic'] = bic

	#Calculate maximum cofactor p-value
	cof_pvals = []
	cof_chrom_pos_pval_list = []
	for i, snp in enumerate(cof_snps):
		t_cofactors = sp.matrix(cof_snps) * H_sqrt_inv.T
		t_cofactors = t_cofactors.tolist()
		del t_cofactors[i]
		lmm.set_factors(t_cofactors)
		res = lmm._emmax_f_test_([snp], H_sqrt_inv)
		cof_pval = res['ps'][0]
		cof_pvals.append(cof_pval)
		cof_chrom_pos_pval_list.append((cof_chr_pos_list[i][0], cof_chr_pos_list[i][1], -math.log10(cof_pval)))
	for i, pval in zip(cof_indices, cof_pvals):
		print r['ps'][i], pval
		r['ps'][i] = pval
	if len(cof_pvals):
		step_dict['max_cof_pval'] = max(cof_pvals)
	else:
		step_dict['max_cof_pval'] = 0.0
	#step_dict['cofactor_snps'] = cof_snps
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)

	if plot_prefix != None:
		res = gr.Result(scores=list(r['ps']), snps_data=sd)
		png_file_name = '%s_%s_manhattan.png' % (plot_prefix, _cofactors_to_string_(cof_chrom_pos_pval_list))
		res.neg_log_trans()
		res.filter_percentile(0.1)
		res.plot_manhattan(png_file=png_file_name, plot_bonferroni=True, highlight_markers=cof_chrom_pos_pval_list)

	return {'stats':step_dict, 'res':r}




def emmax_step_wise(phenotypes, K, sd=None, num_steps=10, file_prefix=None, allow_interactions=False,
		interaction_pval_thres=0.01, forward_backwards=True, local=False, cand_gene_list=None,
		plot_xaxis=True, with_qq_plots=True, sign_threshold=None, log_qq_max_val=5,
		highlight_loci=None, save_pvals=False, snp_priors=None, K2=None, snp_choose_criteria='pval',
		emma_num=0, **kwargs):
	"""
	Run step-wise EMMAX forward-backward.
	"""
	import gwaResults as gr
	if local:
		with_qq_plots = False


	if sd:
	 	kwargs['snps'] = sd.getSnps()
	 	kwargs['positions'] = sd.getPositions()
	 	kwargs['chromosomes'] = sd.get_chr_list()
	 	d = sd.get_mafs()
	 	kwargs['macs'] = d['mafs']
	 	kwargs['mafs'] = d['marfs']
 	if snp_priors:
 		kwargs['snp_priors'] = snp_priors[:]

	snps = kwargs['snps'][:]
	positions = kwargs['positions'][:]
	chromosomes = kwargs['chromosomes'][:]
	mafs = kwargs['mafs'][:]
	macs = kwargs['macs'][:]

	chr_pos_list = zip(chromosomes, positions)
       	lmm = LinearMixedModel(phenotypes)
	lmm.add_random_effect(K)
	if K2 != None:
		lmm.add_random_effect(K2)
	num_snps = len(snps)

	if snp_priors == None:
		snp_priors = [1.0 / num_snps] * num_snps

	if not sign_threshold: #Then use Bonferroni threshold
		sign_threshold = 1.0 / (num_snps * 20.0)

	print "Running EMMAX stepwise"
	s1 = time.time()
 	step_info_list = []
	cofactors = []  #A list of the loci found, together with their statistics.
	cofactor_snps = []
	interactions = []
 	step_i = 0
 	num_par = 2 #mean and variance scalar
 	num_pher_0 = 0

	if K2 != None: #Then first estimate K
		res = lmm.get_estimates_3()
		K = res['opt_k']
	eig_L = lmm._get_eigen_L_(K=K)
	eig_R = lmm._get_eigen_R_(K=K)

	reml_res = lmm.get_REML(eig_L=eig_L, eig_R=eig_R)
	ml_res = lmm.get_ML(eig_L=eig_L, eig_R=eig_R)
	H_sqrt_inv = reml_res['H_sqrt_inv']
	ll = ml_res['max_ll']
	rss = float(reml_res['rss'])
	reml_mahalanobis_rss = float(reml_res['mahalanobis_rss'])
	criterias = {'ebics':[], 'mbics':[], 'bonf':[], 'mbonf':[]}
	(bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lmm.n) #Calculate the BICs
	criterias['ebics'].append(extended_bic)
	criterias['mbics'].append(modified_bic)
	max_cofactor_pval = 0
	criterias['mbonf'].append(max_cofactor_pval)
	criterias['bonf'].append(0)

	criterias['min_cof_ppa'] = [1] #Posterior probability of association
	cof_snp_priors = []
	ppa_cofactors = []

	action = 'None'
	print '\nStep %d: action=%s, num_par=%d, p_her=%0.4f, ll=%0.2f, rss=%0.2f, reml_m_rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f, num_snps=%d' % \
		(step_i, action, num_par, reml_res['pseudo_heritability'], ll, rss, reml_mahalanobis_rss, \
		 bic, extended_bic, modified_bic, num_snps)
	print 'Cofactors:', _cofactors_to_string_(cofactors)
	quantiles_dict = {'log':[], 'norm':[], 'labels':[]}

	for step_i in range(1, num_steps + 1):
		emmax_res = lmm._emmax_f_test_(snps, H_sqrt_inv, snp_priors=snp_priors, emma_num=emma_num)
		if step_i == 1:
			first_emmax_res = emmax_res
		min_pval_i = sp.argmin(emmax_res['ps'])
		min_pval = emmax_res['ps'][min_pval_i]
		mahalnobis_rss = emmax_res['rss'][min_pval_i]
		min_pval_chr_pos = chr_pos_list[min_pval_i]
		max_ppa_i = sp.argmax(emmax_res['ppas'])
		max_ppa = emmax_res['ppas'][max_ppa_i]
		print 'Min p-value:', min_pval
		criterias['bonf'].append(min_pval)
		print 'Min Mahalanobis RSS:', mahalnobis_rss
		step_info = {'pseudo_heritability':reml_res['pseudo_heritability'], 'rss':rss, \
			'reml_mahalanobis_rss': reml_res['mahalanobis_rss'], 'mahalanobis_rss':mahalnobis_rss,
			'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic, 'mbonf':max_cofactor_pval,
			'cofactors':map(tuple, cofactors[:]), 'cofactor_snps':cofactor_snps[:],
			'min_pval':min_pval, 'min_pval_chr_pos': min_pval_chr_pos, 'interactions':interactions,
			'max_ppa':max_ppa, 'max_ppa_pval':emmax_res['ps'][max_ppa_i],
			'max_ppa_chr_pos':chr_pos_list[max_ppa_i], 'ppa_cofactors':map(tuple, ppa_cofactors[:])}
		ppas = emmax_res['ppas'].tolist()

		if snp_choose_criteria == 'pval':
			snp_i = min_pval_i
		elif snp_choose_criteria == 'ppas':
			snp_i = max_ppa_i


		ex_pvals = emmax_res['ps'].tolist()
		if save_pvals:
			step_info['ps'] = ex_pvals

		#Plot gwas results per step 
		if file_prefix:
			_plot_manhattan_and_qq_(file_prefix, step_i - 1, ex_pvals,
					quantiles_dict, positions=positions, chromosomes=chromosomes, mafs=mafs, macs=macs, plot_bonferroni=True,
					highlight_markers=cofactors, cand_genes=cand_gene_list, plot_xaxis=plot_xaxis, log_qq_max_val=log_qq_max_val,
					with_qq_plots=with_qq_plots, highlight_loci=highlight_loci, write_pvals=save_pvals,
					ppas=ppas, highlight_ppa_markers=ppa_cofactors)



		step_info['kolmogorov_smirnov'] = agr.calc_ks_stats(ex_pvals)
		step_info['pval_median'] = agr.calc_median(ex_pvals)
		print step_info['kolmogorov_smirnov'], step_info['pval_median']
		step_info_list.append(step_info)


		#Adding the new SNP as a cofactor
		lmm.add_factor(snps[snp_i])
		cofactor_snps.append(snps[snp_i])

		if K2 != None: #Again first estimate K
			res = lmm.get_estimates_3()
			K = res['opt_k']
			eig_L = lmm._get_eigen_L_(K=K)


		eig_R = lmm._get_eigen_R_(X=lmm.X, K=K)
		reml_res = lmm.get_REML(eig_L=eig_L, eig_R=eig_R)
		ml_res = lmm.get_ML(eig_L=eig_L, eig_R=eig_R)
		H_sqrt_inv = reml_res['H_sqrt_inv']
		ll = ml_res['max_ll']
		rss = float(reml_res['rss'])
		reml_mahalanobis_rss = float(reml_res['mahalanobis_rss'])
		num_par += 1
		action = '+'
		cof_snp_priors.append(snp_priors[snp_i])
		ppa_cofactors.append([chromosomes[snp_i], positions[snp_i], max_ppa])
		cofactors.append([chromosomes[snp_i], positions[snp_i], min_pval])


		#Re-estimate the p-value of the cofactors... with the smallest in the list.
		cofactor_pvals = []
		if snp_priors != None:
			cofactor_ppas = [] #Posterior probabilities of association
		for i, snp in enumerate(cofactor_snps):
			t_cofactors = cofactor_snps[:]
			del t_cofactors[i]
			lmm.set_factors(t_cofactors)
			res = lmm._emmax_f_test_([snp], H_sqrt_inv, snp_priors=[cof_snp_priors[i]], emma_num=0)
			cofactor_ppas.append(res['ppas'][0])
			ppa_cofactors[i][2] = res['ppas'][0]
			pval = res['ps'][0]
			cofactor_pvals.append(pval)
			cofactors[i][2] = -math.log10(pval)
		lmm.set_factors(cofactor_snps)
		max_cofactor_pval = max(cofactor_pvals)
		criterias['mbonf'].append(max_cofactor_pval)
		if snp_priors != None:
			criterias['min_cof_ppa'].append(min(cofactor_ppas))


		#Remove the found SNP from considered SNPs
		del snps[snp_i]
		del snp_priors[snp_i]
		del positions[snp_i]
		del chromosomes[snp_i]
		del chr_pos_list[snp_i]
		del mafs[snp_i]
		del macs[snp_i]
		num_snps -= 1


		(bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lmm.n) #Calculate the BICs
		criterias['ebics'].append(extended_bic)
		criterias['mbics'].append(modified_bic)

		print '\nStep %d: action=%s, num_par=%d, p_her=%0.4f, ll=%0.2f, rss=%0.2f, reml_m_rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f, num_snps=%d' % \
			(step_i, action, num_par, reml_res['pseudo_heritability'], ll, rss, reml_mahalanobis_rss, \
			 bic, extended_bic, modified_bic, num_snps)

		print 'Cofactors:', _cofactors_to_string_(cofactors)
		print ppa_cofactors
#		if reml_res['pseudo_heritability'] < 0.01 and num_pher_0 < 1:
#			num_pher_0 += 1
#		elif reml_res['pseudo_heritability'] < 0.01:
		if reml_res['pseudo_heritability'] < 0.01:
			print 'Breaking early, since pseudoheritability is close to 0.'
			break

	emmax_res = lmm._emmax_f_test_(snps, H_sqrt_inv, snp_priors=snp_priors, emma_num=emma_num) #FINISH!!!
	min_pval_i = sp.argmin(emmax_res['ps'])
	min_pval = emmax_res['ps'][min_pval_i]
	mahalnobis_rss = emmax_res['rss'][min_pval_i]
	min_pval_chr_pos = chr_pos_list[min_pval_i]
	max_ppa_i = sp.argmax(emmax_res['ppas'])
	ppas = emmax_res['ppas'].tolist()
	print 'Min p-value:', min_pval
	print 'Min Mahalanobis RSS:', mahalnobis_rss
	step_info = {'pseudo_heritability':reml_res['pseudo_heritability'], 'rss':rss, \
		'reml_mahalanobis_rss': reml_res['mahalanobis_rss'], 'mahalanobis_rss':mahalnobis_rss,
		'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic, 'mbonf':max_cofactor_pval,
		'cofactors':map(tuple, cofactors[:]), 'cofactor_snps':cofactor_snps[:],
		'min_pval':min_pval, 'min_pval_chr_pos': min_pval_chr_pos, 'interactions':interactions,
		'max_ppa':emmax_res['ppas'][max_ppa_i], 'max_ppa_pval':emmax_res['ps'][max_ppa_i],
		'max_ppa_chr_pos':chr_pos_list[max_ppa_i], 'ppa_cofactors':map(tuple, ppa_cofactors[:])}


	ex_pvals = emmax_res['ps'].tolist()
	if save_pvals:
		step_info['ps'] = ex_pvals
	#Now plotting!
	print "Generating plots"
	if file_prefix:
		_plot_manhattan_and_qq_(file_prefix, step_i, ex_pvals, quantiles_dict, positions=positions,
					chromosomes=chromosomes, mafs=mafs, macs=macs, plot_bonferroni=True,
					highlight_markers=cofactors, cand_genes=cand_gene_list, plot_xaxis=plot_xaxis,
					log_qq_max_val=log_qq_max_val, with_qq_plots=with_qq_plots,
					highlight_loci=highlight_loci, write_pvals=save_pvals, ppas=ppas,
					highlight_ppa_markers=ppa_cofactors)
		#Plot posterior probabilities of association


	step_info['kolmogorov_smirnov'] = agr.calc_ks_stats(ex_pvals)
	step_info['pval_median'] = agr.calc_median(ex_pvals)
	print step_info['kolmogorov_smirnov'], step_info['pval_median']
	step_info_list.append(step_info)

	max_num_cofactors = len(cofactors)

	#Now backward stepwise.
	if forward_backwards:
		print 'Starting backwards..'
		while len(cofactor_snps) > 1:
			step_i += 1
			f_stats = sp.zeros(len(cofactor_snps))
			ppas = sp.zeros(len(cofactor_snps))
			for i, snp in enumerate(cofactor_snps):
				t_cofactors = cofactor_snps[:]
				del t_cofactors[i]
				lmm.set_factors(t_cofactors)
				res = lmm._emmax_f_test_([snp], H_sqrt_inv, snp_priors=[cof_snp_priors[i]], emma_num=0)
				ppas[i] = res['ppas'][0]
				cofactors[i][2] = -math.log10(res['ps'][0])
				f_stats[i] = res['f_stats'][0]
			if snp_choose_criteria == 'pval':
				i_to_remove = f_stats.argmin()
			elif snp_choose_criteria == 'ppas':
				i_to_remove = ppas.argmin()
			del ppa_cofactors[i_to_remove]
			del cofactor_snps[i_to_remove]
			del cofactors[i_to_remove]
			lmm.set_factors(cofactor_snps)
			num_snps += 1


			#Re-estimating the REML and ML.
			if K2 != None: #Again first estimate K
				res = lmm.get_estimates_3()
				K = res['opt_k']
				eig_L = lmm._get_eigen_L_(K=K)
			eig_R = lmm._get_eigen_R_(X=lmm.X, K=K)
			reml_res = lmm.get_REML(eig_L=eig_L, eig_R=eig_R)
			ml_res = lmm.get_ML(eig_L=eig_L, eig_R=eig_R)
			ll = ml_res['max_ll']
			H_sqrt_inv = reml_res['H_sqrt_inv']
			rss = float(reml_res['rss'])
			reml_mahalanobis_rss = float(reml_res['mahalanobis_rss'])
			num_par -= 1
			action = '-'


			#Update the p-values
			cofactor_pvals = []
			if snp_priors != None:
				cofactor_ppas = [] #Posterior probabilities of association
			for i, snp in enumerate(cofactor_snps):
				t_cofactors = cofactor_snps[:]
				del t_cofactors[i]
				lmm.set_factors(t_cofactors)
				res = lmm._emmax_f_test_([snp], H_sqrt_inv, snp_priors=[cof_snp_priors[i]], emma_num=0)
				cofactor_ppas.append(res['ppas'][0])
				ppa_cofactors[i][2] = res['ppas'][0]
				pval = res['ps'][0]
				cofactor_pvals.append(pval)
				cofactors[i][2] = -math.log10(pval)
			max_cofactor_pval = max(cofactor_pvals)
			criterias['mbonf'].append(max_cofactor_pval)
			criterias['min_cof_ppa'].append(min(cofactor_ppas))

			#Calculate the BICs
			(bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lmm.n)
			criterias['ebics'].append(extended_bic)
			criterias['mbics'].append(modified_bic)
			print '\nStep %d: action=%s, num_par=%d, p_her=%0.4f, ll=%0.2f, rss=%0.2f, reml_m_rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f, num_snps=%d' % \
				(step_i, action, num_par, reml_res['pseudo_heritability'], ll, rss,
				reml_mahalanobis_rss, bic, extended_bic, modified_bic, num_snps)

			print 'Cofactors:', _cofactors_to_string_(cofactors)
			print ppa_cofactors

			step_info = {'pseudo_heritability':reml_res['pseudo_heritability'], 'rss':rss, \
				'reml_mahalanobis_rss': reml_res['mahalanobis_rss'], 'll':ll, 'bic':bic,
				'e_bic':extended_bic, 'm_bic':modified_bic, 'mbonf':max_cofactor_pval,
				'cofactors':map(tuple, cofactors[:]), 'cofactor_snps':cofactor_snps[:],
				'mahalanobis_rss':None, 'min_pval':None, 'min_pval_chr_pos':None,
				'kolmogorov_smirnov':None, 'pval_median':None,
				'ppa_cofactors': map(tuple, ppa_cofactors[:])}
			step_info_list.append(step_info)
			print step_info['kolmogorov_smirnov'], step_info['pval_median']

	opt_dict, opt_indices = _analyze_opt_criterias_(criterias, sign_threshold, max_num_cofactors, file_prefix, with_qq_plots, lmm,
				step_info_list, chr_pos_list, quantiles_dict, plot_bonferroni=True, cand_genes=cand_gene_list,
				plot_xaxis=plot_xaxis, log_qq_max_val=log_qq_max_val, eig_L=eig_L, type='emmax',
				highlight_loci=highlight_loci, write_pvals=save_pvals, ** kwargs)

	for step_i in opt_indices:
		for h in ['mahalanobis_rss', 'min_pval', 'min_pval_chr_pos', 'kolmogorov_smirnov', 'pval_median']:
			step_info_list[step_i][h] = opt_indices[step_i][h]

	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)

	if file_prefix:
		_plot_stepwise_stats_(file_prefix, step_info_list, sign_threshold, type='emmax')

	res_dict = {'step_info_list':step_info_list, 'first_emmax_res':first_emmax_res, 'opt_dict':opt_dict}

	return res_dict





def lm_step_wise(phenotypes, sd=None, num_steps=10,
		file_prefix=None, allow_interactions=False,
		interaction_pval_thres=0.01, forward_backwards=True, local=False, cand_gene_list=None,
		plot_xaxis=True, with_qq_plots=True, sign_threshold=None,
		log_qq_max_val=5, highlight_loci=None, save_pvals=False, **kwargs):
	"""
	Run simple step-wise linear model forward-backward.
	"""
	#import plotResults as pr
	if local:
		with_qq_plots = False

	if sd:
	 	kwargs['snps'] = sd.getSnps()
	 	kwargs['positions'] = sd.getPositions()
	 	kwargs['chromosomes'] = sd.get_chr_list()
	 	d = sd.get_mafs()
	 	kwargs['macs'] = d['mafs']
	 	kwargs['mafs'] = d['marfs']

	snps = kwargs['snps'][:]
	positions = kwargs['positions'][:]
	chromosomes = kwargs['chromosomes'][:]
	mafs = kwargs['mafs'][:]
	macs = kwargs['macs'][:]

	chr_pos_list = zip(chromosomes, positions)
       	lm = LinearModel(phenotypes)
       	num_snps = len(snps)

	if not sign_threshold: #Then use Bonferroni threshold
		sign_threshold = 1.0 / (num_snps * 20.0)

	print "Running step-wise LM"
	s1 = time.time()
 	step_info_list = []
	cofactors = []  #A list of the loci found, together with their statistics.
	cofactor_snps = []
	interactions = []
 	step_i = 0
 	num_par = 2 #mean and variance scalar

	rss = lm.get_rss()
	ll = lm.get_ll(rss)
	criterias = {'ebics':[], 'mbics':[], 'bonf':[], 'mbonf':[]}
	(bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lm.n) #Calculate the BICs
	criterias['ebics'].append(extended_bic)
	criterias['mbics'].append(modified_bic)
	max_cofactor_pval = 0
	criterias['mbonf'].append(max_cofactor_pval)
	criterias['bonf'].append(0)
	action = 'None'
	print '\nStep %d: action=%s, num_par=%d, ll=%0.2f, rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f' % \
		(step_i, action, num_par, ll, rss, bic, extended_bic, modified_bic)
	print 'Cofactors:', _cofactors_to_string_(cofactors)
	quantiles_dict = {'log':[], 'norm':[], 'labels':[]}

	for step_i in range(1, num_steps + 1):
		lm_res = lm.fast_f_test(snps)
		if step_i == 1:
			first_lm_res = lm_res
		min_pval_i = sp.argmin(lm_res['ps'])
		min_pval = lm_res['ps'][min_pval_i]
		min_pval_chr_pos = chr_pos_list[min_pval_i]
		print 'Min p-value:', min_pval
		criterias['bonf'].append(min_pval)
		step_info = {'rss':rss, 'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic,
				'mbonf':max_cofactor_pval, 'cofactors':map(tuple, cofactors[:]),
				'cofactor_snps':cofactor_snps[:], 'min_pval':min_pval, 'min_pval_chr_pos': min_pval_chr_pos,
				 'interactions':interactions}

		lm_pvals = lm_res['ps'].tolist()
		#Plot gwas results per step 
		if file_prefix:
			_plot_manhattan_and_qq_(file_prefix, step_i - 1, lm_pvals, quantiles_dict, positions=positions,
					chromosomes=chromosomes, mafs=mafs, macs=macs, plot_bonferroni=True, highlight_markers=cofactors,
					cand_genes=cand_gene_list, plot_xaxis=plot_xaxis, log_qq_max_val=log_qq_max_val,
					with_qq_plots=with_qq_plots, highlight_loci=highlight_loci, write_pvals=save_pvals)
		if save_pvals:
			step_info['ps'] = lm_pvals


		if cand_gene_list:
			#Calculate candidate gene enrichments.
			pass
		step_info['kolmogorov_smirnov'] = agr.calc_ks_stats(lm_pvals)
		step_info['pval_median'] = agr.calc_median(lm_pvals)
		print step_info['kolmogorov_smirnov'], step_info['pval_median']
		step_info_list.append(step_info)


		#Adding the new SNP as a cofactor
		lm.add_factor(snps[min_pval_i])
		cofactor_snps.append(snps[min_pval_i])
		rss = lm.get_rss()
		ll = lm.get_ll(rss)
		num_par += 1
		action = '+'

		cofactors.append([min_pval_chr_pos[0], min_pval_chr_pos[1], min_pval])


		#Re-estimate the p-value of the cofactors... with the smallest in the list.
		cofactor_pvals = []
		for i, snp in enumerate(cofactor_snps):
			t_cofactors = cofactor_snps[:]
			del t_cofactors[i]
			lm.set_factors(t_cofactors)
			pval = lm.fast_f_test([snp])['ps'][0]
			cofactor_pvals.append(pval)
			cofactors[i][2] = pval
		lm.set_factors(cofactor_snps)
		max_cofactor_pval = max(cofactor_pvals)
		criterias['mbonf'].append(max_cofactor_pval)

		#Remove the found SNP from considered SNPs
		del snps[min_pval_i]
		del positions[min_pval_i]
		del chromosomes[min_pval_i]
		del chr_pos_list[min_pval_i]
		del mafs[min_pval_i]
		del macs[min_pval_i]
		num_snps -= 1


		(bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lm.n) #Calculate the BICs
		criterias['ebics'].append(extended_bic)
		criterias['mbics'].append(modified_bic)

		print '\nStep %d: action=%s, num_par=%d, ll=%0.2f, rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f' % \
			(step_i, action, num_par, ll, rss, bic, extended_bic, modified_bic)
		print 'Cofactors:', _cofactors_to_string_(cofactors)

	lm_res = lm.fast_f_test(snps)
	min_pval_i = sp.argmin(lm_res['ps'])
	min_pval = lm_res['ps'][min_pval_i]
	min_pval_chr_pos = chr_pos_list[min_pval_i]
	print 'Min p-value:', min_pval
	step_info = {'rss':rss, 'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic,
		'mbonf':max_cofactor_pval, 'cofactors':map(tuple, cofactors[:]),
		'cofactor_snps':cofactor_snps[:], 'min_pval':min_pval, 'min_pval_chr_pos': min_pval_chr_pos,
		'interactions':interactions}
	lm_pvals = lm_res['ps'].tolist()
	if save_pvals:
		step_info['ps'] = lm_pvals

	#Now plotting!
	print "Generating plots"
	if file_prefix:
		_plot_manhattan_and_qq_(file_prefix, step_i, lm_pvals, quantiles_dict, positions=positions,
					chromosomes=chromosomes, mafs=mafs, macs=macs, plot_bonferroni=True, highlight_markers=cofactors,
					cand_genes=cand_gene_list, plot_xaxis=plot_xaxis, log_qq_max_val=log_qq_max_val,
					with_qq_plots=with_qq_plots, highlight_loci=highlight_loci, write_pvals=save_pvals)

	max_num_cofactors = len(cofactors)
	step_info['kolmogorov_smirnov'] = agr.calc_ks_stats(lm_pvals)
	step_info['pval_median'] = agr.calc_median(lm_pvals)
	print step_info['kolmogorov_smirnov'], step_info['pval_median']
	step_info_list.append(step_info)

	#Now backward stepwise.
	if forward_backwards:
		print 'Starting backwards..'
		while len(cofactor_snps) > 1:
			step_i += 1
			f_stats = sp.zeros(len(cofactor_snps))
			for i, snp in enumerate(cofactor_snps):
				t_cofactors = cofactor_snps[:]
				del t_cofactors[i]
				lm.set_factors(t_cofactors)
				res = lm.fast_f_test([snp])
				cofactors[i][2] = res['ps'][0]
				f_stats[i] = res['f_stats'][0]
			i_to_remove = f_stats.argmin()
			del cofactor_snps[i_to_remove]
			del cofactors[i_to_remove]
			lm.set_factors(cofactor_snps)
			num_snps += 1

			#Re-estimating the REML and ML.
			rss = lm.get_rss()
			ll = lm.get_ll(rss)
			num_par -= 1
			action = '-'

			#Update the p-values
			cofactor_pvals = []
			for i, snp in enumerate(cofactor_snps):
				t_cofactors = cofactor_snps[:]
				del t_cofactors[i]
				lm.set_factors(t_cofactors)
				res = lm.fast_f_test([snp])
				pval = res['ps'][0]
				cofactor_pvals.append(pval)
				cofactors[i][2] = pval
			max_cofactor_pval = max(cofactor_pvals)
			criterias['mbonf'].append(max_cofactor_pval)

			#Calculate the BICs
			(bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lm.n)
			criterias['ebics'].append(extended_bic)
			criterias['mbics'].append(modified_bic)
			print '\nStep %d: action=%s, num_par=%d, ll=%0.2f, rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f' % \
				(step_i, action, num_par, ll, rss, bic, extended_bic, modified_bic)
			print 'Cofactors:', _cofactors_to_string_(cofactors)

			step_info = {'rss':rss, 'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic,
				'mbonf':max_cofactor_pval, 'cofactors':map(tuple, cofactors[:]),
				'cofactor_snps':cofactor_snps[:], 'min_pval':None, 'min_pval_chr_pos':None,
				'kolmogorov_smirnov':None, 'pval_median':None}
			step_info_list.append(step_info)
			print cofactors

	opt_dict, opt_indices = _analyze_opt_criterias_(criterias, sign_threshold, max_num_cofactors, file_prefix, with_qq_plots, lm,
				step_info_list, chr_pos_list, quantiles_dict, plot_bonferroni=True, cand_genes=cand_gene_list,
				plot_xaxis=plot_xaxis, log_qq_max_val=log_qq_max_val, type='lm', highlight_loci=highlight_loci,
				write_pvals=save_pvals, **kwargs)

	for step_i in opt_indices:
		for h in ['min_pval', 'min_pval_chr_pos', 'kolmogorov_smirnov', 'pval_median']:
			step_info_list[step_i][h] = opt_indices[step_i][h]

	if file_prefix:
		_plot_stepwise_stats_(file_prefix, step_info_list, sign_threshold, type == 'lm')

	res_dict = {'step_info_list':step_info_list, 'first_lm_res':first_lm_res, 'opt_dict':opt_dict}

	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)
	return res_dict





def _get_joint_var_matrix_(trait_corr_mat, joint_ecotypes, ecotype_maps, lmm_list, K):
	gen_var_list = []
	err_var_list = []
	for lmm in lmm_list:
		res = lmm.get_REML()
		print 'Pseudo-heritatbility:', res['pseudo_heritability']
		her_list.append(res['pseudo_heritability'])
		err_var_list.append(res['ve'])
		gen_var_list.append(res['vg'])

	#Construct new variance matrix
	V = sp.zeros((len(joint_ecotypes), len(joint_ecotypes)))
	m_i = 0
	for i, pid1 in enumerate(pids):
		ets_map1 = ecotype_maps[pid1]
		num_ets1 = len(ets_map1)
		m_j = 0
		for j, pid2 in enumerate(pids):
			ets_map2 = ecotype_maps[pid2]
			num_ets2 = len(ets_map2)
			if i == j:
				var_sclice = math.sqrt(gen_var_list[i] * gen_var_list[j]) * \
						K[ets_map1, :][:, ets_map2] + err_var_list[i] * sp.eye(num_ets1)
			else:
				if  her_list[i] > 0 and her_list[j] > 0:
					rho = trait_corr_mat[i, j] / sp.sqrt(her_list[i] * her_list[j])
					if rho > 1.0: rho = 1.0
					if rho < -1.0: rho = 1.0
				else:
					rho = 0
				print i, j, rho
				var_sclice = rho * math.sqrt(gen_var_list[i] * gen_var_list[j]) * \
						K[ets_map1, :][:, ets_map2]
			V[m_i:m_i + num_ets1, m_j:m_j + num_ets2] = var_sclice
			m_j += num_ets2
		m_i += num_ets1

	print'Performing Cholesky decomposition'
	H_sqrt = cholesky(V)
	H_sqrt_inv = sp.mat(H_sqrt.T).I
	return H_sqrt_inv



def linear_model(snps, phenotypes, cofactors=None):
	lm = LinearModel(phenotypes)
	if cofactors:
		for cofactor in cofactors:
			lm.add_factor(cofactor)
	print "Running a standard linear model"
	s1 = time.time()
	res = lm.fast_f_test(snps)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)
	return res



def emmax_two_snps(snps, phenotypes, K, cofactors=None):
	lmm = LinearMixedModel(phenotypes)
	lmm.add_random_effect(K)
	if cofactors:
		for cofactor in cofactors:
			lm.add_factor(cofactor)

	print "Running EMMAX on pairs of SNPs"
	s1 = time.time()
	res = lmm.emmax_two_snps(snps)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)
	print 'pseudo_heritability:', res['pseudo_heritability']
	return res



def linear_model_two_snps(snps, phenotypes, cofactors=None):
	lm = LinearModel(phenotypes)
	if cofactors:
		for cofactor in cofactors:
			lm.add_factor(cofactor)

	print "Running standard regression on pairs of SNPs"
	s1 = time.time()
	res = lm.two_snps_ftest(snps)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)
	return res



def emmax_snp_pair_plot(snps, positions, phenotypes, K, fm_scatter_plot_file=None,
			fm_image_plot_file=None, scatter_plot_file=None,
			image_plot_file=None, fig_format='pdf', show_plots=False,
			display_threshold=0.05):
	"""
	Plots a 2D plot.  
	
	Assumes phenotypes and genotypes are already coordinated.
	"""
	r = emmax_two_snps(snps, phenotypes, K)
	thres = -sp.log10(display_threshold)

	print 'Plotting scatter plot.'
	import pylab
	pylab.figure(figsize=(14, 12))
	start_pos = min(positions)
	end_pos = max(positions)
	x_range = end_pos - start_pos
	pylab.axis([start_pos - 0.025 * x_range, end_pos + 0.025 * x_range, start_pos - 0.025 * x_range, end_pos + 0.025 * x_range])

	p_val_array = sp.array(map(lambda y: map(lambda x:-sp.log10(x), y), r['ps']))

	xs = []
	ys = []
	cs = []
	for i, pos, in enumerate(positions):
		xs.extend([pos] * len(positions))
		ys.extend(positions)
		cs.extend(p_val_array[i, :])

	for i in range(len(xs) - 1, -1, -1):
		if cs[i] < thres:
			del xs[i]
			del ys[i]
			del cs[i]

	l = zip(cs, xs, ys)
	l.sort()
	l = map(list, zip(*l))
	cs = l[0]
	xs = l[1]
	ys = l[2]

	pylab.plot([start_pos, end_pos], [start_pos, end_pos], 'k--')
	pylab.scatter(xs, ys, c=cs, alpha=0.8, linewidths=0)
	pylab.colorbar()
	if show_plots:
		 pylab.show()
	if scatter_plot_file:
		pylab.savefig(fm_scatter_plot_file, format=fig_format)
	pylab.clf()
	pylab.imshow(p_val_array, interpolation='nearest')
	pylab.colorbar()
	if show_plots:
		 pylab.show()
	if image_plot_file:
		pylab.savefig(fm_image_plot_file, format=fig_format)

	pylab.clf()
	#Now Vincent's way..  not full model p-vals.
	p_val_array = sp.array(map(lambda y: map(lambda x:-sp.log10(x), y), r['vincent_ps']))

	xs = []
	ys = []
	cs = []
	for i, pos, in enumerate(positions):
		xs.extend([pos] * len(positions))
		ys.extend(positions)
		cs.extend(p_val_array[i, :])

	for i in range(len(xs) - 1, -1, -1):
		if cs[i] < thres:
			del xs[i]
			del ys[i]
			del cs[i]

	l = zip(cs, xs, ys)
	l.sort()
	l = map(list, zip(*l))
	cs = l[0]
	xs = l[1]
	ys = l[2]

	pylab.plot([start_pos, end_pos], [start_pos, end_pos], 'k--')
	pylab.scatter(xs, ys, c=cs, alpha=0.8, linewidths=0)
	pylab.axis([start_pos - 0.025 * x_range, end_pos + 0.025 * x_range, start_pos - 0.025 * x_range, end_pos + 0.025 * x_range])
	pylab.colorbar()
	if show_plots:
		 pylab.show()
	if scatter_plot_file:
		pylab.savefig(scatter_plot_file, format=fig_format)
	pylab.clf()
	pylab.imshow(p_val_array, interpolation='nearest')
	pylab.colorbar()
	if show_plots:
		 pylab.show()
	if image_plot_file:
		pylab.savefig(image_plot_file, format=fig_format)



def prepare_k(k, k_accessions, accessions):
	if k_accessions == accessions:
		return sp.mat(k)
	indices_to_keep = []
	for acc in accessions:
		try:
			i = k_accessions.index(acc)
			indices_to_keep.append(i)
		except:
			continue
	k = k[indices_to_keep, :][:, indices_to_keep]
	return sp.mat(k)

def scale_k(k, verbose=False):
	c = sp.sum((sp.eye(len(k)) - (1.0 / len(k)) * sp.ones(k.shape)) * sp.array(k))
	scalar = (len(k) - 1) / c
	if verbose:
		print 'Kinship scaled by: %0.4f' % scalar
	k = scalar * k
	return k

def load_kinship_from_file(kinship_file, accessions=None, dtype='double', return_accessions=False, scaled=True):
	assert os.path.isfile(kinship_file), 'File not found.'
	#sys.stdout.write("Loading K.\n")
	#sys.stdout.flush()
	f = open(kinship_file, 'r')
	l = cPickle.load(f)
	f.close()
	k = sp.matrix(l[0], dtype=dtype)
	k_accessions = l[1]
	if accessions:
		k = prepare_k(k, k_accessions, accessions)
	if scaled:
		k = scale_k(k)
	if return_accessions:
		return k, k_accessions
	else:
		return k


def save_kinship_to_file(kinship_file, kinship_mat, k_accessions):
	with open(kinship_file, 'wb') as f:
		cPickle.dump([sp.array(kinship_mat).tolist(), k_accessions], f)



def save_kinship_in_text_format(filename, k, accessions):
	with open(filename, 'w') as f:
		for acc, row in it.izip(accessions, k):
			f.write('%s,%s\n' % (acc, ','.join(map(str, row.tolist()))))



#
#def _check_emma_bug_():
#	"""
#	It seems fixed?
#	"""
#	import dataParsers as dp
#	import phenotypeData as pd
#	filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_raw_112210.csv"
#	phed = pd.parse_phenotype_file(filename)
#	phed.convert_to_averages()
#	sd = dp.parse_numerical_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t72.csv.binary', filter=1)
#	sd.coordinate_w_phenotype_data(phed, 5)
#	phenotypes = phed.get_values(5)
#	ets = phed.get_ecotypes(5)
#	K = load_kinship_from_file('/Users/bjarnivilhjalmsson/Projects/Data/250k/kinship_matrix_cm72.pickled',
#				ets)
#	lmm = LinearMixedModel(phenotypes)
#	lmm.add_random_effect(K)
#	print 'Getting REML'
#	res, lls, dlls, deltas = lmm.get_REML()
#	import pylab
#	pylab.plot(deltas, dlls)
#	pylab.ylabel('dlls')
#	pylab.xlabel('log(delta)')
#	pylab.savefig('/Users/bjarni.vilhjalmsson/tmp/dlls.png', format='png')
#	pylab.clf()
#	pylab.plot(deltas, lls)
#	pylab.ylabel('lls')
#	pylab.xlabel('log(delta)')
#	pylab.savefig('/Users/bjarni.vilhjalmsson/tmp/lls.png', format='png')


def local_vs_global_mm(y, local_k, global_k, K):
	"""
	Local vs. global kinship mixed model.
	"""
	lmm0 = LinearMixedModel(Y=y)
	lmm0.add_random_effect(K)
	eig_L = lmm0._get_eigen_L_()
	h0_res = lmm0.get_estimates(eig_L)

	lmm = LinearMixedModel(Y=y)
	lmm.add_random_effect(global_k)
	lmm.add_random_effect(local_k)
	h1_res = lmm.get_estimates_3()
	lrt_stat = 2 * (h1_res['max_ll'] - h0_res['max_ll'])
	pval = stats.chi2.sf(lrt_stat, 1)
	#print 'p-value:', pval
	res_dict = {'pval':pval, 'perc_var1':h1_res['perc_var1'], 'perc_var2':h1_res['perc_var2'],
			'pseudo_heritability0':h0_res['pseudo_heritability'],
			'pseudo_heritability1':h1_res['pseudo_heritability']}
	return res_dict


def local_vs_global_mm_scan(y, sd, file_prefix='/tmp/temp', window_size=1000000, jump_size=500000, kinship_method='ibd', global_k=None):
	"""
	Local vs. global kinship mixed model.
	"""
	print 'Starting Mixed model, local vs. global kinship scan...'
	print 'window size is %d, and jump size is %d' % (window_size, jump_size)
	import gwaResults as gr
	if global_k == None:
		if kinship_method == 'ibd':
			K = sd.get_ibd_kinship_matrix()
		elif kinship_method == 'ibs':
			K = sd.get_ibs_kinship_matrix()
		else:
			raise NotImplementedError
	else:
		K = global_k
	genome_length = sd.get_genome_length()
	est_num_chunks = genome_length / jump_size
	chromosomes = []
	positions = []
	pvals = []
	perc_variances1 = []
	perc_variances2 = []
	h1_heritabilities = []
	h0_heritabilities = []
	chunk_i = 0
	for sdl, chrom in zip(sd.snpsDataList, sd.chromosomes):
		max_pos = sdl.positions[-1]
		pos = sdl.positions[0]
		for focal_pos in range(sdl.positions[0], sdl.positions[-1], jump_size):
			chunk_i += 1
			d = sd.get_local_n_global_kinships((chrom, focal_pos), window_size,
										global_kinship=K,
										kinship_method=kinship_method)
			local_k = scale_k(d['local_k'])
			global_k = scale_k(d['global_k'])
			if local_k != None and global_k != None:
				#print "Chromosome=%d, position=%d" % (chrom, focal_pos)
				res_dict = local_vs_global_mm(y, local_k, global_k, K)
				chromosomes.append(chrom)
				positions.append(focal_pos)
				perc_variances1.append(res_dict['perc_var1'])
				perc_variances2.append(res_dict['perc_var2'])
				h0_heritabilities.append(res_dict['pseudo_heritability0'])
				h1_heritabilities.append(res_dict['pseudo_heritability1'])
				pvals.append(res_dict['pval'])

				#print 'H0: pseudo_heritability=%0.2f' % (res_dict['h0_res']['pseudo_heritability'])
				#print 'H1: pseudo_heritability=%0.2f, perc_var1=%0.2f, perc_var2=%0.2f' % \
				#		(res_dict['h1_res']['pseudo_heritability'],
				#		res_dict['h1_res']['perc_var1'],
				#		res_dict['h1_res']['perc_var2'])
			if est_num_chunks >= 100 and (chunk_i + 1) % int(est_num_chunks / 100) == 0: #Print dots
				sys.stdout.write('.')
				sys.stdout.flush()

	pval_res = gr.Result(scores=pvals, positions=positions, chromosomes=chromosomes)
	pval_res.neg_log_trans()
	pval_res.plot_manhattan(png_file=file_prefix + '_lrt_pvals.png', percentile=0, plot_bonferroni=True)
	perc_var_res = gr.Result(scores=perc_variances2, positions=positions, chromosomes=chromosomes)
	perc_var_res.plot_manhattan(png_file=file_prefix + '_perc_var_explained.png', percentile=0,
				ylab='% of variance explained')






def _emmax_local_global_kinship_test_(pid):
	import dataParsers as dp
	import phenotypeData as pd
	import env
	sd = dp.load_snps_call_method(75)
	phed = pd.get_phenotypes_from_db([pid])
	phed.convert_to_averages()
	print phed.get_name(pid)
	sd.coordinate_w_phenotype_data(phed, pid)
#	local_k, global_k = sd.get_local_n_global_kinships((4, 269000), 1000000)
#	local_k = scale_k(local_k)
#	global_k = scale_k(global_k)
#	print 'Local kinship matrix:', local_k
#	print 'Global kinship matrix:', global_k
	Y = phed.get_values(pid)
	local_vs_global_mm_scan(Y, sd)
	#local_vs_global_mm(Y, local_k, global_k,)



def _emmax_rna_3_matrices_test_(pid):
	import dataParsers as dp
	import phenotypeData as pd
	import env
	sd = dp.load_snps_call_method(75)
	K_pop, K_accessions = load_kinship_from_file('/Users/bjarni.vilhjalmsson/Projects/Data/rna_seq/K_pop_10.pickled', sd.accessions, return_accessions=True)
	sd.filter_accessions(K_accessions)
	phed = pd.parse_phenotype_file('/Users/bjarni.vilhjalmsson/Projects/Data/phenotypes/rna_seq_vs_081411_10C.csv')
	sd.coordinate_w_phenotype_data(phed, pid)
	phed.plot_histogram(pid, png_file='/Users/bjarni.vilhjalmsson/tmp/rna_3k_test_%s_hist.png' % phed.get_name(pid))

	#load kinship matrix
#	K_pop = load_kinship_from_file('/Users/bjarni.vilhjalmsson/Projects/Data/rna_seq/K_pop_10.pickled', sd.accessions)
#	K_expr = load_kinship_from_file('/Users/bjarni.vilhjalmsson/Projects/Data/rna_seq/K_expr_10.pickled', sd.accessions)
	#Run stepwise EMMAX!!!
	emmax_step_wise(phed.get_values(pid), K_pop, sd=sd, K2=K_expr,
			file_prefix='/Users/bjarni.vilhjalmsson/tmp/rna_3k_test_%s' % phed.get_name(pid))
	# FINISH!!



def perform_human_emmax(pid=1):
	import dataParsers as dp
	import phenotypeData as pd
	import env
	import gwaResults as gr
	import random
	import sys
	s1 = time.time()
	plink_prefix = env.env['phen_dir'] + 'NFBC_20091001/NFBC_20091001'
	sd = dp.parse_plink_tped_file(plink_prefix)
	#sd.sample_snps(0.05)
	individs = sd.accessions[:]
	phed = pd.parse_phenotype_file(env.env['phen_dir'] + 'NFBC_20091001/phenotype.csv')
	#phed.filter_ecotypes(pid, random_fraction=0.2)
	sd.coordinate_w_phenotype_data(phed, pid)
	sd.filter_mac_snps(50)
	K = sd.get_ibd_kinship_matrix(num_dots=10)
	#K = prepare_k(K, individs, sd.accessions)
	phen_vals = phed.get_values(pid)
	phen_name = phed.get_name(pid)
	print 'Working on %s' % phen_name
	sys.stdout.flush()
	file_prefix = env.env['results_dir'] + 'NFBC_emmax_step_%s_pid%d' % (phen_name, pid)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds to load and preprocess the data.' % (mins, secs)
	else:
		print 'Took %f seconds to load and preprocess the data..' % (secs)
	emmax_res = emmax_step_wise(phen_vals, K, sd=sd, num_steps=10, file_prefix=file_prefix, plot_xaxis=False)
	#snps = sd.getSnps()
	#emmax_res = emmax(snps, phen_vals, K)
	#res = gr.Result(scores=emmax_res['ps'].tolist(), snps_data=sd)
	#res.write_to_file(env.env['results_dir'] + 'NFBC_emmax_pid%d.pvals' % pid)
	#res.neg_log_trans()
	#res.plot_manhattan(png_file=file_prefix + '.png', plot_xaxis=False, plot_bonferroni=True)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds in total.' % (mins, secs)
	else:
		print 'Took %f seconds in total.' % (secs)



def test_skin_color():
	s1 = time.time()
	import dataParsers as dp
	import phenotype_parsers as pp
	import env
	import gwaResults as gr
	for pid in [1, 2]:
		#pid = 2
		dir_prefix = env.env['home_dir'] + 'Projects/data/skin_eye_color/'
		#dir_prefix = env.env['home_dir'] + 'Projects/Data/Skin_color/'
		plink_prefix = dir_prefix + 'plink'
		sd = dp.parse_plink_tped_file(plink_prefix)
		phed = pp.load_skin_color_traits()
		sd.coordinate_w_phenotype_data(phed, pid)
		sd.filter_mac_snps(30)
		K = load_kinship_from_file(dir_prefix + 'kinship_ibd.pickled',
						accessions=sd.accessions)
		phen_vals = phed.get_values(pid)
		phen_name = phed.get_name(pid)
		print 'Working on %s' % phen_name
		sys.stdout.flush()
		file_prefix = env.env['results_dir'] + 'CVI_ibd_%s_pid%d' % (phen_name, pid)
		snps = sd.getSnps()
		secs = time.time() - s1
		if secs > 60:
			mins = int(secs) / 60
			secs = secs - mins * 60
			print 'Took %d mins and %f seconds to load and preprocess the data.' % (mins, secs)
		else:
			print 'Took %f seconds to load and preprocess the data..' % (secs)
		emmax_res = emmax_step_wise(phen_vals, K, sd=sd, num_steps=10, file_prefix=file_prefix, plot_xaxis=False)
		#emmax_res = emmax(snps, phen_vals, K)
		#phed.plot_histogram(pid, png_file=file_prefix + '_hist.png', p_her=emmax_res['pseudo_heritability'])
#		res = gr.Result(scores=emmax_res['ps'].tolist(), snps_data=sd)
		#emma_res = emma(snps, phen_vals, K)
#		phed.plot_histogram(pid, png_file=file_prefix + '_hist.png', p_her=emmax_res['pseudo_heritability'])
#		res = gr.Result(scores=emmax_res['ps'].tolist(), snps_data=sd)
#		res.write_to_file(env.env['results_dir'] + 'CVI_emmax_%s_pid%d.pvals' % (phen_name, pid))
#		res.neg_log_trans()
#		res.plot_manhattan(png_file=file_prefix + '.png', plot_xaxis=False, plot_bonferroni=True)
		secs = time.time() - s1
		if secs > 60:
			mins = int(secs) / 60
			secs = secs - mins * 60
			print 'Took %d mins and %f seconds in total.' % (mins, secs)
		else:
			print 'Took %f seconds in total.' % (secs)


def _test_emmax_step_():
	import dataParsers as dp
	import phenotypeData as pd
	import env
	pid = 226
	sd = dp.load_250K_snps()
	phed = pd.parse_phenotype_file(env.env['phen_dir'] + 'phen_raw_112210.csv')
	phed.convert_to_averages()
	print phed.get_name(pid)
	sd.coordinate_w_phenotype_data(phed, pid)
	K = load_kinship_from_file(env.env['data_dir'] + 'kinship_matrix_cm72.pickled', accessions=sd.accessions)
	phen_vals = phed.get_values(pid)
	cpl = sd.getChrPosList()
	cof_chr_pos_list = [(4, 6392280), (4, 6418442)]#cpl[100]]
	emmax_res = emmax_step(phen_vals, sd, K, cof_chr_pos_list)
	pdb.set_trace()


def _test_mtmm_():
	import dataParsers as dp
	import phenotypeData as pd
	import env
	import analyze_gwas_results as agr
	pids = [314, 315, 316, 317]
	result_file_prefix = 'mtmm_%s' % ('_'.join(map(str, pids)))
	trans_types = ['most_normal', 'most_normal', 'most_normal', 'most_normal']
	sd = dp.load_250K_snps()
	phed = pd.parse_phenotype_file(env.env['phen_dir'] + 'phen_raw_112210.csv')
	k, k_ets = dp.load_kinship(return_accessions=True, scaled=False)
	phed.convert_to_averages()
	for pid, tt in zip(pids, trans_types):
		phed.transform(pid, tt)#, method='none')
	unique_ets = phed.filter_ecotypes_2(sd.accessions, pids)
	phed.normalize_values(pids)
	phenotype_list = []
	ecotype_list = []
	for pid in pids:
		print phed.get_name(pid)
		phenotype_list.append(phed.get_values(pid))
		ecotype_list.append(phed.get_ecotypes(pid))
	mtmm = MultipleTraitLMM(phenotype_list, ecotype_list, K=k, K_ets=k_ets, unique_ets=unique_ets)
	sd.filter_accessions(list(unique_ets))
	snps = sd.getSnps()
	r = mtmm.emmax_f_test_fm(snps)
	import gwaResults as gr
	res = gr.Result(snps_data=sd, scores=r['ps'])
	res.plot_manhattan(png_file=env.env['tmp_dir'] + result_file_prefix + '_manhattan.png', neg_log_transform=True)
	log_quantiles = agr.get_log_quantiles(r['ps'])
	agr.simple_log_qqplot([log_quantiles], png_file=env.env['tmp_dir'] + result_file_prefix + '_qq_plot.png')
	res_full = gr.Result(snps_data=sd, scores=r['ps_full'])
	res_full.plot_manhattan(png_file=env.env['tmp_dir'] + result_file_prefix + '_manhattan_full.png', neg_log_transform=True)
	log_quantiles = agr.get_log_quantiles(r['ps_full'])
	agr.simple_log_qqplot([log_quantiles], png_file=env.env['tmp_dir'] + result_file_prefix + '_qq_plot_full.png')
	res_ge = gr.Result(snps_data=sd, scores=r['ps_ge'])
	res_ge.plot_manhattan(png_file=env.env['tmp_dir'] + result_file_prefix + '_manhattan_ge.png', neg_log_transform=True)
	log_quantiles = agr.get_log_quantiles(r['ps_ge'])
	agr.simple_log_qqplot([log_quantiles], png_file=env.env['tmp_dir'] + result_file_prefix + '_qq_plot_ge.png')
	pdb.set_trace()


def test_bayes_factor_enrichment():
	import dataParsers as dp
	import gwaResults as gr
	import phenotypeData as pd
	import env
	cpp_list = []
	with open('/Users/bjarni.vilhjalmsson/Projects/Data/DTF1.scan.tsv') as f:
		print f.next()
		for l in f:
			line = l.split()
			cpp_list.append((int(line[1]), int(line[2]), float(line[4])))

	for pid in [314, 315, 316, 317]:
		sd = dp.load_snps_call_method(75)
		phed = pd.parse_phenotype_file(env.env['phen_dir'] + 'phen_raw_112210.csv')
		#phed = pd.get_phenotypes_from_db([pid])
		phed.convert_to_averages()
		phed.transform(pid, 'most_normal')
		sd.coordinate_w_phenotype_data(phed, pid)
		sd.filter_mac_snps(15)
		phen_vals = phed.get_values(pid)
		phen_name = phed.get_name(pid)
		K = dp.load_kinship(75, accessions=sd.accessions, sd=sd)

		#Candidate genes (TAIR IDs)
		cg_tair_ids = ['AT1G04400', 'AT1G65480', 'AT1G77080', 'AT2G26330', 'AT4G00650', 'AT5G10140', 'AT5G65050', \
				'AT5G65060', 'AT5G65070', 'AT5G65080'] \
				#FT, MAF, ER, FRI, FLC, MAF2-MAF5 (Salome et al. 2011)
		cgs = gr.get_genes_w_tair_id(cg_tair_ids)
		#snp_priors = sd.get_cand_genes_snp_priors(cgs, radius=25000, cg_prior_fold_incr=20)
		snp_priors = sd.get_snp_priors(cpp_list, cand_genes=cgs)
		emmax_step_wise(phen_vals, K, sd, 10, env.env['tmp_dir'] + 'bf_qtl_priors_%s_%d' % (phen_name, pid),
				snp_priors=snp_priors, cand_gene_list=cgs, save_pvals=True)





if __name__ == "__main__":
#	import env
#	kinship_file_name = env.env['data_dir'] + 'kinship_matrix_cm75.pickled'
#	k, k_accessions = cPickle.load(open(kinship_file_name))
#	save_kinship_in_text_format(env.env['data_dir'] + 'kinship_matrix_cm75.csv', k, k_accessions)
	#_emmax_local_global_kinship_test_(5)
	perform_human_emmax()
	#_test_joint_analysis_()
	#_test_phyB_snps_()
