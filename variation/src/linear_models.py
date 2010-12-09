"""
Contains functions to perform various linear regression schemes, such as simple, and mixed models.
"""
#import numpy as np
#from scipy import *
import scipy as sp
sp.alterdot()
from scipy import linalg
from scipy import stats
from scipy import optimize
import snpsdata
import warnings
import pdb
import cPickle
import os
import sys
#import multiprocessing as mp
import time
from pdb import Pdb


anti_decoder = {1:0, 0:1}
get_anti_snp = sp.vectorize(lambda x: anti_decoder[x])  #Creating a vectorized function for anti-snp


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
			self.X = sp.mat(sp.repeat(1, self.n)).T #The intercept
			self.p = 1
		if len(factors) > 0:
			self.X = sp.hstack([self.X, sp.matrix(factors).T])
			self.p += len(factors)





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

#	def fast_f_test(self, snps):
#		"""
#		A standard linear model, using a F-test
#		"""
#		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(self.X, self.Y)
#		num_snps = len(snps)
#		rss_list = sp.repeat(h0_rss, num_snps)
#		h0_betas = map(float, list(h0_betas)) + [0.0]
#		betas_list = [h0_betas] * num_snps
#		var_perc = sp.zeros(num_snps)
#		q = 1
#		n_p = self.n - self.p
#		for i, snp in enumerate(snps):
#			(betas, rss, rank, s) = linalg.lstsq(sp.hstack([self.X, sp.matrix(snp).T]), self.Y)
#			if not rss:
#				print 'No predictability in the marker, moving on...'
#				continue
#			rss_list[i] = rss[0]
#			betas_list[i] = map(float, list(betas))
#		rss_ratio = h0_rss / rss_list
#		var_perc = 1 - 1 / rss_ratio
#		f_stats = (rss_ratio - 1) * n_p / float(q)
#		p_vals = stats.f.sf(f_stats, q, n_p)
#		return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list,
#			'var_perc':var_perc, 'h0_rss':h0_rss}


	def fast_f_test(self, snps, verbose=True, Z=None,
			with_betas=False):
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
		Y = self.Y	#The transformed outputs.
		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
		Y = sp.mat(Y - h0_X * h0_betas, dtype=dtype)
		h0_betas = map(float, list(h0_betas))

		if not with_betas:
			(Q, R) = linalg.qr(h0_X, econ=True)  #Do the QR-decomposition for the Gram-Schmidt process.
			Q = sp.mat(Q)
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
	def __init__(self, Y=None):
		"""
		The fixed effects should be a list of fixed effect lists (SNPs)
		"""
		self.n = len(Y)
		self.y_var = sp.var(Y, ddof=1)
		self.Y = sp.matrix([[y] for y in Y])
		self.X = sp.matrix([[1] for y in Y]) #The intercept
		self.p = 1
		self.beta_est = None
		self.cofactors = []

		#A list of random effect type, and the cov matrix.
		self.random_effects = [('normal', sp.matrix(sp.identity(self.n)))] #The first random effect is the IID error.


	def add_random_effect(self, cov_matrix=None, effect_type='normal'):
		if effect_type != 'normal':
			raise Exception('Currently, only Normal random effects are allowed.')
		self.random_effects.append((effect_type, cov_matrix))

	def _get_eigen_L_(self, K):
		evals, evecs = linalg.eigh(K)
		return {'values':evals, 'vectors':sp.mat(evecs).T}



	def _get_eigen_R_(self, X, hat_matrix=None):
		q = X.shape[1]
		if not hat_matrix:
			X_squared_inverse = linalg.pinv(X.T * X) #(X.T*X).I
			hat_matrix = X * X_squared_inverse * X.T
		S = sp.mat(sp.identity(self.n)) - hat_matrix	#S=I-X(X'X)^{-1}X'
		evals, evecs = linalg.eigh(S * (self.random_effects[1][1] + self.random_effects[0][1]) * S) #eigen of S(K+I)S
		return {'values':map(lambda x: x - 1, evals[q:]), 'vectors':(sp.mat(evecs).T[q:])}   #Because of S(K+I)S?



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


	def get_REML(self, ngrids=50, llim= -5, ulim=10, esp=1e-6):
		"""
		Get REML estimates for the effect sizes, as well as the random effect contributions.
		
		This is EMMA
		"""
		K = self.random_effects[1][1]
		eig_L = self._get_eigen_L_(K)
		#Get the variance estimates..
		return self.get_estimates(eig_L=eig_L, ngrids=ngrids, llim=llim, ulim=ulim, esp=esp, method='REML')


	def get_ML(self, ngrids=50, llim= -5, ulim=10, esp=1e-6):
		"""
		Get REML estimates for the effect sizes, as well as the random effect contributions.
		
		This is EMMA
		"""
		K = self.random_effects[1][1]
		eig_L = self._get_eigen_L_(K)
		#Get the variance estimates..
		return self.get_estimates(eig_L=eig_L, ngrids=ngrids, llim=llim, ulim=ulim, esp=esp, method='ML')


	def get_estimates(self, eig_L=None, xs=[], ngrids=100, llim= -5, ulim=10, esp=1e-6,
				return_pvalue=False, return_f_stat=False, method='REML', verbose=False):
		"""
		Get ML/REML estimates for the effect sizes, as well as the random effect contributions.
		Using the EMMA algorithm.
		
		Methods available are 'REML', and 'ML'		
		"""
		if verbose:
			print 'Retrieving %s variance estimates' % method
		if not eig_L:
			raise Exception
		if len(xs):
			X = sp.hstack([self.X, xs])
		else:
			X = self.X

		eig_R = self._get_eigen_R_(X)
		q = X.shape[1] #number of columns
		n = self.n
		p = n - q
		m = ngrids + 1

		etas = sp.array(eig_R['vectors'] * self.Y)
		sq_etas = etas * etas
		log_deltas = (sp.arange(m, dtype='float64') / ngrids) * (ulim - llim) + llim  #a list of deltas to search
		deltas = sp.exp(log_deltas)
		assert len(deltas) == m, 'Number of deltas is incorrect.'
		eig_vals = sp.array(eig_R['values'])
		assert len(eig_vals) == p, 'Number of eigenvalues is incorrect.'

  		lambdas = sp.reshape(sp.repeat(eig_vals, m), (p, m)) + \
  			sp.reshape(sp.repeat(deltas, p), (m, p)).T
	  	s1 = sp.sum(sq_etas / lambdas, axis=0)
  		if method == 'REML':
	  		s2 = sp.sum(sp.log(lambdas), axis=0)
	  		lls = 0.5 * (p * (sp.log((p) / (2.0 * sp.pi)) - 1 - sp.log(s1)) - s2)
	 		s3 = sp.sum(sq_etas / (lambdas * lambdas), axis=0)
	 		s4 = sp.sum(1 / lambdas, axis=0)
	  		dlls = 0.5 * (p * s3 / s1 - s4)
  		elif method == 'ML':
  			#Xis < -matrix(eig.L$values, n, m) + matrix(delta, n, m, byrow=TRUE)
  			eig_vals_L = sp.array(eig_L['values'])
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
			opt_ll, opt_i = max(zero_intervals)
			opt_delta = 0.5 * (deltas[opt_i - 1] + deltas[opt_i])
			#Newton-Raphson
			if method == 'REML':
				new_opt_delta = optimize.newton(self._redll_, opt_delta, args=(eig_vals, sq_etas),
								tol=esp, maxiter=50)
				if deltas[opt_i - 1] - esp < new_opt_delta < deltas[opt_i] + esp:
					opt_delta = new_opt_delta
					opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
				else:
					opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
					raise Exception('Local maximum outside of suggested area??')
			elif method == 'ML':
				new_opt_delta = optimize.newton(self._dll_, opt_delta, args=(eig_vals, eig_vals_L, sq_etas),
								tol=esp, maxiter=50)
				if deltas[opt_i - 1] - esp < new_opt_delta < deltas[opt_i] + esp:
					opt_delta = new_opt_delta
					opt_ll = self._ll_(opt_delta, eig_vals, eig_vals_L, sq_etas)
				else:
					opt_ll = self._ll_(opt_delta, eig_vals, eig_vals_L, sq_etas)
					raise Exception('Local maximum outside of suggested area??')
			if opt_ll < max_ll:
				opt_delta = deltas[max_ll_i]
		else:
			opt_delta = deltas[max_ll_i]
			opt_ll = max_ll


		l = sp.reshape(sq_etas, (len(sq_etas))) / (eig_vals + opt_delta)
		opt_vg = sp.sum(l) / p  #vg   
		opt_ve = opt_vg * opt_delta  #ve

		H_sqrt_inv = sp.diag(1.0 / sp.sqrt(eig_L['values'] + opt_delta)) * eig_L['vectors']
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

		if len(xs) and return_f_stat:
			h0_X = H_sqrt_inv * self.X
			(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y_t)
			f_stat = ((h0_rss - mahalanobis_rss) / (xs.shape[1])) / (h0_rss / p)

			res_dict['f_stat'] = float(f_stat)
		if return_pvalue:
			p_val = stats.f.sf(f_stat, (xs.shape[1]), p)
			res_dict['p_val'] = float(p_val)

		return res_dict



#	def get_expedited_REMLE(self, eig_L=None, xs=None, ngrids=100, llim= -5, ulim=10, esp=1e-6,
#				return_pvalue=False, return_f_stat=False):
#		"""
#		Get ML estimates for the effect sizes, as well as the random effect contributions.
#		
#		Using the algorithm adopted in EMMA [Kang et al. (2008)].		
#		"""
#		if not eig_L:
#			raise Exception
#		if xs:
#			X = sp.hstack([self.X, matrix([[v] for x in xs for v in x ])])
#		else:
#			X = self.X
#		eig_R = self._get_eigen_R_(X)
#		K = self.random_effects[1][1]
#		t = K.shape[0] #number of rows
#		q = X.shape[1] #number of columns
#		#print 'X, number of columns: %d, number of rows: %d' % (q, t)
#		n = self.n
#		p = n - q
#		assert K.shape[0] == K.shape[1] == X.shape[0] == n, 'Dimensions are wrong.'
#
#		etas = eig_R['vectors'] * self.Y
#		sq_etas = [float(eta) * float(eta) for eta in etas]
#		log_deltas = [float(i) / ngrids * (ulim - llim) + llim  for i in range(ngrids + 1)] #a list of deltas to search
#		assert len(log_deltas) == ngrids + 1, 'Delta list size error.'
#		deltas = map(sp.exp, log_deltas)
#		eig_vals = list(eig_R['values'])
#		assert len(eig_vals) == p, 'Number of eigenvalues is incorrect.'
#
#		#LL <- 0.5*(p*(log((p)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
#
#		c_1 = 0.5 * p * (sp.log((p) / (2.0 * sp.pi)) - 1)
#
#
#		def calc_ll(d):
#			"""
#			Calculates the likelihood, and the derivative given a delta.
#			"""
#			sum_1 = 0
#			sum_2 = 0
#			sum_3 = 0
#			sum_4 = 0
#			for j in range(p):
#				v_1 = eig_vals[j] + d
#				v_2 = sq_etas[j] / v_1
#				sum_1 += v_2
#				sum_2 += sp.log(v_1)
#				sum_3 += v_2 / v_1
#				sum_4 += 1.0 / v_1
#			#LL <- 0.5*((p)*(log((p)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
#			ll = c_1 - (0.5) * ((p) * sp.log(sum_1) + sum_2)  #log-likelihoods (eq. 7 from paper)
#			#dLL <- 0.5*((p)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
#			dll = 0.5 * ((p) * sum_3 / sum_1 - sum_4)  #diffrentiated log-likelihoods (eq. 9 from paper)
#			return ll, dll
#
#
#		last_ll, last_dll = calc_ll(deltas[0])
#		max_ll = last_ll
#		max_ll_i = 0
#		zero_intervals = []
#		lls = [max_ll]
#		dlls = [last_dll]
#		for i, d in enumerate(deltas[1:]):
#			ll, dll = calc_ll(d)
#			lls.append(ll)
#			dlls.append(dll)
#			if ll > max_ll:
#				max_ll = ll
#				max_ll_i = i
#			if last_dll > 0 and dll < 0:
#				zero_intervals.append(((ll + last_ll) * 0.5, i))
#			last_dll = dll
#		#pdb.set_trace()
#
#		if len(zero_intervals) > 0:
#			#print 'local maximum found:',zero_intervals
#			opt_ll, opt_i = max(zero_intervals)
#			opt_delta = 0.5 * (deltas[opt_i] + deltas[opt_i + 1])
#			#Newton-Raphson
#			new_opt_delta = optimize.newton(self._redll_, opt_delta, args=(eig_vals, sq_etas), tol=esp, maxiter=50)
#			if deltas[opt_i] - esp < new_opt_delta < deltas[opt_i + 1] + esp:
#				opt_delta = new_opt_delta
#				opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
#			else:
#				opt_delta = new_opt_delta
#				opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
#				print 'Local maximum outside of suggested area??'
#				print opt_delta, new_opt_delta
#				raise Exception('Local maximum outside of suggested area??')
#			if opt_ll < max_ll:
#				opt_delta = deltas[max_ll_i]
#		else:
#			opt_delta = deltas[max_ll_i]
#			opt_ll = max_ll
#
#
#		#eig_vals.reverse()
#		l = map(lambda x, y: x / (y + opt_delta), sq_etas, eig_vals)
#		opt_vg = sum(l) / (p)  #vg   
#		opt_ve = opt_vg * opt_delta  #ve
#
#		H_sqrt = eig_L['vectors'].T * sp.diag([sp.sqrt(ev + opt_delta) for ev in eig_L['values']])
#		H_inverse = eig_L['vectors'].T * sp.diag([1.0 / (ev + opt_delta) for ev in eig_L['values']]) * eig_L['vectors']
#		H_sqrt_inv = eig_L['vectors'].T * sp.diag([1.0 / sp.sqrt(ev + opt_delta) for ev in eig_L['values']])
#		XX = X.T * (H_inverse * X)
#		iXX = linalg.pinv(XX)
#		beta_est = iXX * X.T * (H_inverse * self.Y)
#		x_beta = X * beta_est
#		residuals = self.Y - x_beta
#		mahalanobis_rss = residuals.T * H_inverse * residuals
#		rss = residuals.T * residuals
#		x_beta_var = sp.var(x_beta, ddof=1)
#		var_perc = x_beta_var / self.y_var
#		res_dict = {'max_ll':opt_ll, 'delta':opt_delta, 'beta':beta_est, 've':opt_ve, 'vg':opt_vg,
#			    'var_perc':var_perc, 'rss':rss, 'mahalanobis_rss':mahalanobis_rss, 'H_sqrt':H_sqrt,
#			    'H_sqrt_inv':H_sqrt_inv, 'pseudo_heritability':1.0 / (1 + opt_delta)}
#		if xs and return_f_stat:
#			q = X.shape[1] - self.X.shape[1]
#			h0_iXX = linalg.pinv(self.X.T * (H_inverse * self.X))
#			h0_beta_est = h0_iXX * self.X.T * (H_inverse * self.Y)
#			h0_residuals = self.Y - self.X * h0_beta_est
#			h0_rss = h0_residuals.T * H_inverse * h0_residuals
#			f_stat = ((h0_rss - rss) / (q)) / (rss / p)
#
#			#Redo this t-statistic.
##			t_stat = self.beta_est[q-1]/sp.sqrt(float(iXX[q-1,q-1])*opt_vg)
#			res_dict['f_stat'] = float(f_stat)
#		if return_pvalue:
#			#p_val = stats.t.sf(abs(t_stat), p)*2
#			p_val = stats.f.sf(f_stat, q, p)
#			res_dict['p_val'] = float(p_val)
#
#		print res_dict
#		return res_dict



	def expedited_REML_t_test(self, snps, ngrids=50, llim= -4, ulim=10, esp=1e-6, verbose=True):
		"""
		Single SNP analysis (Not as fast as the R-version?)
		"""
		assert len(self.random_effects) == 2, "Expedited REMLE only works when we have exactly two random effects."
		K = self.random_effects[1][1]
		eig_L = self._get_eigen_L_(K)
		f_stats = []
		vgs = []
		ves = []
		max_lls = []
		var_perc = []
		betas = []
		p_vals = []
		num_snps = len(snps)

		for i, snp in enumerate(snps):
			res = self.get_estimates(eig_L=eig_L, xs=sp.matrix(snp).T, ngrids=ngrids, llim=llim, ulim=ulim,
						       esp=esp, return_pvalue=True, return_f_stat=True)
			f_stats.append(res['f_stat'])
			vgs.append(res['vg'])
			ves.append(res['ve'])
			max_lls.append(res['max_ll'])
			var_perc.append(res['var_perc'])
			betas.append(map(float, list(res['beta'])))
			p_vals.append(res['p_val'])
			if verbose and num_snps >= 10 and (i + 1) % (num_snps / 10) == 0: #Print dots
				sys.stdout.write('.')
				sys.stdout.flush()


		return {'ps':p_vals, 'f_stats':f_stats, 'vgs':vgs, 'ves':ves, 'var_perc':var_perc,
			'max_lls':max_lls, 'betas':betas}


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
			#print 'Permutation nr. %d' % perm_i
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


	def emmax_f_test(self, snps, Z=None, with_betas=False, method='REML'):
		"""
		EMMAX implementation (in python)
		Single SNPs
		
		With interactions between SNP and possible cofactors.
		"""
		K = self.random_effects[1][1]
		eig_L = self._get_eigen_L_(K)
		#s = time.time()
		res = self.get_estimates(eig_L=eig_L, method=method) #Get the variance estimates..
		#print 'Took % .6f secs.' % (time.time() - s)
		print 'pseudo_heritability:', res['pseudo_heritability']

		r = self._emmax_f_test_(snps, res['H_sqrt_inv'], Z=Z, with_betas=with_betas)
		r['pseudo_heritability'] = res['pseudo_heritability']
		r['max_ll'] = res['max_ll']
		return r




	def _emmax_f_test_(self, snps, H_sqrt_inv, verbose=True, return_transformed_snps=False, Z=None,
			with_betas=False):
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
			(Q, R) = linalg.qr(h0_X, econ=True)  #Do the QR-decomposition for the Gram-Schmidt process.
			Q = sp.mat(Q)
			Q2 = Q * Q.T
			M = sp.mat(H_sqrt_inv.T * (sp.eye(n) - Q2), dtype=dtype)
		else:
			betas_list = [h0_betas] * num_snps

		rss_list = sp.repeat(h0_rss, num_snps)
		if return_transformed_snps:
			t_snps = []
		chunk_size = len(Y)
		for i in range(0, len(snps), chunk_size): #Do the dot-product in chuncks!
			snps_chunk = sp.matrix(snps[i:i + chunk_size])
			if with_betas:
				Xs = snps_chunk * (H_sqrt_inv.T)
			else:
				Xs = sp.mat(snps_chunk, dtype=dtype) * M
			for j in range(len(Xs)):
				if return_transformed_snps:
					t_snps.append(sp.array(Xs[j]).flatten())
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
		if return_transformed_snps:
			res_d['t_snps'] = t_snps
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



#	def emmax_two_snps(self, snps, verbose=True):
#		"""
#		Every pair of SNPs, Vincent's results.
#		"""
#		import util
#		K = self.random_effects[1][1]
#		eig_L = self._get_eigen_L_(K)
#		num_snps = len(snps)
#
#		res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
#		delta = res['delta']
#		pseudo_her = 1.0 / (1 + delta)
#		H_sqr = res['H_sqrt']
#		H_sqrt_inv = H_sqr.I
#		Y = H_sqrt_inv * self.Y	#The transformed outputs.
#		h0_X = H_sqrt_inv * self.X     #The transformed fixed inputs.
#		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
#		full_rss = h0_rss
#
#		#Contructing result containers		
#		f_stat_array = sp.zeros((num_snps, num_snps))
#		rss_array = sp.repeat(full_rss, num_snps * num_snps)
#		rss_array.shape = (num_snps, num_snps)
#		p_val_array = sp.ones((num_snps, num_snps))
#		var_perc_array = sp.zeros((num_snps, num_snps))
#		haplotype_counts = [[{} for j in range(i + 1)] for i in range(num_snps)]
#
#		q = 1 #Testing one SNP, given the other
#		p = len(self.X.T) + q
#		n_p = self.n - p
#		#Fill the diagonal with the single SNP emmax
#		for i, snp in enumerate(snps):
#			hap_set, hap_counts, haplotypes = snpsdata.get_haplotypes([snp], self.n,
#										count_haplotypes=True)
#			d = {'num_haps':hap_counts}
#			for hap, hap_c in zip(hap_set, hap_counts):
#				d[hap] = hap_c
#			haplotype_counts[i][i] = d
#			snp_mat = H_sqrt_inv * sp.matrix(snp).T #Transformed inputs
#			X = sp.hstack([h0_X, snp_mat])
#			(betas, rss, rank, s) = linalg.lstsq(X, Y)
#			if not rss:
#				if verbose: print 'No predictability in the marker, moving on...'
#				continue
#			f_stat = ((h0_rss - rss) / q) / (rss / n_p)
#			p_val = stats.f.sf(f_stat, q, n_p)
#			p_val_array[i, i] = p_val[0]
#			f_stat_array[i, i] = f_stat[0]
#			rss_array[i, i] = rss[0]
#			var_perc_array[i, i] = float(1 - rss / h0_rss)
#
#		#Fill the upper right part with poorest of two tests.
#		p = len(self.X.T) + q + 1 #Adding one SNP as a cofactor.
#		n_p = self.n - p
#		for i, snp1 in enumerate(snps):
#			anti_snp1 = get_anti_snp(snp1)
#			snp_mat = H_sqrt_inv * (sp.matrix(snp1).T) #Transformed inputs
#			h0_X1 = sp.hstack([h0_X, snp_mat])
#			(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X1, Y)
#			if not h0_rss:
#				if verbose: print 'No predictability in the marker, moving on...'
#				continue
#			for j, snp2 in enumerate(snps):
#				if i == j: continue #Skip diagonals..
#
#				#Haplotype counts 
#				if j < i:
#					hap_set, hap_counts, haplotypes = snpsdata.get_haplotypes([snp1, snp2], self.n,
#												count_haplotypes=True)
#					d = {'num_haps':len(hap_counts)}
#					for hap, hap_c in zip(hap_set, hap_counts):
#						d[hap] = hap_c
#					haplotype_counts[i][j] = d
#
#				if sp.any(snp1 != snp2) and sp.any(anti_snp1 != snp2):
#					snp_mat = H_sqrt_inv * (sp.matrix(snp2).T) #Transformed inputs
#					X = sp.hstack([h0_X1, snp_mat])
#					(betas, rss, rank, s) = linalg.lstsq(X, Y)
#					if not rss:
#						if verbose: print 'No predictability in the marker, moving on...'
#						continue
#					f_stat = ((h0_rss - rss) / (q)) / (rss / n_p)
#					p_val = stats.f.sf(f_stat, q, n_p)
#					p_val_array[i, j] = p_val[0]
#					f_stat_array[i, j] = f_stat[0]
#					rss_array[i, j] = rss[0]
#					var_perc_array[i, j] = float(1 - rss / h0_rss)
#				else:
#					rss_array[i, j] = h0_rss[0]
#
#		for i in range(len(snps)):
#			for j in range(len(snps[:i])):
#				if p_val_array[i, j] > p_val_array[j, i]:
#					p_val_array[j, i] = p_val_array[i, j]
#					f_stat_array[j, i] = f_stat_array[i, j]
#					rss_array[j, i] = rss_array[i, j]
#					var_perc_array[j, i] = var_perc_array[i, j]
#
#
#		no_snp_count = 0
#		identical_snp_count = 0
#		num_same_snp = 0
#
#		#Now the interaction.
#		p = len(self.X.T) + q + 2 #Adding two SNP as a cofactor.
#		n_p = self.n - p
#		for i, snp1 in enumerate(snps):
#			for j, snp2 in enumerate(snps[:i]):
#				#t_X = sp.mat(sp.vstack([sp.ones(self.n), snp1, snp2])).T
#				#(t_beta, t_rss, t_rank, t_sigma) = linalg.lstsq(t_X, sp.matrix(pseudo_snp).T)
#				failed = False
#				pseudo_snp = snp1 * snp2 #Interaction
#				if haplotype_counts[i][j]['num_haps'] == 4:
#					#DO ANOVA... to simplify beta interpretation.
#					#if t_rss and t_rss[0] > 1e-20 :
#					snp_mat = H_sqrt_inv * sp.mat(sp.vstack([snp1, snp2])).T #Transformed inputs
#					h0_X1 = sp.hstack([h0_X, snp_mat])
#					snp_mat = H_sqrt_inv * (sp.matrix(pseudo_snp).T) #Transformed inputs
#					X = sp.hstack([h0_X1, snp_mat])
#					(betas, rss, rank, s) = linalg.lstsq(X, Y)
#					if rss:
#						h0_rss = rss_array[j, i]
#						f_stat = ((h0_rss - rss) / (q)) / (rss / n_p)
#						p_val = stats.f.sf(f_stat, q, n_p)
#						p_val_array[i, j] = p_val[0]
#						f_stat_array[i, j] = f_stat[0]
#						rss_array[i, j] = rss[0]
#						var_perc_array[i, j] = float(1 - rss / h0_rss)
#
#					else:
#						failed = True
#
#				else:
#					failed = True
#
#				if failed:
#					anti_snp1 = get_anti_snp(snp1)
#					anti_snp2 = get_anti_snp(snp2)
#					if sp.all(snp1 == snp2) or sp.all(snp1 == anti_snp2):
#						num_same_snp += 1
#					elif len(sp.unique(pseudo_snp)) == 1:
#						no_snp_count += 1
#					elif sp.any(pseudo_snp != snp1) or sp.any(pseudo_snp != snp2) or \
#						sp.any(pseudo_snp != anti_snp1) or sp.any(pseudo_snp != anti_snp2):
#						identical_snp_count += 1
#					p_val_array[i, j] = 1
#					f_stat_array[i, j] = 0
#					rss_array[i, j] = full_rss
#					var_perc_array[i, j] = 0
#
#
#		tot_num = num_snps * (num_snps - 1) / 2.0
#		print 'fail_fraction =%f, no_snp_fraction=%f, identical_pseudo_snp_fraction=%f, num_same_snp=%f'\
#			% ((no_snp_count + identical_snp_count + num_same_snp) / tot_num, \
#			no_snp_count / tot_num, identical_snp_count / tot_num, num_same_snp / tot_num)
#
#		res_dict = {'vincent_ps':p_val_array, 'vincent_f_stats':f_stat_array, 'vincent_rss':rss_array, \
#			    'vincent_var_perc':var_perc_array, 'pseudo_heritability': pseudo_her, \
#			    'haplotype_counts':haplotype_counts}
#
#		#Filling up the diagonal
#		f_stat_array = sp.zeros((num_snps, num_snps))
#		rss_array = sp.repeat(full_rss, num_snps * num_snps)
#		rss_array.shape = (num_snps, num_snps)
#		p_val_array = sp.ones((num_snps, num_snps))
#		var_perc_array = sp.zeros((num_snps, num_snps))
#		for i in range(num_snps):
#			p_val_array[i, i] = res_dict['vincent_ps'][i, i]
#			f_stat_array[i, i] = res_dict['vincent_f_stats'][i, i]
#			rss_array[i, i] = res_dict['vincent_rss'][i, i]
#			var_perc_array[i, i] = res_dict['vincent_var_perc'][i, i]
#
#		q1 = 2 #Testing significance of two SNPs, compared to none
#		p1 = len(self.X.T) + q1
#		n_p1 = self.n - p1
#		q2 = 3 #Testing significance of two SNPs and the interaction, compared to none
#		p2 = len(self.X.T) + q2
#		n_p2 = self.n - p2
#		for i, snp1 in enumerate(snps):
#			for j, snp2 in enumerate(snps[:i]):
#				anti_snp1 = get_anti_snp(snp1)
#				anti_snp2 = get_anti_snp(snp2)
#				pseudo_snp = snp1 * snp2
#				if sp.any(snp1 != snp2) and sp.any(snp1 != anti_snp2):
#					snp_mat = H_sqrt_inv * sp.mat(sp.vstack([snp1, snp2])).T #Transformed inputs
#					X = sp.hstack([h0_X, snp_mat])
#					(betas, rss, rank, s) = linalg.lstsq(X, Y)
#					f_stat = ((full_rss - rss) / (q1)) / (rss / n_p1)
#					p_val = stats.f.sf(f_stat, q1, n_p1)
#					p_val_array[j, i] = p_val[0]
#					f_stat_array[j, i] = f_stat[0]
#					rss_array[j, i] = rss[0]
#					var_perc_array[j, i] = float(1 - rss / full_rss)
#				else:
#					rss = res_dict['vincent_rss'][j, i]
#					f_stat = ((full_rss - rss) / (q1)) / (rss / n_p1)
#					p_val = stats.f.sf(f_stat, q1, n_p1)
#					p_val_array[j, i] = p_val[0]
#					f_stat_array[j, i] = f_stat[0]
#					rss_array[j, i] = rss
#					var_perc_array[j, i] = float(1 - rss / full_rss)
#
#				t_X = sp.mat(sp.vstack([sp.ones(self.n), snp1, snp2])).T
#				(t_beta, t_rss, t_rank, t_sigma) = linalg.lstsq(t_X, sp.matrix(pseudo_snp).T)
#				if t_rss and t_rss[0] > 1e-20 :
#					snp_mat = H_sqrt_inv * sp.mat(sp.vstack([snp1, snp2])).T  #Transformed inputs
#					X = sp.hstack([h0_X, snp_mat])
#					(betas, rss, rank, s) = linalg.lstsq(X, Y)
#					f_stat = ((full_rss - rss) / (q2)) / (rss / n_p2)
#					p_val = stats.f.sf(f_stat, q2, n_p2)
#					p_val_array[i, j] = p_val[0]
#					f_stat_array[i, j] = f_stat[0]
#					rss_array[i, j] = rss[0]
#					var_perc_array[i, j] = float(1 - rss / full_rss)
#				else:
#					rss = rss_array[j, i]
#					f_stat = ((full_rss - rss) / (q2)) / (rss / n_p2)
#					p_val = stats.f.sf(f_stat, q2, n_p2)
#					p_val_array[i, j] = p_val[0]
#					f_stat_array[i, j] = f_stat[0]
#					rss_array[i, j] = rss
#					var_perc_array[i, j] = var_perc_array[j, i]
#
#		res_dict['ps'] = p_val_array
#		res_dict['f_stats'] = f_stat_array
#		res_dict['rss'] = rss_array
#		res_dict['var_perc'] = var_perc_array
#		return res_dict



def get_emma_reml_estimates(y, K):
	lmm = LinearMixedModel(y)
	lmm.add_random_effect(K)
	res = lmm.get_REML()
	res['Y_t'] = res['H_sqrt_inv'] * lmm.Y	#The transformed outputs.
	res['X_t'] = res['H_sqrt_inv'] * lmm.X
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



def emmax(snps, phenotypes, K, cofactors=None, Z=None, with_betas=False):
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

	print "Running EMMAX"
	s1 = time.time()
	res = lmm.emmax_f_test(snps, Z=Z, with_betas=with_betas)
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



def _get_interactions_(isnp, snps, mac_threshold=15):
	isnps = []
	cofactor_indices = []
	anti_isnp = sp.vectorize(lambda x: 1 if x == 0 else 0)(isnp)
	for i, snp in enumerate(snps):
		min_count = min(min(sp.bincount(isnp & snp)), min(sp.bincount(isnp | snp)), \
			min(sp.bincount(snp & anti_isnp)), min(sp.bincount(snp | anti_isnp)))
		if min_count > mac_threshold and min:
			isnps.append(isnp & snp)
			cofactor_indices.append(i)
	return isnps, cofactor_indices



def _log_choose_(n, k):
 	#Srinivasa Ramanujan approximation of log(n!)
 	#log_fact = lambda n: n * sp.log(n) - n + (sp.log(n * (1 + 4 * n * (1 + 2 * n))) / 6) + sp.log(sp.pi) / 2
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



def emmax_step_wise(phenotypes, K, sd=None, all_snps=None, all_positions=None,
		all_chromosomes=None, num_steps=10, file_prefix=None, allow_interactions=False,
		interaction_pval_thres=0.01, forward_backwards=True):
	"""
	Run EMMAX stepwise.. forward, with one possible backward at each step.
	"""
	import plotResults as pr

	def _to_string_(cofactors):
		st = ''
		if len(cofactors) > 0:
			for tup in cofactors[:-1]:
				st += '%d_%d_%f,' % (tup[0], tup[1], -sp.log10(tup[2]))
			st += '%d_%d_%f' % (cofactors[-1][0], cofactors[-1][1], -sp.log10(cofactors[-1][2]))
		return st

	if sd:
	 	all_snps = sd.getSnps()
	 	all_positions = sd.getPositions()
	 	all_chromosomes = sd.get_chr_list()

	snps = all_snps[:]
	positions = all_positions[:]
	chromosomes = all_chromosomes[:]
	chr_pos_list = zip(chromosomes, positions)
       	lmm = LinearMixedModel(phenotypes)
	lmm.add_random_effect(K)

	print "Running EMMAX stepwise"
	s1 = time.time()
 	step_info_list = []
	cofactors = []  #A list of the loci found, together with their statistics.
	cofactor_snps = []
	interactions = []
 	step_i = 0
 	num_par = 2 #mean and variance scalar

	reml_res = lmm.get_REML()
	ml_res = lmm.get_ML()
	H_sqrt_inv = reml_res['H_sqrt_inv']
	ll = ml_res['max_ll']
	rss = float(reml_res['rss'])
	reml_mahalanobis_rss = float(reml_res['mahalanobis_rss'])
	criterias = {'ebics':[], 'mbics':[]}
	(bic, extended_bic, modified_bic) = _calc_bic_(ll, len(snps), num_par, lmm.n) #Calculate the BICs
	criterias['ebics'].append(extended_bic)
	criterias['mbics'].append(modified_bic)
	action = 'None'
	print '\nStep %d: action=%s, num_par=%d, p_her=%0.4f, ll=%0.2f, rss=%0.2f, reml_m_rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f' % \
		(step_i, action, num_par, reml_res['pseudo_heritability'], ll, rss, reml_mahalanobis_rss, \
		 bic, extended_bic, modified_bic)
	print 'Cofactors:', _to_string_(cofactors)

	for step_i in range(1, num_steps + 1):
		emmax_res = lmm._emmax_f_test_(snps, H_sqrt_inv)
		min_pval_i = sp.argmin(emmax_res['ps'])
		min_pval = emmax_res['ps'][min_pval_i]
		mahalnobis_rss = emmax_res['rss'][min_pval_i]
		min_pval_chr_pos = chr_pos_list[min_pval_i]
		print 'Min p-value:', min_pval
		print 'Min Mahalanobis RSS:', mahalnobis_rss
		step_info = {'pseudo_heritability':reml_res['pseudo_heritability'], 'rss':rss, \
			'reml_mahalanobis_rss': reml_res['mahalanobis_rss'], 'mahalanobis_rss':mahalnobis_rss,
			'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic, 'ps': emmax_res['ps'],
			'cofactors':map(tuple, cofactors[:]), 'cofactor_snps':cofactor_snps[:], 'min_pval':min_pval,
			'min_pval_chr_pos': min_pval_chr_pos, 'interactions':interactions}
		step_info_list.append(step_info)

		#Plot gwas results per step 
		if file_prefix:
			pdf_file_name = file_prefix + '_step' + str(step_i - 1) + '.pdf'
			png_file_name = file_prefix + '_step' + str(step_i - 1) + '.png'
			#Manhattan plot
			pr.plot_raw_result(emmax_res['ps'], chromosomes, positions, highlight_markers=cofactors,
					 png_file=png_file_name)
			#QQ plot



		#Adding the new SNP as a cofactor
		lmm.add_factor(snps[min_pval_i])
		cofactor_snps.append(snps[min_pval_i])
		reml_res = lmm.get_REML()
		ml_res = lmm.get_ML()
		H_sqrt_inv = reml_res['H_sqrt_inv']
		ll = ml_res['max_ll']
		rss = float(reml_res['rss'])
		reml_mahalanobis_rss = float(reml_res['mahalanobis_rss'])
		num_par += 1
		action = '+'

		cofactors.append([min_pval_chr_pos[0], min_pval_chr_pos[1], min_pval])


		#Re-estimate the p-value of the cofactors... with the smallest in the list.
		for i, snp in enumerate(cofactor_snps):
			t_cofactors = cofactor_snps[:]
			del t_cofactors[i]
			lmm.set_factors(t_cofactors)
			cofactors[i][2] = lmm._emmax_f_test_([snp], H_sqrt_inv)['ps'][0]
		lmm.set_factors(cofactor_snps)


		#Remove the found SNP from considered SNPs
		del snps[min_pval_i]
		del positions[min_pval_i]
		del chromosomes[min_pval_i]
		del chr_pos_list[min_pval_i]

		#Try adding an interaction.... 
#		if allow_interactions and len(cofactor_snps) > 1:
#			isnps, cofactor_indices = _get_interactions_(cofactor_snps[-1], cofactor_snps[:-1])
#			if isnps:
#				emmax_res = lmm._emmax_f_test_(isnps, H_sqrt_inv)
#				min_pval_i = sp.argmin(emmax_res['ps'])
#				print emmax_res['ps'][min_pval_i], emmax_res['rss'][min_pval_i]
#				if emmax_res['ps'][min_pval_i] < interaction_pval_thres and emmax_res['rss'][min_pval_i] < rss:
#					if lmm.add_factor(snps[min_pval_i]):
#						cofactor_snps.append(isnps[min_pval_i])
#						interactions.append(((cofactors[min_pval_i][0], cofactors[min_pval_i][1]),
#								(cofactors[-1][0], cofactors[-1][1]), emmax_res['ps'][min_pval_i]))
#						cofactors_str += '%d_%d_X_%d_%d,' \
#								% (cofactors[min_pval_i][0], cofactors[min_pval_i][1],
#								cofactors[-1][0], cofactors[-1][1])
#						step_info['interactions'] = interactions
#						action += 'i'
#						reml_res = lmm.get_REML()
#						ml_res = lmm.get_ML()
#						H_sqrt_inv = reml_res['H_sqrt_inv']
#						ll = reml_res['max_ll']
#						rss = float(reml_res['rss'])
#						reml_mahalanobis_rss = float(reml_res['mahalanobis_rss'])
#						num_par += 1
#						print "Just added an interaction:", interactions


		(bic, extended_bic, modified_bic) = _calc_bic_(ll, len(snps), num_par, lmm.n) #Calculate the BICs
		criterias['ebics'].append(extended_bic)
		criterias['mbics'].append(modified_bic)

		print '\nStep %d: action=%s, num_par=%d, p_her=%0.4f, ll=%0.2f, rss=%0.2f, reml_m_rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f' % \
			(step_i, action, num_par, reml_res['pseudo_heritability'], ll, rss, reml_mahalanobis_rss, \
			bic, extended_bic, modified_bic)
		print 'Cofactors:', _to_string_(cofactors)
		if reml_res['pseudo_heritability'] < 0.01:
			print 'Breaking early, since pseudoheritability is close to 0.'
			break

	emmax_res = lmm._emmax_f_test_(snps, H_sqrt_inv)
	min_pval_i = sp.argmin(emmax_res['ps'])
	min_pval = emmax_res['ps'][min_pval_i]
	mahalnobis_rss = emmax_res['rss'][min_pval_i]
	min_pval_chr_pos = chr_pos_list[min_pval_i]
	print 'Min p-value:', min_pval
	print 'Min Mahalanobis RSS:', mahalnobis_rss
	step_info = {'pseudo_heritability':reml_res['pseudo_heritability'], 'rss':rss, \
		'reml_mahalanobis_rss': reml_res['mahalanobis_rss'], 'mahalanobis_rss':mahalnobis_rss,
		'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic, 'ps': emmax_res['ps'],
		'cofactors':map(tuple, cofactors[:]), 'cofactor_snps':cofactor_snps[:], 'min_pval':min_pval,
		'min_pval_chr_pos': min_pval_chr_pos, 'interactions':interactions}
	step_info_list.append(step_info)

	#Now plotting!
	print "Generating plots"
	if file_prefix:
		pdf_file_name = file_prefix + '_step' + str(step_i) + '.pdf'
		png_file_name = file_prefix + '_step' + str(step_i) + '.png'
		pr.plot_raw_result(emmax_res['ps'], chromosomes, positions, highlight_markers=cofactors,
				png_file=png_file_name)

	#Now backward stepwise.
	if forward_backwards:
		print 'Starting backwards..'
		while len(cofactor_snps) > 0:
			step_i += 1
			f_stats = sp.zeros(len(cofactor_snps))
			for i, snp in enumerate(cofactor_snps):
				t_cofactors = cofactor_snps[:]
				del t_cofactors[i]
				lmm.set_factors(t_cofactors)
				res = lmm._emmax_f_test_([snp], H_sqrt_inv)
				cofactors[i][2] = res['ps'][0]
				f_stats[i] = res['f_stats'][0]
			i_to_remove = f_stats.argmin()
			del cofactor_snps[i_to_remove]
			del cofactors[i_to_remove]
			lmm.set_factors(cofactor_snps)


			#Re-estimating the REML and ML.
			reml_res = lmm.get_REML()
			ml_res = lmm.get_ML()
			ll = ml_res['max_ll']
			rss = float(reml_res['rss'])
			reml_mahalanobis_rss = float(reml_res['mahalanobis_rss'])
			num_par -= 1
			action = '-'

			#Update the p-values
			for i, snp in enumerate(cofactor_snps):
				t_cofactors = cofactor_snps[:]
				del t_cofactors[i]
				lmm.set_factors(t_cofactors)
				res = lmm._emmax_f_test_([snp], H_sqrt_inv)
				cofactors[i][2] = res['ps'][0]

			#Calculate the BICs
			(bic, extended_bic, modified_bic) = _calc_bic_(ll, len(snps), num_par, lmm.n)
			criterias['ebics'].append(extended_bic)
			criterias['mbics'].append(modified_bic)
			print '\nStep %d: action=%s, num_par=%d, p_her=%0.4f, ll=%0.2f, rss=%0.2f, reml_m_rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f' % \
				(step_i, action, num_par, reml_res['pseudo_heritability'], ll, rss,
				reml_mahalanobis_rss, bic, extended_bic, modified_bic)
			print 'Cofactors:', _to_string_(cofactors)

			step_info = {'pseudo_heritability':reml_res['pseudo_heritability'], 'rss':rss, \
				'reml_mahalanobis_rss': reml_res['mahalanobis_rss'], 'll':ll, 'bic':bic,
				'e_bic':extended_bic, 'm_bic':modified_bic, 'cofactors':map(tuple, cofactors[:]),
				'cofactor_snps':cofactor_snps[:], 'mahalanobis_rss':None, 'min_pval':None,
				'min_pval_chr_pos':None}
			step_info_list.append(step_info)
			print cofactors


	for c in criterias:
		print 'GWAs for optimal %s criteria:' % c
		i_opt = sp.array(criterias[c]).argmin()
		print "    %d'th step was optimal." % i_opt
		cofactor_snps = step_info_list[i_opt]['cofactor_snps']
		cofactors = step_info_list[i_opt]['cofactors']
		print cofactors
		lmm.set_factors(cofactor_snps)
		reml_res = lmm.get_REML()
		H_sqrt_inv = reml_res['H_sqrt_inv']
		emmax_res = lmm._emmax_f_test_(all_snps, H_sqrt_inv)
		min_pval_i = emmax_res['ps'].argmin()
		min_pval = emmax_res['ps'][min_pval_i]
		mahalnobis_rss = emmax_res['rss'][min_pval_i]
		min_pval_chr_pos = chr_pos_list[min_pval_i]
		print 'Min p-value:', min_pval
		print 'Min Mahalanobis RSS:', mahalnobis_rss

		if file_prefix:
			pdf_file_name = file_prefix + '_step' + str(i_opt) + '_opt_' + c + '.pdf'
			png_file_name = file_prefix + '_step' + str(i_opt) + '_opt_' + c + '.png'
			pr.plot_raw_result(emmax_res['ps'], all_chromosomes, all_positions, highlight_markers=cofactors,
					 png_file=png_file_name)




	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds.' % (mins, secs)
	else:
		print 'Took %f seconds.' % (secs)


	if file_prefix:
		p_her_list = []
		rss_list = []
		reml_mahalanobis_rss_list = []
		mahalanobis_rss_list = []
		ll_list = []
		bic_list = []
		e_bic_list = []
		m_bic_list = []
		min_pval_list = []
		f = open(file_prefix + "_stats.csv", 'w')
		d_keys = ['pseudo_heritability', 'rss', 'reml_mahalanobis_rss', 'mahalanobis_rss', 'll', 'bic', 'e_bic', 'm_bic', 'min_pval']
		f.write(','.join(['step_nr'] + d_keys + ['min_pval_pos_chr', 'cofactors']) + '\n')
		for i, si in enumerate(step_info_list):
			st = ','.join(map(str, [i] + [si[k] for k in d_keys]))
			if si['min_pval_chr_pos']:
				st += ',%d_%d,' % si['min_pval_chr_pos']
			else:
				st += ',,'
			st += _to_string_(si['cofactors'])
			st += '\n'
			f.write(st)
			p_her_list.append(float(si['pseudo_heritability']))
			rss_list.append(float(si['rss']))
			reml_mahalanobis_rss_list.append(float(si['reml_mahalanobis_rss']))
			if si['mahalanobis_rss']:
				mahalanobis_rss_list.append(float(si['mahalanobis_rss']))
			ll_list.append(si['ll'])
			bic_list.append(si['bic'])
			e_bic_list.append(si['e_bic'])
			m_bic_list.append(si['m_bic'])
			if si['min_pval']:
				min_pval_list.append(float(si['min_pval']))
		f.close()
		import pylab
		pylab.figure(figsize=(6, 4))
		#pylab.axes([0.05, 0.05, 0.92, 0.95])
		pylab.plot(range(len(p_her_list)), p_her_list, 'o-')
		pylab.ylabel('Pseudo-heritability')
		pylab.savefig(file_prefix + '_stats_p_her.pdf', format='pdf')
		pylab.clf()
		pylab.plot(range(len(rss_list)), rss_list, 'o-')
		pylab.ylabel('RSS')
		pylab.savefig(file_prefix + '_stats_rss.pdf', format='pdf')
		pylab.clf()
		pylab.plot(range(len(reml_mahalanobis_rss_list)), reml_mahalanobis_rss_list, 'o-')
		pylab.ylabel('REML Mahalanobis RSS')
		pylab.savefig(file_prefix + '_stats_reml_mahalanobis_rss.pdf', format='pdf')
		pylab.clf()
		pylab.plot(range(len(mahalanobis_rss_list)), mahalanobis_rss_list, 'o-')
		pylab.ylabel('Mahalanobis RSS')
		pylab.savefig(file_prefix + '_stats_mahalanobis_rss.pdf', format='pdf')
		pylab.clf()
		pylab.plot(range(len(ll_list)), ll_list, 'o-')
		pylab.ylabel('Log likelihood')
		pylab.savefig(file_prefix + '_stats_ll.pdf', format='pdf')
		pylab.clf()
		pylab.plot(range(len(min_pval_list)), map(lambda x:-sp.log10(x), min_pval_list), 'o-')
		pylab.ylabel('Min. p-value')
		pylab.savefig(file_prefix + '_stats_pval.pdf', format='pdf')
		pylab.clf()
		pylab.plot(range(len(bic_list)), bic_list, 'o-')
		pylab.ylabel('BIC')
		pylab.savefig(file_prefix + '_stats_bic.pdf', format='pdf')
		pylab.clf()
		pylab.plot(range(len(e_bic_list)), e_bic_list, 'o-')
		pylab.ylabel('Extended BIC')
		pylab.savefig(file_prefix + '_stats_ebic.pdf', format='pdf')
		pylab.clf()
		pylab.plot(range(len(m_bic_list)), m_bic_list, 'o-')
		pylab.ylabel('Modified BIC')
		pylab.savefig(file_prefix + '_stats_mbic.pdf', format='pdf')
		pylab.clf()
		max_rss = max(rss_list)
		rss_array = sp.array(rss_list) / max_rss
		p_her_array = rss_array * sp.array(p_her_list)
		genetic_variance = p_her_array + (1 - rss_array)
		variance_explained = (1 - rss_array)
		pylab.figure(figsize=(10, 6))
		pylab.fill_between([0, step_i], 0, 1, color='#DD3333', alpha=0.8, label='Variance explained')
		pylab.fill_between(sp.arange(step_i + 1), 0, genetic_variance, color='#22CC44', alpha=0.8, label='Genetic variance')
		pylab.fill_between(sp.arange(step_i + 1), 0, variance_explained, color='#2255AA', alpha=0.8, label='Variance explained')
		pylab.ylabel('Percentage of variance')
		pylab.xlabel('Step number')
		pylab.axvline(x=num_steps, c='k', linestyle=':')
		pylab.legend(loc=1, ncol=3, shadow=True)

		pylab.savefig(file_prefix + '_stats_variances.png', format='png')
		pylab.savefig(file_prefix + '_stats_variances.pdf', format='pdf')

	return step_info_list



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



def _filter_k_(k, indices_to_keep):
	new_k = sp.zeros((len(indices_to_keep), len(indices_to_keep)))
	for i in range(len(indices_to_keep)):
		for j in range(len(indices_to_keep)):
			new_k[i, j] = k[indices_to_keep[i], indices_to_keep[j]]
	k = new_k
	return k


def filter_k_for_accessions(k, k_accessions, accessions):
	indices_to_keep = []
	for i, acc in enumerate(k_accessions):
		if acc in accessions:
			indices_to_keep.append(i)
	return _filter_k_(k, indices_to_keep)

def load_kinship_from_file(kinship_file, accessions=None):
	assert os.path.isfile(kinship_file), 'File not found.'
	#sys.stdout.write("Loading K.\n")
	#sys.stdout.flush()
	f = open(kinship_file, 'r')
	l = cPickle.load(f)
	f.close()
	k = l[0]
	k_accessions = l[1]
	return filter_k_for_accessions(k, k_accessions, accessions)


def _test_stepwise_emmax_():
	import dataParsers as dp
	import phenotypeData as pd
	#filename = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_expression_101011.txt"
	#filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_raw_112210.csv"
	filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/b_dilkes_metabolites.csv"
	mac_threshold = 15
	for pid, trans in [(0, 'log_trans')]:#[ (5, False)]:#, (2, False)]:#(617, False), (619, False), (621, False), (623, False), (625, False), (627, False), (629, False), (631, False), (633, False)]:#, (5, False), (226, True), (264, False), (1025, False)]:#[(264, False), (265, False), (266, False), (267, False)]:
		phed = pd.parse_phenotype_file(filename)
		phed.convert_to_averages()
		if trans == 'sqrt_trans':
			phed.sqrt_transform(pid)
		elif trans == 'log_trans':
			phed.log_transform(pid)
		phen_name = phed.get_name(pid)
		sd = dp.parse_numerical_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t72.csv.binary', filter=1)
		sd.coordinate_w_phenotype_data(phed, pid)
		if mac_threshold:
			sd.filter_mac_snps(mac_threshold) #Filter MAF SNPs!
		phenotypes = phed.get_values(pid)
		K = load_kinship_from_file('/Users/bjarnivilhjalmsson/Projects/Data/250k/kinship_matrix_cm72.pickled',
					phed.get_ecotypes(pid))

		info_list = emmax_step_wise(phenotypes, K, sd=sd, \
					file_prefix='/Users/bjarni.vilhjalmsson/tmp/emmax_stepwise_dilkes' \
					+ str(pid) + '_' + phen_name, num_steps=60, allow_interactions=True,
					interaction_pval_thres=0.001)



def _test_joint_analysis_():
	import dataParsers as dp
	import phenotypeData as pd
	filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_raw_112210.csv"
	phed = pd.parse_phenotype_file(filename)
	phed.convert_to_averages()

	mac_threshold = 15
	pids = []

	join_phen_vals = []
	joint_ecotypes = []
	for pid in pids:
		sd = dp.parse_numerical_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t72.csv.binary', filter=1)
		sd.coordinate_w_phenotype_data(phed, pid)
		if mac_threshold:
			sd.filter_mac_snps(mac_threshold) #Filter MAF SNPs!
		snps = sd.getSnps()
		phenotypes = phed.get_values(pid)
		K = load_kinship_from_file('/Users/bjarnivilhjalmsson/Projects/Data/250k/kinship_matrix_cm72.pickled',
					phed.get_ecotypes(pid))
		emma(snps, phenotypes, K)



	#1. Get pseudoheritabilities for all phenotypes

	#2. Get correlations between traits.. (if correlation method is used).

	#3. Figure out all ecotypes which are in 250K and for which we have phenotypic values.

	#4. Create a new phenotype, and SNP indicator vector.

	#5. Run analysis..







if __name__ == "__main__":
	_test_stepwise_emmax_()
