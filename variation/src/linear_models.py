"""
Contains functions to perform various linear regression schemes, such as simple, and mixed models.
"""
#import numpy as np
from scipy import *
from scipy import linalg
from scipy import stats
from scipy import optimize
import math
import warnings
import pdb
import cPickle
import os
import sys
import multiprocessing as mp
import time


class LinearModel(object):
        """
        A simple linear model
        """
        def __init__(self,Y=None):
                """
                The fixed effects should be a list of fixed effect lists (SNPs)
                """
                self.n = len(Y)
                self.Y = matrix(Y).T
                self.X = matrix([1]*len(Y)).T #The intercept
                self.p = 1
                self.beta_est = None
                self.cofactors = []


        def add_factor(self, x, lin_depend_thres = 1e-8):
                """
                Adds an explanatory variable to the X matrix.
                """
                #Checking whether this new cofactor in linearly independent.
                new_x = matrix(x).T
                (beta,rss,rank,sigma) = linalg.lstsq(self.X,new_x)
                if float(rss)<lin_depend_thres:
                        warnings.warn('A factor was found to be linearly dependent on the factors already in the X matrix.  Hence skipping it!')
                        return
                self.X = hstack([self.X,new_x])
                self.cofactors.append(x)
                self.p += 1
        
        
        def get_hat_matrix(self):
                self.X_squared_inverse = (self.X.T*self.X).I
                self.hat_matrix = self.X*self.X_squared_inverse*self.X.T
                return self.hat_matrix
                
        
        
        def least_square_estimate(self): 
                """
                Via Normal equations, get LSEs
                """
                self.X_squared_inverse = (self.X.T*self.X).I
                self.beta_est = self.X_squared_inverse*self.X.T*self.Y
                return self.beta_est
        
        def get_residuals(self):
                """
                Calculates and returns the residuals as a column vector.
                """
                #if self.beta_est==None:
                self.least_square_estimate()
                self.residuals = self.Y - self.X*self.beta_est
                return self.residuals
        
        
        def general_f_test(self,A,c):
                """
                A general F-test implementation.
                Where the hypothesis is A*beta=c, a constraint.
                
                Here A is a matrix, and c a column vector                
                """
                #if not self.residuals:
                self.get_residuals()
                q,p = shape(A)
                assert p == self.p, 'Shape of A is wrong!'
                B = (A*self.beta_est-c)
                M = A*(self.X_squared_inverse)*A.T
                f_stat = (B.T*M.I*B)/((self.residuals.T*self.residuals)/(self.n-self.p))
                p_value = 1-stats.f.cdf(f_stat,q,self.n-self.p)
                return p_value,f_stat
        
        def fast_f_test(self,snps):
                """
                A standard linear model, using a F-test
                """                
                (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(self.X,self.Y)
                num_snps = len(snps)
                rss_list = repeat(h0_rss,num_snps)
                betas_list = [h0_betas]*num_snps
                var_perc = zeros(num_snps)
                q = 1
                n_p = self.n-self.p
                missing_indices = []
                for i, snp in enumerate(snps):
                        (betas, rss, rank, s) = linalg.lstsq(hstack([self.X,matrix(snp).T]),self.Y)
                        if not rss:
                                print 'No predictability in the marker, moving on...'
                                missing_indices.append(i)
                                continue                               
                        rss_list[i] = rss[0]                
                        betas_list[i] = map(float,list(betas))
                rss_ratio = h0_rss/rss_list
                var_perc = 1-1/rss_ratio
                f_stats = (rss_ratio-1)*n_p/float(q)
                p_vals = stats.f.sf(f_stats,q,n_p)
                return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list, 
                        'var_perc':var_perc, 'missing_indices':missing_indices}

                
        

        def test_explanatory_variable(self,x):
                """
                Returns a p-value for whether adding this variable to the model explains the model better.
                
                Hopefully a sped-up version of a specific F-test.
                """
                
                #THIS CAN BE SPED-UP MORE IF WE CHECK WHETHER self.X IS A VECTOR.  
                #AND USE t-TEST. 
                res_1 = self.get_residuals()
                
                X_mat = hstack([self.X,matrix([[v] for v in x])])
                X_sq = X_mat.T*X_mat
                try:
                        X_sq_inv = X_sq.I
                except Exception, err_str:
                        print err_str
                        raise Exception('Andskotinn!!')
                
                res_2 = self.Y - X_mat*X_sq_inv*X_mat.T*self.Y                
                rss_1 = res_1.T*res_1                
                rss_2 = res_2.T*res_2
                f_stat = (rss_1-rss_2)/(rss_2/(self.n-self.p+1))
                p_value = stats.f.sf(f_stat,1,self.n-self.p+1)
                return p_value,f_stat                
        
        


                
                
class LinearMixedModel(LinearModel):
        def __init__(self,Y=None):
                """
                The fixed effects should be a list of fixed effect lists (SNPs)
                """
                self.n = len(Y)
                self.y_var = var(Y,ddof=1)
                self.Y = matrix([[y] for y in Y])
                self.X = matrix([[1] for y in Y]) #The intercept
                self.p = 1
                self.beta_est = None
                self.cofactors = []
                
                #A list of random effect type, and the cov matrix.
                self.random_effects = [('normal',matrix(identity(self.n)))] #The first random effect is the IID error.
                        

        def add_random_effect(self,cov_matrix=None,effect_type='normal'):
                if effect_type!='normal':
                        raise Exception('Currently, only Normal random effects are allowed.')
                self.random_effects.append((effect_type,cov_matrix))

#        def least_square_estimate(self): 
#                raise Exception("LSE not applicable for mixed models")


        def _get_eigen_L_(self,K):
                evals,evecs = linalg.eigh(K)
                #FIXME, eigenvalues and vectors are in opposite order compared with R.        
                return {'values':evals,'vectors':mat(evecs).T}


        
        def _get_eigen_R_(self,X,hat_matrix=None,complete=True):
                q = X.shape[1]
                if not hat_matrix:
                        X_squared_inverse = linalg.pinv(X.T*X) #(X.T*X).I
                        hat_matrix = X*X_squared_inverse*X.T
                S = mat(identity(self.n))-hat_matrix        #S=I-X(X'X)^{-1}X'
                evals,evecs = linalg.eigh(S*(self.random_effects[1][1]+self.random_effects[0][1])*S) #eigen of S(K+I)S
                return {'values':map(lambda x: x-1, evals[q:]),'vectors':(mat(evecs).T[q:])}   #Because of S(K+I)S?



        def _ll_(self,delta,eig_vals,sq_etas):
                num_eig_vals = len(eig_vals)
                c_1 = 0.5*num_eig_vals*(math.log(num_eig_vals/(2.0*math.pi))-1)
                sum_1 = 0
                sum_2 = 0
                for j in range(num_eig_vals):
                        v_1 = eig_vals[j]+delta
                        v_2 = sq_etas[j]/v_1
                        sum_1 += v_2
                        sum_2 += math.log(v_1)
                return c_1-0.5*(num_eig_vals*math.log(sum_1)+sum_2)  #log-likelihoods (eq. 7 from paper)

 
        def _dll_(self,delta,eig_vals,sq_etas):
                num_eig_vals = len(eig_vals)
                sum_1 = 0
                sum_2 = 0
                sum_3 = 0
                for j in range(num_eig_vals):
                        v_1 = eig_vals[j]+delta
                        v_2 = sq_etas[j]/v_1
                        sum_1 += v_2
                        sum_2 += v_2/v_1
                        sum_3 += 1.0/v_1 
                return num_eig_vals*sum_2/sum_1-sum_3  #diffrentiated log-likelihoods (*2) (eq. 9 from paper)

        def get_expedited_REMLE(self, eig_L=None, xs=None, ngrids=50, llim=-4, ulim=10, esp=1e-6, 
                                return_pvalue=False, return_f_stat=False):
                """
                Get REML estimates for the effect sizes, as well as the random effect contributions.
                
                Using the EMMA algorithm.                
                """
                if not eig_L:
                        raise Exception
#                if xs:
#                        X = hstack([self.X,matrix([[v] for x in xs for v in x ])])
#                else:
#                        X = self.X
		X = self.X
                eig_R = self._get_eigen_R_(X)
                K = self.random_effects[1][1]
                t = K.shape[0] #number of rows
                q = X.shape[1] #number of columns
                print 'X, number of columns: %d, number of rows: %d'%(q,t)
                n = self.n
                p = n-q
                assert K.shape[0]==K.shape[1]==X.shape[0]==n,'Dimensions are wrong.'

                etas = eig_R['vectors']*self.Y
                sq_etas = [float(eta)*float(eta) for eta in etas]
                log_deltas =  [float(i)/ngrids*(ulim-llim)+llim  for i in range(ngrids+1)] #a list of deltas to search
                assert len(log_deltas)==ngrids+1,'Delta list size error.'
                deltas = map(math.exp,log_deltas)
                eig_vals = list(eig_R['values'])
                assert len(eig_vals)==p,'Number of eigenvalues is incorrect.'

                #LL <- 0.5*(p*(log((p)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
                
                c_1 = 0.5*p*(math.log((p)/(2.0*math.pi))-1)
                
               
                def calc_ll(d):
                        """
                        Calculates the likelihood, and the derivative given a delta.
                        """
                        sum_1 = 0
                        sum_2 = 0
                        sum_3 = 0
                        sum_4 = 0
                        for j in range(p):
                                v_1 = eig_vals[j]+d
                                v_2 = sq_etas[j]/v_1
                                sum_1 += v_2
                                sum_2 += math.log(v_1)
                                sum_3 += v_2/v_1
                                sum_4 += 1.0/v_1 
                        #LL <- 0.5*((p)*(log((p)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
                        ll = c_1-(0.5)*((p)*math.log(sum_1)+sum_2)  #log-likelihoods (eq. 7 from paper)
                        #dLL <- 0.5*((p)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
                        dll = 0.5*((p)*sum_3/sum_1-sum_4)  #diffrentiated log-likelihoods (eq. 9 from paper)
                        return ll,dll

                     
                last_ll,last_dll = calc_ll(deltas[0])
                max_ll = last_ll
                max_ll_i = 0
                zero_intervals = []
                lls = [max_ll]
                dlls = [last_dll] 
                for i, d in enumerate(deltas[1:]):
                        ll,dll = calc_ll(d)
                        lls.append(ll)
                        dlls.append(dll)
                        if ll>max_ll:
                                max_ll = ll
                                max_ll_i=i
                        if last_dll>0 and dll<0:
                                zero_intervals.append(((ll+last_ll)*0.5,i))
                        last_dll = dll

                if len(zero_intervals)>0: 
                        #print 'local maximum found:',zero_intervals
                        opt_ll,opt_i = max(zero_intervals)
                        opt_delta = 0.5*(deltas[opt_i]+deltas[opt_i+1])
                        #Newton-Raphson
                        new_opt_delta = optimize.newton(self._dll_, opt_delta, args=(eig_vals,sq_etas), tol=esp,maxiter=50)
                        if deltas[opt_i]-esp<new_opt_delta<deltas[opt_i+1]+esp:
                                opt_delta = new_opt_delta                               
                                opt_ll = self._ll_(opt_delta,eig_vals,sq_etas)
                        else:
                                opt_delta = new_opt_delta                               
                                opt_ll = self._ll_(opt_delta,eig_vals,sq_etas)
                                print 'Local maximum outside of suggested area??'
                                print opt_delta, new_opt_delta
                                raise Exception('Local maximum outside of suggested area??')
                        if opt_ll<max_ll:
                                opt_delta = deltas[max_ll_i]
                else:
                        opt_delta = deltas[max_ll_i]
                        opt_ll = max_ll

                        
                #eig_vals.reverse()
                l = map(lambda x,y: x/(y+opt_delta),sq_etas,eig_vals)
                opt_vg = sum(l)/(p)  #vg   
                opt_ve = opt_vg*opt_delta  #ve
                 
                H_sqrt = eig_L['vectors'].T*diag([math.sqrt(ev+opt_delta) for ev in eig_L['values']])
                H_inverse = eig_L['vectors'].T*diag([1.0/(ev+opt_delta) for ev in eig_L['values']])*eig_L['vectors']
                XX = X.T*(H_inverse*X)
                iXX = linalg.pinv(XX)
                beta_est = iXX*X.T*(H_inverse*self.Y)
                x_beta = X*beta_est
                residuals = self.Y-x_beta
                rss = residuals.T*H_inverse*residuals
                x_beta_var = var(x_beta,ddof=1)
                var_perc = x_beta_var/self.y_var
                res_dict = {'max_ll':opt_ll, 'delta':opt_delta, 'beta':beta_est, 've':opt_ve, 'vg':opt_vg,
                            'var_perc':var_perc, 'rss':rss, 'H_sqrt':H_sqrt}
                if xs and return_f_stat:
                        q=X.shape[1]-self.X.shape[1]
                        h0_iXX = linalg.pinv(self.X.T*(H_inverse*self.X))
                        h0_beta_est = h0_iXX*self.X.T*(H_inverse*self.Y)
                        h0_residuals = self.Y-self.X*h0_beta_est
                        h0_rss = h0_residuals.T*H_inverse*h0_residuals
                        f_stat = ((h0_rss-rss)/(q))/(rss/p)
                        
                        #Redo this t-statistic.
#                        t_stat = self.beta_est[q-1]/math.sqrt(float(iXX[q-1,q-1])*opt_vg)
                        res_dict['f_stat']=float(f_stat)
                if return_pvalue:
                        #p_val = stats.t.sf(abs(t_stat), p)*2
                        p_val = stats.f.sf(f_stat,q,p)
                        res_dict['p_val']=float(p_val)
                
                return res_dict


        
        def expedited_REML_t_test(self, snps, ngrids=50, llim=-4, ulim=10, esp=1e-6):
                """
                Single SNP analysis (Not as fast as the R-version?)
                """
                assert len(self.random_effects)==2,"Expedited REMLE only works when we have exactly two random effects."
                K = self.random_effects[1][1]
                eig_L = self._get_eigen_L_(K)
                f_stats = []
                vgs = []
                ves = []
                max_lls = []
                var_perc = []
                betas = []
                p_vals = []
                               
                for snp in snps:
                        res = self.get_expedited_REMLE(eig_L=eig_L, xs=[snp], ngrids=ngrids, llim=llim, ulim=ulim, 
                                                       esp=esp, return_pvalue=True, return_f_stat=True)
                        f_stats.append(res['f_stat'])                
                        vgs.append(res['vg'])                
                        ves.append(res['ve'])                
                        max_lls.append(res['max_ll'])                
                        var_perc.append(res['var_perc'])                
                        betas.append(map(float,list(res['beta'])))    
                        p_vals.append(res['p_val'])            
                        
                return {'ps':p_vals,'f_stats':f_stats,'vgs':vgs,'ves':ves,'var_perc':var_perc,
                        'max_lls':max_lls,'betas':betas}
        

        def emmax_f_test_w_interactions(self,snps,int_af_threshold=15):
                """
                EMMAX implementation (in python)
                Single SNPs
                
                With interactions between SNP and possible cofactors.
                """
                assert len(self.random_effects)==2,"Expedited REMLE only works when we have exactly two random effects."
                p_0 = len(self.X.T)  
                n = self.n
                                

                K = self.random_effects[1][1]
                eig_L = self._get_eigen_L_(K)
                res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
                delta = res['delta']
                print 'pseudo_heritability:',1.0/(1+delta)
                H_sqr = res['H_sqrt']
                H_sqrt_inv = H_sqr.I
                Y = H_sqrt_inv*self.Y        #The transformed outputs.
                h0_X = H_sqrt_inv*self.X
                (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X,Y)
                h0_betas = map(float,list(h0_betas))
                f_stats = []
                rss_list = []
                betas_list = []
                p_vals = []
                var_perc = []
                cofactors = array(self.cofactors)
                num_cof = len(cofactors)
                low_int_af_count = 0
                for snp in snps:
                        snp_a = array(snp)
                        if sum(snp_a)>int_af_threshold:
                                interactions = cofactors*snp_a
                                interactions = interactions[sum(interactions,axis=1)>int_af_threshold]
                                low_int_af_count += num_cof-len(interactions)                                        
                                X = hstack([h0_X,H_sqrt_inv*matrix(snp_a).T,H_sqrt_inv*matrix(interactions).T])
                        else:
                                low_int_af_count += num_cof
                                X = hstack([h0_X,H_sqrt_inv*matrix(snp_a).T]) 
                        (betas, rss, p, sigma) = linalg.lstsq(X,Y)
                        q = p - p_0
                        n_p = n-p
                        if not rss:
                                if q==0:
                                        print 'No predictability in the marker, moving on...'
                                        p_vals.append(1)            
                                        f_stats.append(0)                
                                        rss_list.append(h0_rss)                
                                        betas_list.append(h0_betas)
                                        var_perc.append(0)
                                        continue
                                else:
                                        res = (Y-X*betas)
                                        rss = res.T*res
                        f_stat = ((h0_rss-rss)/q)/(rss/n_p)
                        p_val = stats.f.sf(f_stat,q,n_p)
                        p_vals.append(p_val[0])            
                        f_stats.append(f_stat[0])                
                        rss_list.append(rss[0])                
                        betas_list.append(map(float,list(betas)))
                        var_perc.append(float(1-rss/h0_rss))
                
                print 'low_int_af_count:',low_int_af_count
                      
                return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list, 
                        'delta':delta, 'pseudo_heritability': 1.0/(1+delta), 'var_perc':var_perc}


        
        def emmax_f_test(self,snps):
                """
                EMMAX implementation (in python)
                Single SNPs
                
                With interactions between SNP and possible cofactors.
                """
                q = 1  # Single SNP is being tested
                p = len(self.X.T)+q
                n = self.n
                n_p = n-p                                

                K = self.random_effects[1][1]
                eig_L = self._get_eigen_L_(K)
                res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
                delta = res['delta']
                print 'pseudo_heritability:',1.0/(1+delta)
                H_sqr = res['H_sqrt']
                H_sqrt_inv = H_sqr.I
                Y = H_sqrt_inv*self.Y        #The transformed outputs.
                h0_X = H_sqrt_inv*self.X
                (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X,Y)
                h0_betas = map(float,list(h0_betas))
                num_snps = len(snps)
                rss_list = repeat(h0_rss,num_snps)
                betas_list = [h0_betas]*num_snps
                var_perc = zeros(num_snps)
                H = array(H_sqrt_inv.T)
		missing_indices = []
                for i, snp in enumerate(snps):
                	X_1 = dot(snp,H)
                        (betas, rss, p, sigma) = linalg.lstsq(hstack([h0_X,matrix(X_1).T]),Y)
                        if not rss:
                                print 'No predictability in the marker, moving on...'
                                continue
                               	missing_indices.append(i)
                        rss_list[i] = rss[0]                
                        betas_list[i] = map(float,list(betas))
                        if (i+1)%(num_snps/10)==0:
                        	sys.stdout.write('.')
                        	sys.stdout.flush()
                sys.stdout.write('\n')
                rss_ratio = h0_rss/rss_list
                var_perc = 1-1/rss_ratio
                f_stats = (rss_ratio-1)*n_p/float(q)
                p_vals = stats.f.sf(f_stats,q,n_p)

                      
                return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list, 
                        'delta':delta, 'pseudo_heritability': 1.0/(1+delta), 'var_perc':var_perc,
                        'missing_indices':missing_indices}



	def emmax_f_test_2(self,snps):
		"""
		EMMAX implementation (in python)
		Single SNPs
		
		With interactions between SNP and possible cofactors.
		"""
		q = 1  # Single SNP is being tested
		p = len(self.X.T)+q
		n = self.n
		n_p = n-p				

		K = self.random_effects[1][1]
		eig_L = self._get_eigen_L_(K)
		res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
		delta = res['delta']
		print 'pseudo_heritability:',1.0/(1+delta)
		H_sqr = res['H_sqrt']
		H_sqrt_inv = H_sqr.I
		Y = H_sqrt_inv*self.Y	#The transformed outputs.
		h0_X = H_sqrt_inv*self.X
		(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X,Y)
		h0_betas = map(float,list(h0_betas))
		num_snps = len(snps)
		rss_list = repeat(h0_rss,num_snps)
		betas_list = [h0_betas]*num_snps
		var_perc = zeros(num_snps)
		chunk_size = len(Y)
		#ndot = lambda x,y: linalg.blas.get_blas_funcs(['gemm'])[0](1,x,y)
		#H = array(H_sqrt_inv.T)
		for i in range(0,len(snps),chunk_size):
			snps_chunk = matrix(snps[i:i+chunk_size])
			Xs = snps_chunk*(H_sqrt_inv.T)
        		#snps_chunk = array(snps[i:i+chunk_size])
			#Xs = ndot(snps_chunk,H)
        		for j in range(len(Xs)):
        			(betas, rss, p, sigma) = linalg.lstsq(hstack([h0_X,matrix(Xs[j]).T]),Y)
        			if not rss:
        				print 'No predictability in the marker, moving on...'
        				continue
        			rss_list[i+j] = rss[0]		
        			betas_list[i+j] = map(float,list(betas))
        			if (i+j+1)%(num_snps/10)==0:
        				sys.stdout.write('.')
        				sys.stdout.flush()
		sys.stdout.write('\n')
		rss_ratio = h0_rss/rss_list
		var_perc = 1-1/rss_ratio
		f_stats = (rss_ratio-1)*n_p/float(q)
		p_vals = stats.f.sf(f_stats,q,n_p)

		      
		return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list, 
			'delta':delta, 'pseudo_heritability': 1.0/(1+delta), 'var_perc':var_perc}


        
        def emmax_simple_f_test(self,snps):
                """
                EMMAX implementation (in python)
                Single SNPs
                
                With interactions between SNP and possible cofactors.
                """
                q = 1  # Single SNP is being tested
                p = len(self.X.T)+q
                n = self.n
                n_p = n-p                                

                K = self.random_effects[1][1]
                eig_L = self._get_eigen_L_(K)
                res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
                delta = res['delta']
                print 'pseudo_heritability:',1.0/(1+delta)
                H_sqr = res['H_sqrt']
                H_sqrt_inv = H_sqr.I
                Y = H_sqrt_inv*self.Y        #The transformed outputs.
                h0_X = H_sqrt_inv*self.X
                (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X,Y)
                num_snps = len(snps)
                rss_list = repeat(h0_rss,num_snps)
                for i, snp in enumerate(snps):
                        (betas, rss, p, sigma) = linalg.lstsq(hstack([h0_X,H_sqrt_inv*(matrix(snp).T)]),Y)
                        if not rss:
                                print 'No predictability in the marker, moving on...'
                                continue
                        rss_list[i] = rss[0]                
			if (i+1)%(num_snps/10)==0:
				sys.stdout.write('.')
				sys.stdout.flush()
		rss_ratio = h0_rss/rss_list
                f_stats = (rss_ratio-1)*n_p/float(q)
                p_vals = stats.f.sf(f_stats,q,n_p)

                      
                return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'delta':delta, 'pseudo_heritability': 1.0/(1+delta)}

        
        def emmax_f_test_parallel(self,snps,num_proc=2):
                """
                EMMAX implementation (in python)
                Single SNPs
                
                With interactions between SNP and possible cofactors.
                """
                assert len(self.random_effects)==2,"Expedited REMLE only works when we have exactly two random effects."
                q = 1  # Single SNP is being tested
                p = len(self.X.T)+q
                n = self.n
                n_p = n-p                                

                K = self.random_effects[1][1]
                eig_L = self._get_eigen_L_(K)
                res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
                delta = res['delta']
                print 'pseudo_heritability:',1.0/(1+delta)
                H_sqr = res['H_sqrt']
                H_sqrt_inv = H_sqr.I
                Y = H_sqrt_inv*self.Y        #The transformed outputs.
                h0_X = H_sqrt_inv*self.X
                (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X,Y)

                
                h0_betas = map(float,list(h0_betas))
                num_snps = len(snps)
                var_perc = zeros(num_snps)

                rss_ratio = mp.Array('d', range(num_snps))
                rss_list = mp.Array('d', [h0_rss]*num_snps)
                beta1_arr = mp.Array('d', range(num_snps))
                beta0_arr = mp.Array('d', range(num_snps))
                def _emmax_mp_loop_(self,snps,j):
                        for i, snp in enumerate(snps):
                                (betas, rss, p, sigma) = linalg.lstsq(hstack([h0_X,H_sqrt_inv*(matrix(snp).T)]),Y)
                                if not rss:
                                        print 'No predictability in the marker, moving on...'
                                        continue
                                rss_ratio[i+j] = h0_rss/rss             
                                rss_list[i+j] = rss[0]                
                                betas_list[i+j] = map(float,list(betas))

                snp_chunk_size = num_snps/num_proc 
                processes = []
                for i in range(num_proc-1):
                        snp_chunk = snps[i:i+snp_chunk_size]
                        betas_list_ = [h0_betas]*num_snps
                        rss_ratio_ = ones(num_snps)
                        rss_list_ = repeat(h0_rss,num_snps)
                        p = mp.Process(target=_emmax_mp_loop_, args=(snp_chunk,i))
                        p.start()
                        processes.append(p)
                        
                var_perc = 1-1/rss_ratio
                f_stats = (rss_ratio-1)*n_p/float(q)
                p_vals = stats.f.sf(f_stats,q,n_p)

                      
                return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list, 
                        'delta':delta, 'pseudo_heritability': 1.0/(1+delta), 'var_perc':var_perc}
                    
                        
         

        def emmax_two_snps(self,snps):
                """
                Every pair of SNPs, Vincent's results.
                """
                import util
                K = self.random_effects[1][1]
                eig_L = self._get_eigen_L_(K)
                num_snps = len(snps)
                
                res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
                delta = res['delta']
                pseudo_her = 1.0/(1+delta)
                H_sqr = res['H_sqrt']
                H_sqrt_inv = H_sqr.I
                Y = H_sqrt_inv*self.Y        #The transformed outputs.
                h0_X = H_sqrt_inv*self.X     #The transformed fixed inputs.
                (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X,Y)
                full_rss = h0_rss
                
                #Contructing result containers                
                f_stat_array = zeros((num_snps,num_snps))
                rss_array = zeros((num_snps,num_snps))
                p_val_array = zeros((num_snps,num_snps))
                var_perc_array = zeros((num_snps,num_snps))
                betas_list= [[0]*num_snps for s in snps] 
                
                q = 1 #Testing one SNP, given the other
                p = len(self.X.T)+q 
                n_p = self.n-p
                #Fill the diagonal with the single SNP emmax
                for i, snp in enumerate(snps):                        
                        snp_mat = H_sqrt_inv*matrix(snp).T #Transformed inputs
                        X = hstack([h0_X,snp_mat]) 
                        (betas, rss, rank, s) = linalg.lstsq(X,Y)
                        f_stat = ((h0_rss-rss)/q)/(rss/n_p)
                        p_val = stats.f.sf(f_stat,q,n_p)
                        p_val_array[i,i] = p_val[0]    
                        f_stat_array[i,i] = f_stat[0]               
                        rss_array[i,i] = rss[0]            
                        var_perc_array[i,i] = float(1-rss/h0_rss)
                        betas_list[i][i] = map(float,list(betas))
                
                #Fill the upper right part with poorest of two tests.
                p = len(self.X.T)+q+1 #Adding one SNP as a cofactor.
                n_p = self.n-p
                for i, snp1 in enumerate(snps):
                        anti_snp1 = util.anti_binary_snp(snp1)
                        snp_mat = H_sqrt_inv*(matrix(snp1).T) #Transformed inputs
                        h0_X1 = hstack([h0_X,snp_mat])                                                 
                        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X1,Y)
                        for j, snp2 in enumerate(snps):
                                if i==j: continue #Skip diagonals..
                                
                                if snp1!=snp2 and anti_snp1!=snp2: 
                                        snp_mat = H_sqrt_inv*(matrix(snp2).T) #Transformed inputs
                                        X = hstack([h0_X1,snp_mat]) 
                                        (betas, rss, rank, s) = linalg.lstsq(X,Y)
                                        f_stat = ((h0_rss-rss)/(q))/(rss/n_p)
                                        p_val = stats.f.sf(f_stat,q,n_p)
                                        p_val_array[i,j] = p_val[0]            
                                        f_stat_array[i,j] = f_stat[0]               
                                        rss_array[i,j] = rss[0]            
                                        var_perc_array[i,j] = float(1-rss/h0_rss)
                                        betas_list[i][j] = map(float,list(betas))
                                else:
                                        p_val_array[i,j] = 1            
                                        f_stat_array[i,j] = 0               
                                        rss_array[i,j] = h0_rss[0]        
                                        var_perc_array[i,j] = 0
                                        betas_list[i][j] = map(float,list(h0_betas))+[0.0]

                for i in range(len(snps)):
                        for j in range(len(snps[:i])):
                                if p_val_array[i,j]>p_val_array[j,i]:
                                        p_val_array[j,i] = p_val_array[i,j]
                                        f_stat_array[j,i] = f_stat_array[i,j]
                                        rss_array[j,i] = rss_array[i,j]
                                        var_perc_array[j,i] = var_perc_array[i,j]
                                        betas_list[j][i] = betas_list[i][j]
                                        
                
                no_snp_count = 0
                identical_snp_count = 0
                num_same_snp = 0 
                                        
                #Now the interaction.
                p = len(self.X.T)+q+2 #Adding two SNP as a cofactor.
                n_p = self.n-p
                for i, snp1 in enumerate(snps):
                        for j, snp2 in enumerate(snps[:i]):
                                anti_snp1 = util.anti_binary_snp(snp1)
                                anti_snp2 = util.anti_binary_snp(snp2)
                                pseudo_snp = [v1*v2 for v1,v2 in zip(snp1,snp2)]
                                if snp1 != snp2 and snp1 != anti_snp2 and len(set(pseudo_snp))>1 \
                                        and pseudo_snp!=snp1 and pseudo_snp!=snp2 and pseudo_snp!=anti_snp1 \
                                        and pseudo_snp!=anti_snp2:
                                        
                                        snp_mat = H_sqrt_inv*matrix(zip(snp1,snp2)) #Transformed inputs
                                        h0_X1 = hstack([h0_X,snp_mat])                                                 
                                        #(h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X1,Y)                                
                                        snp_mat = H_sqrt_inv*(matrix(pseudo_snp).T) #Transformed inputs
                                        X = hstack([h0_X1,snp_mat]) 
                                        (betas, rss, rank, s) = linalg.lstsq(X,Y)
                                        h0_rss = rss_array[j,i]                                        
                                        f_stat = ((h0_rss-rss)/(q))/(rss/n_p)
                                        p_val = stats.f.sf(f_stat,q,n_p)
                                        p_val_array[i,j] = p_val[0]            
                                        f_stat_array[i,j] = f_stat[0]               
                                        rss_array[i,j] = rss[0]            
                                        var_perc_array[i,j] = float(1-rss/h0_rss)
                                        betas_list[i][j] = map(float,list(betas))
                                else:
                                        if snp1 == snp2 or snp1 == anti_snp2:
                                                num_same_snp += 1
                                        elif len(set(pseudo_snp))==1:
                                                no_snp_count+=1
                                        elif pseudo_snp!=snp1 or pseudo_snp!=snp2 or \
                                                pseudo_snp!=anti_snp1 or pseudo_snp!=anti_snp2:
                                                identical_snp_count+=1
                                        p_val_array[i,j] = 1            
                                        f_stat_array[i,j] = 0               
                                        rss_array[i,j] = full_rss        
                                        var_perc_array[i,j] = 0
                                        betas_list[i][j] = map(float,list(h0_betas))+[0.0]
                
                
                tot_num = num_snps*(num_snps-1)/2.0
                print 'fail_fraction =%f, no_snp_fraction=%f, identical_pseudo_snp_fraction=%f, num_same_snp=%f'\
                        %((no_snp_count+identical_snp_count+num_same_snp)/tot_num,no_snp_count/tot_num,identical_snp_count/tot_num,num_same_snp/tot_num) 

                res_dict = {'vincent_ps':p_val_array, 'vincent_f_stats':f_stat_array, 'vincent_rss':rss_array, \
                            'vincent_betas': betas_list, 'vincent_var_perc':var_perc_array, \
                            'pseudo_heritability': pseudo_her}
                
                
                #Filling up the diagonal
                f_stat_array = zeros((num_snps,num_snps))
                rss_array = zeros((num_snps,num_snps))
                p_val_array = zeros((num_snps,num_snps))
                var_perc_array = zeros((num_snps,num_snps))
                betas_list= [[0]*num_snps for s in snps] 
                for i in range(num_snps):
                        p_val_array[i,i] = res_dict['vincent_ps'][i,i]
                        f_stat_array[i,i] = res_dict['vincent_f_stats'][i,i]
                        rss_array[i,i] = res_dict['vincent_rss'][i,i]
                        betas_list[i][i] = res_dict['vincent_betas'][i][i]
                        var_perc_array[i,i] = res_dict['vincent_var_perc'][i,i]
                 
                 
                q1 = 2 #Testing significance of two SNPs, compared to none
                p1 = len(self.X.T)+q1 
                n_p1 = self.n-p1
                q2 = 3 #Testing significance of two SNPs and the interaction, compared to none
                p2 = len(self.X.T)+q2 
                n_p2 = self.n-p2
                for i, snp1 in enumerate(snps):
                        for j, snp2 in enumerate(snps[:i]):
                                anti_snp1 = util.anti_binary_snp(snp1)
                                anti_snp2 = util.anti_binary_snp(snp2)
                                pseudo_snp = [v1*v2 for v1,v2 in zip(snp1,snp2)]
                                if snp1 != snp2 and snp1 != anti_snp2:                                        
                                        snp_mat = H_sqrt_inv*matrix(zip(snp1,snp2)) #Transformed inputs
                                        X = hstack([h0_X,snp_mat])                                                 
                                        (betas, rss, rank, s) = linalg.lstsq(X,Y)                                
                                        f_stat = ((full_rss-rss)/(q1))/(rss/n_p1)
                                        p_val = stats.f.sf(f_stat,q1,n_p1)
                                        p_val_array[j,i] = p_val[0]            
                                        f_stat_array[j,i] = f_stat[0]               
                                        rss_array[j,i] = rss[0]            
                                        var_perc_array[j,i] = float(1-rss/full_rss)
                                        betas_list[j][i] = map(float,list(betas))
                                else:
                                        rss = res_dict['vincent_rss'][j,i]
                                        f_stat = ((full_rss-rss)/(q1))/(rss/n_p1)
                                        p_val = stats.f.sf(f_stat,q1,n_p1)
                                        p_val_array[j,i] = p_val[0]            
                                        f_stat_array[j,i] = f_stat[0]               
                                        rss_array[j,i] = rss           
                                        var_perc_array[j,i] = float(1-rss/full_rss)
                                        betas_list[j][i] =res_dict['vincent_betas'][j,i]
                                
                                if  snp1 != snp2 and snp1 != anti_snp2 and len(set(pseudo_snp))>1 and pseudo_snp!=snp1 \
                                        and pseudo_snp!=snp2 and pseudo_snp!=anti_snp1 and pseudo_snp!=anti_snp2:
                                        snp_mat = H_sqrt_inv*matrix(zip(snp1,snp2,pseudo_snp)) #Transformed inputs
                                        X = hstack([h0_X,snp_mat])                                                 
                                        (betas, rss, rank, s) = linalg.lstsq(X,Y)                                
                                        f_stat = ((full_rss-rss)/(q2))/(rss/n_p2)
                                        p_val = stats.f.sf(f_stat,q2,n_p2)
                                        p_val_array[i,j] = p_val[0]            
                                        f_stat_array[i,j] = f_stat[0]               
                                        rss_array[i,j] = rss[0]            
                                        var_perc_array[i,j] = float(1-rss/full_rss)
                                        betas_list[i][j] = map(float,list(betas))
                                else:
                                        rss = rss_array[j,i]
                                        f_stat = ((full_rss-rss)/(q2))/(rss/n_p2)
                                        p_val = stats.f.sf(f_stat,q2,n_p2)
                                        p_val_array[i,j] = p_val[0]            
                                        f_stat_array[i,j] = f_stat[0]               
                                        rss_array[i,j] = rss            
                                        var_perc_array[i,j] = var_perc_array[j,i]
                                        betas_list[i][j] = betas_list[j][i]
                                         
                res_dict['ps']=p_val_array
                res_dict['f_stats']=f_stat_array
                res_dict['rss']=rss_array
                res_dict['betas']=betas_list
                res_dict['var_perc']=var_perc_array
                
                
                return res_dict



#        def emmax_two_snps_3(self,snps):
#                """
#                Every pair of SNPs against the null of no SNP.
#                """
#                try:
#                        import psyco
#                        psyco.full()
#                except ImportError:
#                        pass
#                assert len(self.random_effects)==2,"Expedited REMLE only works when we have exactly two random effects."
#
#                K = self.random_effects[1][1]
#                eig_L = self._get_eigen_L_(K)
#                q = 2 #Testing two SNPs
#                p = len(self.X.T)+q
#                n = self.n
#                n_p = n-p
#
#                res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
#                delta = res['delta']
#                H_sqr = res['H_sqrt']
#                H_sqrt_inv = H_sqr.I
#                Y = H_sqrt_inv*self.Y        #The transformed outputs.
#                h0_X = H_sqrt_inv*self.X
#                (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X,Y)
#                
#              
#                #First testing without interactions
#                f_stats = []
#                rss_list = []
#                betas_list = []
#                p_vals = []
#                var_perc = []
#                for i, snp1 in enumerate(snps):
#                        f_l = []
#                        rss_l = []
#                        betas_l = []
#                        p_vals_l = []
#                        var_perc_l = []
#                        for snp2 in snps[:i]:
#                                X = hstack([h0_X,H_sqrt_inv*matrix([[v1,v2] for v1,v2 in zip(snp1,snp2)])]) 
#                                res = _emmax_snp_test_(Y,X,h0_rss,n,p,q)
#                                p_vals_l.append(res['p_val'])            
#                                f_l.append(res['f_stat'])                
#                                rss_l.append(res['rss'])                
#                                betas_l.append(res['betas'])
#                                var_perc_l.append(res['var_perc'])
#                        f_stats.append(f_l)
#                        rss_list.append(rss_l)
#                        betas_list.append(betas_l)
#                        p_vals.append(p_vals_l)
#                        var_perc.append(var_perc_l)
#                        
#                return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list, 
#                        'delta':delta, 'pseudo_heritability': 1.0/(1+delta), 'var_perc':var_perc}
#                        
#         
#        def emmax_two_snps_2(self,snps):
#                """
#                Every pair of SNPs, compared to the best model contain one of the SNPs.
#                """
#                try:
#                        import psyco
#                        psyco.full()
#                except ImportError:
#                        pass
#                assert len(self.random_effects)==2,"Expedited REMLE only works when we have exactly two random effects."
#
#                K = self.random_effects[1][1]
#                eig_L = self._get_eigen_L_(K)
#                X = self.X
#                
#                q = 1 #Testing one SNP, given the other
#                p = len(X.T)+q+1 
#                n = self.n
#                n_p = n-p
#                
#                #Then also without interactions, but given the other SNP.
#                f_stats = []
#                rss_list = []
#                betas_list = []
#                p_vals = []
#                var_perc = []
#                pseudo_her_list = []
#                for i, snp1 in enumerate(snps):                        
#                        f_l = []
#                        rss_l = []
#                        betas_l = []
#                        p_vals_l = []
#                        var_perc_l = []
#
#                        self.X = hstack([X,matrix([[v1] for v1 in snp1])]) #Add the cofactor 
#                        res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
#                        delta = res['delta']
#                        pseudo_her_list.append(1.0/(1+delta))
#                        H_sqr = res['H_sqrt']
#                        H_sqrt_inv = H_sqr.I
#                        Y = H_sqrt_inv*self.Y        #The transformed outputs.
#                        h0_X = H_sqrt_inv*self.X
#                        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X,Y)
#                        for j, snp2 in enumerate(snps):                            
#                                if i==j: 
#                                        p_vals_l.append(0)            
#                                        f_l.append(0)                
#                                        rss_l.append(0)                
#                                        betas_l.append(0)
#                                        var_perc_l.append(0)
#                                else:
#                                        X1 = hstack([h0_X,H_sqrt_inv*matrix([[v2] for v2 in snp2])]) 
#                                        res = self._emmax_snp_test_(Y,X1,h0_rss,n,p,q)
#                                        p_vals_l.append(res['p_val'])            
#                                        f_l.append(res['f_stat'])                
#                                        rss_l.append(res['rss'])                
#                                        betas_l.append(res['betas'])
#                                        var_perc_l.append(res['var_perc'])
#                        f_stats.append(f_l)
#                        rss_list.append(rss_l)
#                        betas_list.append(betas_l)
#                        p_vals.append(p_vals_l)
#                        var_perc.append(var_perc_l)
#                
#                
#                best_f_stats = []
#                best_rss_list = []
#                best_betas_list = []
#                best_p_vals = []
#                best_var_perc = []
#                fixed_choice = []
#                for i in range(len(snps)):
#                        for j in range(i):
#                                if p_vals[i][j]>p_vals[j][i]:
#                                        bi=i
#                                        bj=j
#                                else:
#                                        bi=j
#                                        bj=i
#                                        
#                                best_f_stats.append(f_stats[bi][bj])
#                                best_rss_list.append(rss_list[bi][bj])
#                                best_betas_list.append(betas_list[bi][bj])
#                                best_p_vals.append(p_vals[bi][bj])
#                                best_var_perc.append(var_perc[bi][bj])
#                                fixed_choice.append(bi)
#                
#                
#                return {'ps':best_p_vals, 'f_stats':best_f_stats, 'rss':best_rss_list, \
#                                 'betas': best_betas_list, 'pseudo_heritability': pseudo_her_list, \
#                                 'var_perc':best_var_perc, 'fixed_choice':fixed_choice}
#                
#
#
#        def emmax_two_snp_interaction(self,snps): 
#                """
#                An interaction of a pair of SNPs compared to the model not containing the interaction term.
#                """
#                try:
#                        import psyco
#                        psyco.full()
#                except ImportError:
#                        pass
#                assert len(self.random_effects)==2,"Expedited REMLE only works when we have exactly two random effects."
#
#                K = self.random_effects[1][1]
#                eig_L = self._get_eigen_L_(K)
#                q = 1 #Testing one interaction
#                p = len(self.X.T)+q+2 #Given 2 SNPs
#                n = self.n
#                n_p = n-p
#
#                #Then also without interactions, but given the other SNP.
#                f_stats = []
#                rss_list = []
#                betas_list = []
#                p_vals = []
#                var_perc = []
#                pseudo_her_list = []
#                for i, snp1 in enumerate(snps):                        
#                        f_l = []
#                        rss_l = []
#                        betas_l = []
#                        p_vals_l = []
#                        var_perc_l = []
#                        pseudo_her_l
#                        for snp2 in snps[:i]:                            
#                                self.X = hstack([X,matrix([[v1,v2] for v1,v2 in zip(snp1,snp2)])]) #Add the cofactor 
#                                res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
#                                delta = res['delta']
#                                pseudo_her_l.append(1.0/(1+delta))
#                                H_sqr = res['H_sqrt']
#                                H_sqrt_inv = H_sqr.I
#                                Y = H_sqrt_inv*self.Y        #The transformed outputs.
#                                h0_X = H_sqrt_inv*self.X
#                                (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X,Y)
#                                X = hstack([h0_X,H_sqrt_inv*matrix([[v1*v2] for v1,v2 in zip(snp1,snp2)])]) 
#                                res = _emmax_snp_test_(Y,X,h0_rss,n,p,q)
#                                p_vals_l.append(res['p_val'])            
#                                f_l.append(res['f_stat'])                
#                                rss_l.append(res['rss'])                
#                                betas_l.append(res['betas'])
#                                var_perc_l.append(res['var_perc'])
#                        f_stats.append(f_l)
#                        rss_list.append(rss_l)
#                        betas_list.append(betas_l)
#                        p_vals.append(p_vals_l)
#                        var_perc.append(var_perc_l)
#                        pseudo_her_list.append(pseudo_her_l)
#                        
#                return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list,\
#                        'pseudo_heritability': pseudo_her_list, 'var_perc':var_perc}


        
        
def emmax(snps,phenotypes,K,cofactors=None,with_interactions=False,int_af_threshold=15):
        """
        Run EMMAX
        """
        lmm = LinearMixedModel(phenotypes)
        lmm.add_random_effect(K)
        if cofactors:
                for cofactor in cofactors:
                        lmm.add_factor(cofactor)
 
        print "Running EMMAX" 
        s1 = time.time()
        if with_interactions:
                res = lmm.emmax_f_test_w_interactions(snps,int_af_threshold=int_af_threshold)
        else:
                res = lmm.emmax_f_test_2(snps)        
        secs = time.time()-s1
        if secs>60:
                mins = int(secs)/60
                secs = secs - mins*60
                print 'Took %d mins and %f seconds.'%(mins,secs)
        else:
                print 'Took %f seconds.'%(secs)
        return res        
        
        
def linear_model(snps,phenotypes,cofactors=None):
        lm = LinearModel(phenotypes)
        if cofactors:
                for cofactor in cofactors:
                        lm.add_factor(cofactor)
        print "Running a standard linear model" 
        s1 = time.time()
        res = lm.fast_f_test(snps)
        secs = time.time()-s1
        if secs>60:
                mins = int(secs)/60
                secs = secs - mins*60
                print 'Took %d mins and %f seconds.'%(mins,secs)
        else:
                print 'Took %f seconds.'%(secs)
        return res        
        
def emmax_two_snps(snps,phenotypes,K,cofactors=None):
        lmm = LinearMixedModel(phenotypes)
        lmm.add_random_effect(K)
        if cofactors:
                for cofactor in cofactors:
                        lm.add_factor(cofactor)
 
        print "Running EMMAX on pairs of SNPs" 
        s1 = time.time()
        res = lmm.emmax_two_snps(snps)
        secs = time.time()-s1
        if secs>60:
                mins = int(secs)/60
                secs = secs - mins*60
                print 'Took %d mins and %f seconds.'%(mins,secs)
        else:
                print 'Took %f seconds.'%(secs)
        print 'pseudo_heritability:',res['pseudo_heritability']
        return res        
        

def emmax_snp_pair_plot(snps,positions,phenotypes,K,fm_scatter_plot_file=None,
                        fm_image_plot_file=None,scatter_plot_file=None,
                        image_plot_file=None,fig_format='pdf',show_plots=False,
                        display_threshold=0.05):
        """
        Plots a 2D plot.  
        
        Assumes phenotypes and genotypes are already coordinated.
        """
        r = emmax_two_snps(snps,phenotypes,K)
        thres = -log10(display_threshold)
        
        print 'Plotting scatter plot.'
        import pylab
        pylab.figure(figsize = (14,12))
        start_pos = min(positions)
        end_pos = max(positions)
        x_range = end_pos-start_pos
        pylab.axis([start_pos-0.025*x_range,end_pos+0.025*x_range,start_pos-0.025*x_range,end_pos+0.025*x_range])
        
        p_val_array = array(map(lambda y: map(lambda x: -log10(x),y),r['ps']))
        
        xs = []
        ys = []
        cs = []
        for i, pos, in enumerate(positions):
                xs.extend([pos]*len(positions))
                ys.extend(positions)
                cs.extend(p_val_array[i,:])

        for i in range(len(xs)-1,-1,-1):
                if cs[i]<thres:
                        del xs[i]
                        del ys[i]
                        del cs[i]
        
        l = zip(cs,xs,ys)
        l.sort()
        l = map(list,zip(*l))
        cs = l[0]
        xs = l[1]
        ys = l[2]
                 
        pylab.plot([start_pos,end_pos],[start_pos,end_pos],'k--')
        pylab.scatter(xs,ys, c=cs, alpha=0.8,linewidths=0)
        pylab.colorbar()
        if show_plots:
                 pylab.show()
        if scatter_plot_file:
                pylab.savefig(fm_scatter_plot_file,format=fig_format)
        pylab.clf()
        pylab.imshow(p_val_array,interpolation='nearest')
        pylab.colorbar()
        if show_plots:
                 pylab.show()
        if image_plot_file:
                pylab.savefig(fm_image_plot_file,format=fig_format)        
                
        pylab.clf()
        #Now Vincent's way..  not full model p-vals.
        p_val_array = array(map(lambda y: map(lambda x: -log10(x),y),r['vincent_ps']))
        
        xs = []
        ys = []
        cs = []
        for i, pos, in enumerate(positions):
                xs.extend([pos]*len(positions))
                ys.extend(positions)
                cs.extend(p_val_array[i,:])

        for i in range(len(xs)-1,-1,-1):
                if cs[i]<thres:
                        del xs[i]
                        del ys[i]
                        del cs[i]
        
        l = zip(cs,xs,ys)
        l.sort()
        l = map(list,zip(*l))
        cs = l[0]
        xs = l[1]
        ys = l[2]
                 
        pylab.plot([start_pos,end_pos],[start_pos,end_pos],'k--')
        pylab.scatter(xs,ys, c=cs, alpha=0.8,linewidths=0)
        pylab.axis([start_pos-0.025*x_range,end_pos+0.025*x_range,start_pos-0.025*x_range,end_pos+0.025*x_range])
        pylab.colorbar()
        if show_plots:
                 pylab.show()
        if scatter_plot_file:
                pylab.savefig(scatter_plot_file,format=fig_format)
        pylab.clf()
        pylab.imshow(p_val_array,interpolation='nearest')
        pylab.colorbar()
        if show_plots:
                 pylab.show()
        if image_plot_file:
                pylab.savefig(image_plot_file,format=fig_format)        



#def _test_emma_():
#        import time
#        #from rpy import r 
#        import rpy2.robjects as ro
#        import rpy2.robjects.numpy2ri
#        import util
#        r = ro.r
#        r.source("emma_fast.R")
#
#        phenotypes = [0,1,2,3,4,5,6,7,8,9,10]
#        snp = [0,0,0,0,0,1,1,1,1,1,1]
#        #K = diag([1]*len(phenotypes))
##        K = matrix([[1,0.5,0.1,0.2,0,0,0,0,0,0,0],
##                    [0.5,1,0.1,0.1,0,0,0,0,0,0,0],
##                    [0.1,0.1,1,0.1,0,0,0,0,0,0,0],
##                    [0.2,0.1,0.1,1,0,0,0,0,0,0,0],
##                    [0,0,0,0,1,0,0,0,0,0,0],
##                    [0,0,0,0,0,1,0,0,0,0,0],
##                    [0,0,0,0,0,0,1,0,0,0,0],
##                    [0,0,0,0,0,0,0,1,0,0,0.6],
##                    [0,0,0,0,0,0,0,0,1,0.3,0.6],
##                    [0,0,0,0,0,0,0,0,0.3,1,0.6],
##                    [0,0,0,0,0,0,0,0.6,0.6,0.6,1]])
##        print K.shape
#        phen_r = ro.conversion.py2ri(array([phenotypes]))
#        snps_r = ro.conversion.py2ri(array([snp]))
#        k_r = ro.conversion.py2ri(array(K))
#        phen_var_r = ro.conversion.py2ri(var(phenotypes,ddof=1)) 
#        reml_t = r['emma.REML.t']
#        s1 = time.time()
#        for i in range(1000):
#                res = reml_t(phen_r,snps_r,k_r,phen_var_r)
#        print time.time()-s1
#
#        lmm = LinearMixedModel(phenotypes)
#        lmm.add_random_effect(K)
#        lmm.add_factor(snp)
#        
#        #print lmm
#        #print eig_L
#        #print eig_R
#        #timeit 'lmm.get_expedited_REMLE(eig_L,eig_R)'
#        s1 = time.time()
#        for i in range(1000):
#                eig_L = lmm._get_eigen_L_(K)
#                eig_R = lmm._get_eigen_R_()
#                res = lmm.get_expedited_REMLE(eig_L,eig_R,return_pvalue=False)
#        print time.time()-s1
#        print res
#        print K.shape
#
#
#def _test_emma_2_():
#        #from rpy import r 
#        import rpy2.robjects as ro
#        import rpy2.robjects.numpy2ri
#        r = ro.r
#        r.source("emma_fast.R")
#        reml_t = r['emma.REML.t']
#        import gwa
#        import dataParsers as dp
#        import phenotypeData as pd
#        import util
#        sd = dp.parse_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t54.csv',format=0,filter=0.1)
#        filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_all_raw_070810.tsv"
#        phed = pd.readPhenotypeFile(filename)  
#        pid = 1
#        sd.coordinate_w_phenotype_data(phed,pid)
#        phenotypes = phed.getPhenVals(pid)
#        print phenotypes
#        snps = sd.getSnps()
#        chr_pos_list = sd.getChrPosList()
#        K = gwa.retrieve_kinship(phed.accessions,'/Users/bjarnivilhjalmsson/Projects/Data/250k/kinship_matrix_cm54.pickled')
#
#
##
##        phen_var_r = ro.conversion.py2ri(var(phenotypes,ddof=1)) 
##        phen_r = ro.conversion.py2ri(array([phenotypes]))
##        k_r = ro.conversion.py2ri(array(K))
##        snps_r = ro.conversion.py2ri(array(snps))
##        s1 = time.time()
##        res = reml_t(phen_r,snps_r,k_r,phen_var_r)
##        print util.r_list_2_dict(res)
##        print time.time()-s1
#
#        lmm = LinearMixedModel(phenotypes)
#        lmm.add_random_effect(K)
#        #lmm.add_factor(snps[0])
#        #lmm.add_factor(snps[-1])
#        #lmm.add_factor(snps[1])
##        s1 = time.time()
##        print lmm.expedited_REML_t_test(snps)
##        print time.time()-s1
#
#        s1 = time.time()
#        res = lmm.emmax_f_test(snps)
#        secs = time.time()-s1
#        if secs>60:
#                mins = int(secs)/60
#                secs = secs - mins*60
#                print 'Took %d mins and %f seconds.'%(mins,secs)
#        else:
#                print 'Took %f seconds.'%(secs)
#        
#        r = map(list,zip(*chr_pos_list))
#        print res['pseudo_heritability']
#        print res['betas'][0:50]
#        print res['var_perc'][0:50]
#        print res['ps'][0:50]
#        chromosomes = r[0]
#        positions = r[1]
#        import plotResults
#        plotResults.plot_raw_result(res['ps'],chromosomes,positions,pdf_file='/Users/bjarni.vilhjalmsson/tmp/EMMAX_FT10.pdf')
#        #print lmm
#        #print eig_L
#        #print eig_R
#        #timeit 'lmm.get_expedited_REMLE(eig_L,eig_R)'
##        for i in range(1000):
##                eig_L = lmm._get_eigen_L_(K)
##                eig_R = lmm._get_eigen_R_()
##                res = lmm.get_expedited_REMLE(eig_L,eig_R,return_pvalue=False)
##        print res
##        print K.shape


def _filter_k_(k,indices_to_keep):
        new_k = zeros((len(indices_to_keep),len(indices_to_keep)))
        for i in range(len(indices_to_keep)):
                for j in range(len(indices_to_keep)):
                        new_k[i,j]=k[indices_to_keep[i],indices_to_keep[j]]
        k = new_k
        return k


def filter_k_for_accessions(k,k_accessions,accessions):
        indices_to_keep = []
        for i, acc in enumerate(k_accessions):
                if acc in accessions:
                        indices_to_keep.append(i)                
        return _filter_k_(k, indices_to_keep)

def load_kinship_from_file(kinship_file,accessions=None):
        assert os.path.isfile(kinship_file), 'File not found.'
        sys.stdout.write("Loading K.\n")
        sys.stdout.flush()
        f = open(kinship_file,'r')
        l = cPickle.load(f)
        f.close()
        k = l[0]
        k_accessions = l[1]
        return filter_k_for_accessions(k,k_accessions,accessions)
                        



def calc_kinship_old(snps):
        """
        Requires EMMA to be installed.
        
        Uses rpy2 
        """
        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri
        import env
        r_source = robjects.r('source')
        a = array(snps)
        r_source(env.env['script_dir']+"emma_fast.R")
        r_emma_kinship = robjects.r['emma.kinship']
        return array(r_emma_kinship(a))
        #from rpy import r 
        #r.source(script_dir+"emma_fast.R")
        #return r.emma_kinship(a)


def calc_kinship(snps):
        """
        Requires EMMA to be installed.
        
        Uses the old rpy
        """
	import env
	a = array(snps)
	from rpy import r 
        r.source(env.env['script_dir']+"emma_fast.R")
        return r.emma_kinship(a)


def _test_two_snps_emma_():
        import gwa
        import dataParsers as dp
        import phenotypeData as pd
        import util
        filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_all_raw_070810.tsv"
        phed = pd.readPhenotypeFile(filename)  
        pid = 5
        filter_accessions = phed.getNonNAEcotypes(pid)
        sd = dp.parse_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t54.csv',format=0,filter_accessions=filter_accessions)
        sd.coordinate_w_phenotype_data(phed,pid)
        phenotypes = phed.getPhenVals(pid)
        snps = sd.getSnps()
        chr_pos_list = sd.getChrPosList()
        K = load_kinship_from_file('/Users/bjarnivilhjalmsson/Projects/Data/250k/kinship_matrix_cm54.pickled',phed.accessions)        
        
        chromosome = 5
        #r = sd.get_region_pos_snp_dict(4,250000,290000)
        r = sd.get_region_pos_snp_dict(chromosome,3160000,3200000)
        snps = r['snps']
        positions = r['positions']
        print 'Found %d SNPs'%len(snps)
 
        
        res = emmax(snps,phenotypes,K)
        emmax_scores = map(lambda x: -math.log10(x),res['ps'])
        res = linear_model(snps,phenotypes)
        lm_scores = map(lambda x: -math.log10(x),res['ps'])

        import env
        tmp_dir = env.env['tmp_dir']
        import regionPlotter as rp
        rp.plot_simple_region([lm_scores,emmax_scores],[positions,positions[:]],chromosome,caption='FLC_region, FT 10C',png_file=tmp_dir+'FLC_loci_FT10.png',tair_file=tmp_dir+'FLC_loci_FT10_tair.txt')        

        fm_scat_plot = tmp_dir+'FLC_loci_FT10_fm_scatter.pdf'
        scat_plot = tmp_dir+'FLC_loci_FT10_vincent_scatter.pdf'
        fm_im_plot = tmp_dir+'FLC_loci_FT10_fm_image.pdf'
        im_plot = tmp_dir+'FLC_loci_FT10_vincent_image.pdf'
        emmax_snp_pair_plot(snps,positions,phenotypes,K,
                            scatter_plot_file=scat_plot, fm_scatter_plot_file=fm_scat_plot, 
                            fm_image_plot_file=fm_im_plot, image_plot_file=im_plot)
  
  
  
  
def _test_cofactor_emma_():
        import gwa
        import dataParsers as dp
        import phenotypeData as pd
        import util
        filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_all_raw_070810.tsv"
        phed = pd.readPhenotypeFile(filename)  
        pid = 5
        filter_accessions = phed.getNonNAEcotypes(pid)
        sd = dp.parse_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t54.csv',format=0,filter_accessions=filter_accessions,filter=1.0)
        sd.coordinate_w_phenotype_data(phed,pid)
        phenotypes = phed.getPhenVals(pid)
        snps = sd.getSnps()
        K = load_kinship_from_file('/Users/bjarnivilhjalmsson/Projects/Data/250k/kinship_matrix_cm54.pickled',phed.accessions)        
        chr_pos_list = sd.getChrPosList()
        r = map(list,zip(*chr_pos_list))
        chromosomes = r[0]
        positions = r[1]
        flc_snp = sd.get_snp_at(5,3188328)
        import random
        res = emmax(snps,phenotypes,K,cofactors=[flc_snp],with_interactions=True)
        import gwaResults
        res = gwaResults.Result(scores=res['ps'], snps_data=sd)
        res.filter_attr("mafs",15)
        res.neg_log_trans()
        res.plot_manhattan(png_file='/Users/bjarni.vilhjalmsson/tmp/EMMAX_FT10_FLC_loci_w_interactions.png',
                           percentile=90, type="pvals", ylab="$-$log$_{10}(p)$", plot_bonferroni=True)
#        import plotResults
#        plotResults.plot_raw_result(res['ps'],chromosomes,positions,pdf_file='/Users/bjarni.vilhjalmsson/tmp/EMMAX_LD_FLC_loci.pdf')
#        res = emmax(snps,phenotypes,K)
#        plotResults.plot_raw_result(res['ps'],chromosomes,positions,pdf_file='/Users/bjarni.vilhjalmsson/tmp/EMMAX_LD.pdf')
#        res = linear_model(snps,phenotypes)
#        lm_scores = map(lambda x: -math.log10(x),res['ps'])
#
#        import env
#        tmp_dir = env.env['tmp_dir']
#        import regionPlotter as rp
#        rp.plot_simple_region([lm_scores,emmax_scores],[positions,positions[:]],chromosome,caption='FLC_region, FT 10C',png_file=tmp_dir+'FLC_loci_FT10.png',tair_file=tmp_dir+'FLC_loci_FT10_tair.txt')        
#
#        fm_scat_plot = tmp_dir+'FLC_loci_FT10_fm_scatter.pdf'
#        scat_plot = tmp_dir+'FLC_loci_FT10_vincent_scatter.pdf'
#        fm_im_plot = tmp_dir+'FLC_loci_FT10_fm_image.pdf'
#        im_plot = tmp_dir+'FLC_loci_FT10_vincent_image.pdf'
#        emmax_snp_pair_plot(snps,positions,phenotypes,K,
#                            scatter_plot_file=scat_plot, fm_scatter_plot_file=fm_scat_plot, 
#                            fm_image_plot_file=fm_im_plot, image_plot_file=im_plot)
  


def _test_kinship_():
        import dataParsers as dp
        import random
        sd = dp.parse_binary_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t54.csv.binary',filter=0.1)
        snps = sd.getSnps()
        #snps = random.sample(sd.getSnps(),10000)
        print 'Calculating kinship'
        K = calc_kinship(snps)
        print K
  
  
if __name__ == "__main__":
        _test_kinship_()
