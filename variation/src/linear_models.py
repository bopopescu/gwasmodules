"""
Contains functions to perform various linear regression schemes, such as simple, and mixed models.
"""
from scipy import *
from scipy import linalg
from scipy import stats
from scipy import optimize
import math
        

class LinearModel(object):
        """
        A simple linear model
        """
        def __init__(self,Y=None):
                """
                The fixed effects should be a list of fixed effect lists (SNPs)
                """
                self.n = len(Y)
                self.Y = matrix([[y] for y in Y])
                self.X = matrix([[1] for y in Y]) #The intercept
                self.p = 1
                self.beta_est = None


        def add_factor(self,x):
                """
                Adds an explanatory variable to the X matrix.
                """
                self.X = hstack([self.X,matrix([[v] for v in x])])
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
        
        
        def get_generalized_transformation(self,variance_matrix):
                pass
                
                
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
                
                #A list of random effect type, and the cov matrix.
                self.random_effects = [('normal',matrix(identity(self.n)))] #The first random effect is the IID error.
                        

        def add_random_effect(self,cov_matrix=None,effect_type='normal'):
                if effect_type!='normal':
                        raise Exception('Currently, only Normal random effects are allowed.')
                self.random_effects.append((effect_type,cov_matrix))

        def least_square_estimate(self): 
                raise Exception("LSE not applicable for mixed models")


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
                if xs:
                        X = hstack([self.X,matrix([[v] for x in xs for v in x ])])
                else:
                        X = self.X
                eig_R = self._get_eigen_R_(X)
                K = self.random_effects[1][1]
                t = K.shape[0] #number of rows
                q = X.shape[1] #number of columns
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
                try:
                        import psyco
                        psyco.full()
                except ImportError:
                        pass
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
        


        def emmax_f_test(self,snps):
                """
                EMMAX implementation (in python)
                Single SNPs
                """
                q = 1  # Single SNPs
                p = self.n-q
                try:
                        import psyco
                        psyco.full()
                except ImportError:
                        pass
                assert len(self.random_effects)==2,"Expedited REMLE only works when we have exactly two random effects."
                K = self.random_effects[1][1]
                eig_L = self._get_eigen_L_(K)
                res = self.get_expedited_REMLE(eig_L=eig_L) #Get the variance estimates..
                delta = res['delta']
                H_sqr = res['H_sqrt']
                H_sqrt_inv = H_sqr.I
                Y = H_sqrt_inv*self.Y        #The transformed outputs.
                h0_X = H_sqrt_inv*self.X
                (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X,Y)
                #Y = Y - h0_X*h0_betas
                #h0_betas = float(h0_betas)
                f_stats = []
                rss_list = []
                betas_list = []
                p_vals = []
                var_perc = []
                for snp in snps:
                        snp_mat = H_sqrt_inv*matrix([[v] for v in snp]) #Transformed inputs
                        X = hstack([h0_X,snp_mat]) 
                        #X = H_sqrt_inv*matrix([[v] for v in snp]) #Transformed inputs 
                        (betas, rss, rank, s) = linalg.lstsq(X,Y)
                        f_stat = ((h0_rss-rss)/(q))/(rss/p)
                        p_val = stats.f.sf(f_stat,q,p)
                        p_vals.append(p_val[0])            
                        f_stats.append(f_stat[0])                
                        rss_list.append(rss[0])                
                        #betas_list.append([h0_betas]+map(float,list(betas)))
                        betas_list.append(map(float,list(betas)))
                        var_perc.append(float(1-rss/h0_rss))

                        
                return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list, 
                        'delta':delta, 'pseudo_heritability': 1.0/(1+delta), 'var_perc':var_perc}
                        
         

                
                

#class GeneralizedLinearModel(LinearMixedModel):
               
        

def _test_emma_():
        import time
        #from rpy import r 
        import rpy2.robjects as ro
        import rpy2.robjects.numpy2ri
        import util
        r = ro.r
        r.source("emma_fast.R")

        phenotypes = [0,1,2,3,4,5,6,7,8,9,10]
        snp = [0,0,0,0,0,1,1,1,1,1,1]
        #K = diag([1]*len(phenotypes))
#        K = matrix([[1,0.5,0.1,0.2,0,0,0,0,0,0,0],
#                    [0.5,1,0.1,0.1,0,0,0,0,0,0,0],
#                    [0.1,0.1,1,0.1,0,0,0,0,0,0,0],
#                    [0.2,0.1,0.1,1,0,0,0,0,0,0,0],
#                    [0,0,0,0,1,0,0,0,0,0,0],
#                    [0,0,0,0,0,1,0,0,0,0,0],
#                    [0,0,0,0,0,0,1,0,0,0,0],
#                    [0,0,0,0,0,0,0,1,0,0,0.6],
#                    [0,0,0,0,0,0,0,0,1,0.3,0.6],
#                    [0,0,0,0,0,0,0,0,0.3,1,0.6],
#                    [0,0,0,0,0,0,0,0.6,0.6,0.6,1]])
#        print K.shape
        phen_r = ro.conversion.py2ri(array([phenotypes]))
        snps_r = ro.conversion.py2ri(array([snp]))
        k_r = ro.conversion.py2ri(array(K))
        phen_var_r = ro.conversion.py2ri(var(phenotypes,ddof=1)) 
        reml_t = r['emma.REML.t']
        s1 = time.time()
        for i in range(1000):
                res = reml_t(phen_r,snps_r,k_r,phen_var_r)
        print time.time()-s1

        lmm = LinearMixedModel(phenotypes)
        lmm.add_random_effect(K)
        lmm.add_factor(snp)
        
        #print lmm
        #print eig_L
        #print eig_R
        #timeit 'lmm.get_expedited_REMLE(eig_L,eig_R)'
        s1 = time.time()
        for i in range(1000):
                eig_L = lmm._get_eigen_L_(K)
                eig_R = lmm._get_eigen_R_()
                res = lmm.get_expedited_REMLE(eig_L,eig_R,return_pvalue=False)
        print time.time()-s1
        print res
        print K.shape


def _test_emma_2_():
        import time
        #from rpy import r 
        import rpy2.robjects as ro
        import rpy2.robjects.numpy2ri
        r = ro.r
        r.source("emma_fast.R")
        reml_t = r['emma.REML.t']
        import gwa
        import dataParsers as dp
        import phenotypeData as pd
        import util
        sd = dp.parse_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t54.csv',format=0,filter=0.1)
        filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_all_raw_070810.tsv"
        phed = pd.readPhenotypeFile(filename)  
        pid = 1
        sd.coordinate_w_phenotype_data(phed,pid)
        phenotypes = phed.getPhenVals(pid)
        print phenotypes
        snps = sd.getSnps()
        chr_pos_list = sd.getChrPosList()
        K = gwa.retrieve_kinship(phed.accessions,'/Users/bjarnivilhjalmsson/Projects/Data/250k/kinship_matrix_cm54.pickled')

#
#        phen_var_r = ro.conversion.py2ri(var(phenotypes,ddof=1)) 
#        phen_r = ro.conversion.py2ri(array([phenotypes]))
#        k_r = ro.conversion.py2ri(array(K))
#        snps_r = ro.conversion.py2ri(array(snps))
#        s1 = time.time()
#        res = reml_t(phen_r,snps_r,k_r,phen_var_r)
#        print util.r_list_2_dict(res)
#        print time.time()-s1

        lmm = LinearMixedModel(phenotypes)
        lmm.add_random_effect(K)
        #lmm.add_factor(snps[0])
        #lmm.add_factor(snps[-1])
        #lmm.add_factor(snps[1])
#        s1 = time.time()
#        print lmm.expedited_REML_t_test(snps)
#        print time.time()-s1

        s1 = time.time()
        res = lmm.emmax_f_test(snps)
        secs = time.time()-s1
        if secs>60:
                mins = int(secs)/60
                secs = secs - mins*60
                print 'Took %d mins and %f seconds.'%(mins,secs)
        else:
                print 'Took %f seconds.'%(secs)
        
        r = map(list,zip(*chr_pos_list))
        print res['pseudo_heritability']
        print res['betas'][0:50]
        print res['var_perc'][0:50]
        print res['ps'][0:50]
        chromosomes = r[0]
        positions = r[1]
        import plotResults
        plotResults.plot_raw_result(res['ps'],chromosomes,positions,pdf_file='/Users/bjarni.vilhjalmsson/tmp/EMMAX_FT10.pdf')
        #print lmm
        #print eig_L
        #print eig_R
        #timeit 'lmm.get_expedited_REMLE(eig_L,eig_R)'
#        for i in range(1000):
#                eig_L = lmm._get_eigen_L_(K)
#                eig_R = lmm._get_eigen_R_()
#                res = lmm.get_expedited_REMLE(eig_L,eig_R,return_pvalue=False)
#        print res
#        print K.shape


if __name__ == "__main__":
        _test_emma_2_()
