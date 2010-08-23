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
                p_value = 1-stats.f.cdf(f_stat,1,self.n-self.p+1)
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

        
        def _get_eigen_R_(self,hat_matrix=None,complete=True):
                if not hat_matrix:
                        hat_matrix =  self.get_hat_matrix()
                S = mat(identity(self.n))-hat_matrix        #S=I-X(X'X)^{-1}X'
                evals,evecs = linalg.eigh(S*(self.random_effects[1][1]+self.random_effects[0][1])*S) #eigen of S(K+I)S
                #print mat(evecs),'\n'
                #print (mat(evecs).T[self.p:]),'\n'
                #print evals[self.p:],'\n'
                return {'values':map(lambda x: x-1, evals[self.p:]),'vectors':(mat(evecs).T[self.p:])}   #Because of S(K+I)S?


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

        def get_expedited_REMLE(self, x, eig_L, eig_R, ngrids=100, llim=-4, ulim=10, esp=1e-6, return_pvalue=True):
                """
                Get REML estimates for the effect sizes, as well as the random effect contributions.
                
                Using the EMMA algorithm.                
                """
                X = hstack([self.X,matrix([[v] for v in x])])
                K = self.random_effects[1][1]
                t = K.shape[0] #number of rows
                q = X.shape[1] #number of columns
                n = self.n
                #assert K.shape[0]==K.shape[1]==X.shape[0]==n,'Dimensions are wrong.'

                etas = eig_R['vectors']*self.Y
                sq_etas = [float(eta)*float(eta) for eta in etas]
                log_deltas =  [float(i)/ngrids*(ulim-llim)+llim  for i in range(0,ngrids+1)] #a list of deltas to search
                #assert len(log_deltas)==ngrids+1,'Delta list size error.'
                deltas = map(math.exp,log_deltas)
                eig_vals = list(eig_R['values'])
                #assert len(eig_vals)==n-q,'Number of eigenvalues is incorrect.'
                
                #LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
                
                c_1 = 0.5*(n-q)*(math.log((n-q)/(2.0*math.pi))-1)
                
               
                def calc_ll(d):
                        """
                        Calculates the likelihood, and the derivative given a delta.
                        """
                        sum_1 = 0
                        sum_2 = 0
                        sum_3 = 0
                        sum_4 = 0
                        for j in range(n-q):
                                v_1 = eig_vals[j]+d
                                v_2 = sq_etas[j]/v_1
                                sum_1 += v_2
                                sum_2 += math.log(v_1)
                                sum_3 += v_2/v_1
                                sum_4 += 1.0/v_1 
                        ll = c_1-(0.5)*((n-q)*math.log(sum_1)+sum_2)  #log-likelihoods (eq. 7 from paper)
                        dll = 0.5*((n-q)*sum_3/sum_1-sum_4)  #diffrentiated log-likelihoods (eq. 9 from paper)
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
                        new_opt_delta = optimize.newton(self._dll_, opt_delta, args=(eig_vals,sq_etas), tol=esp,maxiter=50)
                        if deltas[opt_i]<new_opt_delta<deltas[opt_i+1]:
                                opt_delta = new_opt_delta
                                opt_ll = self._ll_(opt_delta,eig_vals,sq_etas)
                        else:
                                raise Exception('Local maximum outside of suggested area??')
                        #Implement Newton-Raphson algorithm here?
                        if opt_ll<max_ll:
                                opt_delta = deltas[max_ll_i]
                else:
                        opt_delta = deltas[max_ll_i]
                        opt_ll = max_ll

                        
                #eig_vals.reverse()
                l = map(lambda x,y: x/(y+opt_delta),sq_etas,eig_vals)
                opt_vg = sum(l)/(n-q)  #vg   
                opt_ve = opt_vg*opt_delta  #ve
                 
                H_inverse = eig_L['vectors'].T*diag([1/(ev+opt_delta) for ev in eig_L['values']])*eig_L['vectors']
                XX = X.T*(H_inverse*X)
                iXX = XX.I
                self.beta_est = iXX*X.T*(H_inverse*self.Y)
                x_beta = X*self.beta_est
                x_beta_var = var(x_beta,ddof=1)
                var_perc = x_beta_var/self.y_var
                t_stat = self.beta_est[q-1]/math.sqrt(float(iXX[q-1,q-1])*opt_vg)
                if return_pvalue:
                        p_val = stats.t.sf(abs(t_stat), n-q)*2
                        res_dict = {'max_ll':opt_ll, 'delta':opt_delta, 've':opt_ve, 'vg':opt_vg,
                                    'var_perc':var_perc, 't_stat':float(t_stat), 'p_val':float(p_val)}
                else:
                        res_dict = {'max_ll':opt_ll, 'delta':opt_delta, 've':opt_ve, 'vg':opt_vg,
                            'var_perc':var_perc, 't_stat':float(t_stat)}
                
                return res_dict

        
        def expedited_REML_t_test(self,snps):
                """
                Single SNP analysis
                """
                assert len(self.random_effects)==2,"Expedited REMLE only works when we have exactly two random effects."
                K = self.random_effects[1][1]
                eig_L = self._get_eigen_L_(K)
                eig_R = self._get_eigen_R_()
                abs_t_stats = []
                for snp in snps:
                        res = self.get_expedited_REMLE(snp,eig_L,eig_R,return_pvalue=False)
                        abs_t_stats.append(abs(res['t_stat']))                
                

                p_vals = map(lambda x: x*2, stats.t.sf(abs_t_stats, len(eig_L)))
                return p_vals,res
                
                

class GeneralizedLinearModel(object):
       
       def __init__(self):
               pass         

def simple_lm(snps,phenValues):
        """
        Linear model implementation.
        """
        pass


def emmax_r(snps,phenotypeValues):
        pass

def _test_emma_():
        import time
        #from rpy import r 
        import rpy2.robjects as ro
        import rpy2.robjects.numpy2ri
        r = ro.r
        r.source("emma_fast.R")

        phenotypes = [0,1,2,3,4,5,6,7,8,9,10]
        snp = [0,0,0,0,0,1,1,1,1,1,1]
        #K = diag([1]*len(phenotypes))
        K = matrix([[1,0.5,0.1,0.2,0,0,0,0,0,0,0],
                    [0.5,1,0.1,0.1,0,0,0,0,0,0,0],
                    [0.1,0.1,1,0.1,0,0,0,0,0,0,0],
                    [0.2,0.1,0.1,1,0,0,0,0,0,0,0],
                    [0,0,0,0,1,0,0,0,0,0,0],
                    [0,0,0,0,0,1,0,0,0,0,0],
                    [0,0,0,0,0,0,1,0,0,0,0],
                    [0,0,0,0,0,0,0,1,0,0,0.6],
                    [0,0,0,0,0,0,0,0,1,0.3,0.6],
                    [0,0,0,0,0,0,0,0,0.3,1,0.6],
                    [0,0,0,0,0,0,0,0.6,0.6,0.6,1]])
        print K.shape
        phen_r = ro.conversion.py2ri(array([phenotypes]))
        snps_r = ro.conversion.py2ri(array([snp]))
        k_r = ro.conversion.py2ri(array(K))
        phen_var_r = ro.conversion.py2ri(var(phenotypes,ddof=1)) 
        reml_t = r['emma.REML.t']
        s1 = time.time()
        for i in range(1000):
                res = reml_t(phen_r,snps_r,k_r,phen_var_r)
        print time.time()-s1
        print res        

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
        sd = dp.parse_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t54.csv',format=0,filter=0.0005)
        filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_all_070110.tsv"
        phed = pd.readPhenotypeFile(filename)  
        sd.coordinate_w_phenotype_data(phed,5)
        phenotypes = phed.getPhenVals(5)
        phen_var_r = ro.conversion.py2ri(var(phenotypes,ddof=1)) 
        snps = sd.getSnps()
        K = gwa.retrieve_kinship(phed.accessions,'/Users/bjarnivilhjalmsson/Projects/Data/250k/kinship_matrix_cm54.pickled')

#        phenotypes = [0,1,2,3,4,5,6,7,8,9,10]
#        snp = [0,0,0,0,0,1,1,1,1,1,1]
#        K = diag([1]*len(phenotypes))
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
        print K.shape
        s1 = time.time()
        for snp in snps:
                phen_r = ro.conversion.py2ri(array([phenotypes]))
                snps_r = ro.conversion.py2ri(array([snp]))
                k_r = ro.conversion.py2ri(array(K))
                res = reml_t(phen_r,snps_r,k_r,phen_var_r)
        print res
        print time.time()-s1

        lmm = LinearMixedModel(phenotypes)
        lmm.add_random_effect(K)
        s1 = time.time()
        print lmm.expedited_REML_t_test(snps)
        print time.time()-s1
        
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
