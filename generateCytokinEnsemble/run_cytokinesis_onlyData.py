# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import os
import pickle

class generateCyto(object):
    
    def __init__(self,fit_until_this_t=14,nrealiz=5):

        self.fit_until_this_t = fit_until_this_t
        self.nrealiz=nrealiz

        ##################################### fixed parameters here, variable parameters are passed in run()
        self.seed = 123456
        self.var_start_sim = 0
        self.Tabsorb = 0
        self.x_upper_bound = 1700
        self.TminLength = 0
        self.deltaTswitch = 5
        self.varDeltaTswitch =0
        self.maxDeltaTswitch = self.deltaTswitch
        self.Tswitch = 400
        self.x0_for_reverse = 0       
        self.max_rtime = fit_until_this_t        
        self.traFreq =10
        # wfreq only for cov
        self.wfreq =200
        self.total_time = 4000
        
        self.dt = 0.05
        
        self.D1 = 0.04
        self.D2 = 0.04      
        self.D = self.D2

        self.Nb = 40
        self.A = 0.12e-3 
        self.r0 = 14.5
   
        self.K = 10
        self.xi = 3.75e-4        
        
        self.gamma1pre = 2*self.K*self.xi*10
        self.gamma1 = 0
        self.gamma2pre = 20*2*np.pi*self.A*self.xi*self.Nb 
        self.gamma2 = 20*2*np.pi*self.A*self.xi*self.Nb      
        
    def run(self,params):
        kp1,kd1,kp,kd,alpha = params
        
        self.kp1 = kp1
        self.kp = kp
        self.kd1 = kd1
        self.kd = kd
        self.alpha = alpha  
        
        ##################
        self.c_start_sim = kp1/kd1 
        print("c_start_sim", self.c_start_sim      )  
        
        self.x_start_sim = self.r0 - self.gamma2pre*self.c_start_sim**alpha/self.gamma1pre
        print("R0 * c0",self.x_start_sim*self.c_start_sim)
        self.rc_fac = self.x_start_sim*self.c_start_sim        
        
        
        ###########################################      
        os.system('./src_4_cytokinesis_PNAS2/prog %f %f %f %f %d %f %f %f %d %d %f %f %f %f %d %d %f %d %d %d %d %f %d %f %f %f %f %f %f %f %f %f %f' % (
        self.xi, self.K,self.x_start_sim,self.A, self.Nb, self.kp,self.kd,self.c_start_sim,
        self.nrealiz,self.seed,self.dt,self.D1,self.D2,
        self.x0_for_reverse,self.total_time,self.Tswitch,self.deltaTswitch,
        self.TminLength,self.wfreq,self.traFreq,self.max_rtime,self.x_upper_bound,self.Tabsorb,
        self.var_start_sim,self.varDeltaTswitch,self.maxDeltaTswitch,
        self.gamma1,self.gamma2,self.alpha,self.kp1,self.kd1,self.gamma1pre,self.gamma2pre))  
        ########################################### 
        
        ##### load data into class for later saving
        numfile_x = 'num_trj.txt'
        x_num = np.loadtxt(numfile_x, unpack=1)
        
        self.x_data = x_num
        


#########################  M  A I  N  ########################################
##############################################################################
##############################################################################

### instantiate class
nrealiz = 10000
names = ['drift','neg_two']#['drift','bessel','neg_two']
#kp1,kd1,kp,kd,alpha
p0 = [400,10,2000,20,1]
p1 = [0.5,0.01,0,0,1]
p2 = [0.7,0.1,0,0,2]
all_params = [p0,p2] #[p0,p1,p2]

gC = generateCyto(fit_until_this_t=300,nrealiz=nrealiz) 

### produce data
for params,name in zip(all_params,names):
    gC.run(params)
    save_name = "x_"+name+"highres10000.pickle"
    pickle.dump(gC.x_data,open(save_name,"wb"))    


