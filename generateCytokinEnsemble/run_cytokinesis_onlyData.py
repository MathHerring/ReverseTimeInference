# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import pylab as plt
import os
import scipy.special as sp
#import scipy.stats as st

import sys
import pickle


class PureDriftAndDiffusion(object):
    
    def __init__(self,fit_until_this_t=14,nrealiz=5):


        self.fit_until_this_t = fit_until_this_t
        self.nrealiz=nrealiz



        #####################################
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
        
        ### paper prameters ####
        """
        kd = 0.1 1/sec
        kp = 2.5 1/(s mu m)
        Nb = 40
        A = 1.2*10^3 nN mum^2
        xi = 3.75 10^-4 mum
        c0 = 25 1/mum
        K = 10 nN mum^-1
        R0 = 14.5 mum
        """        
        
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
        
        datafile_hittingT='hittingTime.txt'
        self.hitting_times = np.genfromtxt(datafile_hittingT)
        self.x_start_mean_values = np.ones(self.nrealiz)*self.x_start_sim
        
        numfile_x = 'num_trj.txt'
        x_num = np.loadtxt(numfile_x, unpack=1)
        
        self.x_data = x_num
        
        numfile_mean = 'num_mean.txt'
        self.mean = np.loadtxt(numfile_mean, unpack=1)

        numfile_var = 'num_var.txt'
        self.var = np.loadtxt(numfile_var, unpack=1)
        
        numfile_skew = 'num_skew.txt'
        self.skew = np.loadtxt(numfile_skew, unpack=1)         

        numfile_meanConc = 'num_meanConc.txt'
        self.meanConc = np.loadtxt(numfile_meanConc, unpack=1)  
        
        numfile_conc = 'num_conc.txt'
        self.conc = np.loadtxt(numfile_conc, unpack=1)
        






#########################  M  A I  N  ########################################
##############################################################################
##############################################################################

bool_init = True
bool_run = True


if(bool_init):


    nrealiz = 10000
    
    names = ['drift','neg_two']#['drift','bessel','neg_two']
    #kp1,kd1,kp,kd,alpha
    p0 = [400,10,2000,20,1]
    p1 = [0.5,0.01,0,0,1]
    p2 = [0.7,0.1,0,0,2]
    all_params = [p0,p2] #[p0,p1,p2]

    ##########################################################################    
    # Instantiate Classes###################################################  
    ##########################################################################  
    PD = PureDriftAndDiffusion(fit_until_this_t=300,nrealiz=nrealiz) #fit_until_this_t=1000


    """
    kp1= params[0]
    kd1= params[1]
    kp = params[2]
    kd = params[3]
    alph = params[4] 

    if(name=='bessel'):
        gam=PD.gamma2*PD.rc_fac**alph

    if(process=='NESS_drift'):
        gam = PD.gamma2*(kp/kd)**alph - PD.gamma1*PD.r0

    if(process=='NESS_neg_two'):
        gam=PD.gamma2*PD.rc_fac**alph
    """

            

    ##########################################################################
    # Produce Data from Classes ###############################################
    ##########################################################################
if(bool_run):

    for params,name in zip(all_params,names):
        PD.run(params)
        save_name = "x_"+name+"highres10000.pickle"
        pickle.dump(PD.x_data,open(save_name,"wb"))    


