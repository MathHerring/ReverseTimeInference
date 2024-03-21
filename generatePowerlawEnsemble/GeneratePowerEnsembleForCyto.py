#!/usr/bin/env python2
# -*- coding: utf-8 -*-


from __future__ import absolute_import, unicode_literals,division#, print_function
import numpy as np
import scipy as sp

import os

###############################################################################


## need to set parameters in this class
class PureFwdPowerlaw(object):

    def __init__(self,fit_until_this_t=300,nrealiz=5,gamma=1.0,D=0.2,alpha=-1.0,beta=0.0):


        self.fit_until_this_t = fit_until_this_t
        self.nrealiz=nrealiz
        self.gamma = gamma
        self.D = D
        self.alpha = alpha
        self.beta = beta



        #####################################
        self.var_start_sim = 0
        self.x_upper_bound = 1000 #3
        self.max_rtime = fit_until_this_t
        self.traFreq = 1#100
        self.wfreq = 50#100#100
        self.total_time=2000#80#max(40,fit_until_this_t)
        self.x0_for_reverse = 0
        #self.D = 0.2
        self.dt = 0.005
        self.dt_eff = self.dt*self.traFreq
        self.seed = 12345
        self.x_start_sim = 20


    def run(self):

        ###########################################      
        os.system('./src_5_pureFwdPowerlaw_notrj/prog %f %d %d %f %f %f %d %d %d %d %f %f %f %f %f' % (
        self.x_start_sim,self.nrealiz,self.seed,self.dt,self.D,
        self.x0_for_reverse,self.total_time,self.wfreq,
        self.traFreq,self.max_rtime,self.x_upper_bound,
        self.var_start_sim,self.gamma,self.alpha,self.beta)
        )
        ###########################################
        
         
##############################################################################
#gamma,D,alpha,beta
params=[1.0,0.2,-1.0,0.0]
fit_until_this_t=30
fit_until_this_t_ref = fit_until_this_t
nrealiz = 2000
gamma,D,alpha,beta = params
max_load_fwd = 50
t_theory = np.linspace(0,3,300)
x0 = 1 ### ATTENTION: must be set above in the c++ simulation, this one is only for saving the data and must be set manually

#gamma = 0.16228566887687898

simulationkey = True 
evaluationkey = False 


if simulationkey == True:
    alpha_ar = [1]#[0.0, -1.0, -2.0]
    gamma_ar = [1]#[0.03, 0.31, 2.16]
    D_ar = [0.02]#[0.04,0.04,0.04]

    for i in np.arange(len(alpha_ar)):
        alpha = alpha_ar[i]
        gamma = gamma_ar[i]
        D = D_ar[i]

        #gamma = 0.0
        print('#####################################')
        print('alpha = ', alpha, ' ,gamma = ', gamma, ',D = ', D)
        print('#####################################')
        PFP = PureFwdPowerlaw(fit_until_this_t,nrealiz,gamma,D,alpha,beta)
        PFP.run()

        '''
        # load output of c++ simulation and save in python
        # can also load the raw trajectories and more, for info look into c++ code. the 'num_...' files are the c++ outputs. their naming should already be telling as well.
        md = np.loadtxt('num_mean.txt')
        mv = np.loadtxt('num_var.txt')
        cov_ar=np.loadtxt("num_cov.txt")
        #hittingTime=np.loadtxt("hittingTime.txt")
        
        np.savetxt("meannegtwo_gamma{}_D{}_alpha{}_beta{}_N{}.txt".format(int(gamma*10),int(D*10), int(alpha*10), int(beta*10), nrealiz), md)
        np.savetxt("varnegtwo_gamma{}_D{}_alpha{}_beta{}_N{}.txt".format(int(gamma*10),int(D*10), int(alpha*10), int(beta*10), nrealiz), mv)
        np.savetxt("covnegtwo_gamma{}_D{}_alpha{}_beta{}_N{}.txt".format(int(gamma*10),int(D*10), int(alpha*10), int(beta*10), nrealiz), cov_ar)
        '''


