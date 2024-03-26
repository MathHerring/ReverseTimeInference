#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from __future__ import absolute_import, unicode_literals,division
import numpy as np
import os

###############################################################################
class generatePowerEnsemble(object):

    def __init__(self,fit_until_this_t=300,nrealiz=5,gamma=1.0,D=0.2,alpha=-1.0,beta=0.0):
        
        ##################################### these parameters are passed as argument
        self.fit_until_this_t = fit_until_this_t
        self.nrealiz=nrealiz
        self.gamma = gamma
        self.D = D
        self.alpha = alpha
        self.beta = beta

        ##################################### need to set these parameters in this class
        self.var_start_sim = 0
        self.x_upper_bound = 1000 
        self.max_rtime = fit_until_this_t
        self.traFreq = 1
        self.wfreq = 50
        self.total_time=2000
        self.x0_for_reverse = 0
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
fit_until_this_t=30
nrealiz = 2000

alpha_ar = [1]#[0.0, -1.0, -2.0]
gamma_ar = [1]#[0.03, 0.31, 2.16]
D_ar = [0.02]#[0.04,0.04,0.04]

for i in np.arange(len(alpha_ar)):
    alpha = alpha_ar[i]
    gamma = gamma_ar[i]
    D = D_ar[i]

    print('#####################################')
    print('alpha = ', alpha, ' ,gamma = ', gamma, ',D = ', D)
    print('#####################################')
    gPE = generatePowerEnsemble(fit_until_this_t,nrealiz,gamma,D,alpha)
    gPE.run()

    '''
    # load output of c++ simulation and save in python
    # can also load the raw trajectories and more, for info look into c code. the 'num_...' files are the c outputs. their naming should already be telling as well.
    md = np.loadtxt('num_mean.txt')
    mv = np.loadtxt('num_var.txt')
    cov_ar=np.loadtxt("num_cov.txt")
    #hittingTime=np.loadtxt("hittingTime.txt")
    
    np.savetxt("meannegtwo_gamma{}_D{}_alpha{}_beta{}_N{}.txt".format(int(gamma*10),int(D*10), int(alpha*10), int(beta*10), nrealiz), md)
    np.savetxt("varnegtwo_gamma{}_D{}_alpha{}_beta{}_N{}.txt".format(int(gamma*10),int(D*10), int(alpha*10), int(beta*10), nrealiz), mv)
    np.savetxt("covnegtwo_gamma{}_D{}_alpha{}_beta{}_N{}.txt".format(int(gamma*10),int(D*10), int(alpha*10), int(beta*10), nrealiz), cov_ar)
    '''


