#ifndef __TYPES_H__
#define __TYPES_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b)) 

typedef struct {
	double   alpha;
	double   xi;
	double   K;
	double   D1;
	double   D2;
	double   gamma1pre;
	double   gamma2pre;	
	double   gamma1;
	double   gamma2;
	double   A;
	int      Nb;
	double   kp1;
	double   kd1;
	double   kp;
	double   kd;
	double   c_start_sim;

	double   tau1;
	double   tau2;
	int      nrealiz;        
	int      seed;
	double   dt;
	double   x_start_sim;
	double   var_start_sim;
	double   x_start_mean;
	double   x0_for_reverse;
	double   x_upper_bound;
	int      total_time;
	int      max_rtime; 
	int      Tswitch;
	double   deltaTswitch;
	double   varDeltaTswitch;
	double   maxDeltaTswitch;
	int      TminLength;
	int      Tabsorb;
	int      wfreq;
	int      traFreq;
	int      N;
	double   sigma;
				
	double*  x;
	double*  c;
	int      t_hit_min;
	int      t_hit_max;
	         
	double*  xmean;
	double*  xmeanConc;
	double*  x2mean;
	double*  x3mean;
	int*     Norm;
	double** cov;
	int      t_rev_max;
	int      t_rev_max_wfreq;  
} CE;


#endif
