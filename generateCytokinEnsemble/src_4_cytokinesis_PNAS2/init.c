#include "types.h"
#include "init.h"
#include "randomnumbers.h"

void CE_init(CE *ce){
  
  //ce->gamma1           = 1.0 / ce->tau1			 ;
  //ce->gamma2           = 1.0 / ce->tau2			 ;
  ce->N               = (int) (ce->total_time / ce->dt)  	 ;
  ce->sigma           = sqrt(0.5*ce->D1 / (ce->gamma1pre))	 ;
  ce->t_hit_max       = 0   				 	 ; 
  ce->t_hit_min       = ce->N				 	 ;
  //ce->Tabsorb         = ce->Tswitch			 	 ;
  ce->t_rev_max       = (int) (ce->max_rtime / ce->dt)   	 ;
  if(ce->t_rev_max > ce->N) ce->t_rev_max = ce->N        	 ;
  ce->t_rev_max_wfreq = (int) (ce->t_rev_max/ce->wfreq)  	 ;
  
  initRandSeedStructs(ce->seed);
  
  ce->xmean  = (double  *) calloc( ce->t_rev_max		   ,sizeof(double));
  ce->xmeanConc  = (double  *) calloc( ce->t_rev_max		   ,sizeof(double));  
  ce->x2mean = (double  *) calloc( ce->t_rev_max		   ,sizeof(double));
  ce->x3mean = (double  *) calloc( ce->t_rev_max		   ,sizeof(double));  
  ce->Norm   = (int     *) calloc( ce->t_rev_max	 	   ,sizeof(int   ));
  ce->cov    = (double **) calloc( (int)(ce->t_rev_max_wfreq +1.0) ,sizeof(double));
  for(int i=0.0; i <= ce->t_rev_max_wfreq; i++){
    ce->cov[i] = (double *) calloc((int)(ce->t_rev_max_wfreq +1.0) ,sizeof(double));  
  }
  
}
