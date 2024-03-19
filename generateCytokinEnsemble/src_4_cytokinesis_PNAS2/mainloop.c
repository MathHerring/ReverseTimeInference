#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "types.h"
#include "init.h"
#include "randomnumbers.h"

//#include "randomnumbers.h"

#define M_PI 3.14159265358979323846

#define PRINT_TRJ 1
#define PRINT_CONC 1
#define PRINT_TRJ_FWD 0
#define PRINT_HITTING_TIME 1
#define PRINT_DELTA_T_SWITCH 1
#define UPPER_BOUND 1
#define READ_XMEAN 0

#define SIGMOIDAL_SWITCH 0
#define OU_DRIFT 0
#define OU_FORCE_SWITCHOFF 0
#define PURE_DRIFT 0
#define CYTOKINESIS 0
#define CYTOKINESIS_FOR_FLY 1

CE* CE_create() {
        CE* tmp = (CE*) calloc((size_t)1,sizeof(CE));
        return tmp;
}

void CE_destroy(CE* ce) {
        if (ce != NULL) free(ce);
}

static void run(CE *ce){
	printf("start src_4\n");
  
        int t_hit = 0;
        double countNonTerm = 0.0;
        double countToShort = 0.0;
	double xmean2power;
	double var;
	double std;
	double delta_x;
	double gamma2;
	// to avoid problems with numerics -- reflection already at + secure_dist
	double secure_dist = 0.1;
	FILE *fileMean;
	FILE *fileMeanConc;	
	FILE *fileVar;
	FILE *fileCov;
	FILE *fileTrj;
	FILE *fileConc;	
	FILE *fileTrj_fwd;
	FILE *fileSkew;
	FILE *fileHittingTime;
	FILE *fileDeltaTswitch;
	
	fileMean = fopen("num_mean.txt", "w+");
	fileMeanConc = fopen("num_meanConc.txt", "w+");	
	fileVar  = fopen("num_var.txt" , "w+");
	fileCov  = fopen("num_cov.txt" , "w+");
	if(PRINT_TRJ)fileTrj = fopen("num_trj.txt", "w+");
	if(PRINT_CONC)fileConc = fopen("num_conc.txt", "w+");
	if(PRINT_TRJ_FWD)fileTrj_fwd = fopen("num_trj_fwd.txt", "w+");
	if(PRINT_HITTING_TIME) fileHittingTime = fopen("hittingTime.txt", "w+");
	if(PRINT_DELTA_T_SWITCH) fileDeltaTswitch = fopen("deltaTswitch.txt", "w+");
	// brauche freadf ....
	fileSkew = fopen("num_skew.txt", "w+");

	
	double* x_mean_read;
	FILE *fileReadXmean;
	
	
	
	if(READ_XMEAN){
	  x_mean_read  = (double  *) calloc( ce->nrealiz,sizeof(double));
	  //fileReadXmean  = fopen("start_values.txt" , "r+");
	  fileReadXmean  = fopen("start_values_qstat.txt" , "r+");
	  for(int j=0; j< ce->nrealiz; j++){
	    fscanf(fileReadXmean,"%lf\n",&(x_mean_read[j]) );
	    //printf("%f\n",x_mean_read[j]);
	  }
	}
	
	for(int j=0; j < ce->nrealiz; j++){
	    if((j%100)==0) printf("ensemble idx: %i\n",j);
	    if ((countNonTerm/ce->nrealiz) >= 0.9) {
	      printf("\n");
	      printf("Termiated as not collapsing\n");
	      break;
	    }
	    for(int i=0; i < (ce->N-1); i++){
                // init distribution
                if (i== 0){
                    ce->x = (double *) calloc(ce->N,sizeof(double));
		    ce->c = (double *) calloc(ce->N,sizeof(double));
                    ce->x[i]= ce->sigma*gasdev() + ce->x_start_sim;
		    gamma2 = ce->gamma2; //+ ( (double)rand() / (double)((unsigned)RAND_MAX + 1) )*1.0*ce->gamma2;
		    //gamma2 = ce->gamma2 + ( (double)rand() / (double)((unsigned)RAND_MAX + 1) )*1.0*ce->gamma2;
		    if(READ_XMEAN){
		      if(x_mean_read[j] <= 0.0){
		        printf("Exiting the program....\n");
		        printf("....invalid or to view init coords\n");			
			exit(0);
		      }
		      ce->x_start_mean = x_mean_read[j];
		      ce->x[i] = x_mean_read[j];
		    }else{
		      ce->x_start_mean = 0.0;
		      while(ce->x_start_mean <= 0.0){
			//printf("select starting position: %f\n",ce->x_start_mean);
			//ce->x_start_mean = sqrt(ce->var_start_sim)*gasdev() + ce->x_start_sim;
			ce->x_start_mean = ( (double)rand() / (double)((unsigned)RAND_MAX + 1) )*ce->var_start_sim + ce->x_start_sim;
			printf("select starting position: %f\n",ce->x[0]);
			ce->deltaTswitch = ce->maxDeltaTswitch - ( (double)rand() / (double)((unsigned)RAND_MAX + 1) )*ce->varDeltaTswitch;
			//printf("select deltaTswitch: %f  varDeltaTswitch:  %f\n",ce->deltaTswitch, ce->varDeltaTswitch);
		      }
		      //ce->x[i] = ce->x_start_mean;
		      ce->c[i] = ce->c_start_sim;
		    }
		}
                // boundaries
                // upper
                if (UPPER_BOUND && ce->x[i] >= ce->x_upper_bound)ce->x[i] = ce->x_upper_bound - fabs(ce->x[i] - ce->x_upper_bound);
		//  lower
		// old is plainly wrong:
                //if (ce->x[i]<= ce->x0_for_reverse + ((i < (int)(1.0*ce->Tswitch/ce->dt)) ? 0.0 : secure_dist ) ){
		//better
                //if (ce->x[i]<= ce->x0_for_reverse + ((i < (int)(1.0*ce->Tabsorb/ce->dt)) ? secure_dist : 0.0 ) ){
                // best -- secure_dist seems obsolet with multiplicative noise excludet ???
		if (ce->x[i]<= ce->x0_for_reverse){		  
		    //printf("first x[i] smaller 0:    %f\n",ce->x[i]);
		    //printf("should be zero:    %f\n",((i < (int)(1.0*ce->Tswitch/ce->dt)) ? 0.0 : secure_dist));
		    //printf("criterion:    %i\n",(int)(1.0*ce->Tswitch/ce->dt));
                    // refelcting                    
                    if(i < (int)(1.0+1.0*ce->Tabsorb/ce->dt)){
                        ce->x[i] = ce->x0_for_reverse + secure_dist + fabs(ce->x0_for_reverse - ce->x[i]);
		    }
                    // absorbing    
                    else{
		        //printf("trajNR : %i     x[t_hit-1]:   %f x[t_hit] :    %f\n",j,ce->x[i-1],ce->x[i]);
                        ce->x[i]= ce->x0_for_reverse;
 			//printf("x0 for reverse:    %f\n",ce->x0_for_reverse);
                        t_hit = i;
			//printf("trajNR : %i     t_hit : %i TabsorbInt:    %i \n",j,t_hit,(int)(1.0+1.0*ce->Tabsorb/ce->dt));
                        if (ce->t_hit_min > t_hit) ce->t_hit_min = t_hit;            
                        if (t_hit > ce->t_hit_max) ce->t_hit_max = t_hit;   
                        break;
		    }
		}

		if(CYTOKINESIS){
		    double switchFactor;
		    if((i==0) &&(j==0) ) printf("CYTOKINESIS\n");

		    // position
		    // switch constant A from off to on
		    if(ce->deltaTswitch > 0){
		      switchFactor = (i*ce->dt - ce->Tswitch)/ce->deltaTswitch;
		      switchFactor = 1.0 / (1.0 + exp(-switchFactor) );
		    }else{
		      if((i*ce->dt - ce->Tswitch) < 0 ) switchFactor = 0;
		      else switchFactor = 1;		      
		    }
		    //if((j==0) && ((i%ce->wfreq)==0) ) printf("%f      %f\n",switchFactor,ce->deltaTswitch);
		    // reduce c**2 to c
		    ce->x[i+1] = ce->x[i] 
		    - ce->xi * ( 
				ce->K*(ce->x[i]- ce->x_start_mean)
				+ 2.0 * M_PI * ce->A*switchFactor * ce->Nb * ce->c[i] // * ce->c[i]
				)* ce->dt
		    + sqrt(ce->D1 *ce->dt)*gasdev(); 
		    
		    // concentration
		    
		    ce->c[i+1] = ce->c[i]
		    + ( ce->kp - ce->kd * ce->c[i] ) * ce->dt 
		    - (ce->x[i+1]-ce->x[i])/(0.5*(ce->x[i+1] + ce->x[i]) ) * ce->c[i]
		    + sqrt(ce->D2 *ce->dt)*gasdev();
		    
		    //ce->c[i+1] = ce->kp / (ce->kd + (ce->x[i+1]-ce->x[i])/(0.5*ce->dt(ce->x[i+1] + ce->x[i]) )  );
		    //printf("trajNR : %i      x[t_hit] :    %f\n",j,ce->x[i]);

		}else if(CYTOKINESIS_FOR_FLY){
				    double switchFactor;
		    if((i==0) &&(j==0) ) printf("CYTOKINESIS FOR FLY -- PAPER PLOTS\n");

		    // position
		    // switch constant A from off to on
		    if(ce->deltaTswitch > 0){
		      switchFactor = (i*ce->dt - ce->Tswitch)/ce->deltaTswitch;
		      switchFactor = 1.0 / (1.0 + exp(-switchFactor) );
		    }else{
		      if((i*ce->dt - ce->Tswitch) < 0 ) switchFactor = 0;
		      else switchFactor = 1;
		    }
		    
		    // radius 
		    delta_x =   - ( ( (1.0-switchFactor)*ce->gamma1pre + switchFactor*ce->gamma1) * (ce->x[i]- ce->x_start_mean)
					+ ((1.0-switchFactor)*ce->gamma2pre + switchFactor*ce->gamma2) * pow(ce->c[i],ce->alpha) //* ce->c[i]
				  )* ce->dt
				+ sqrt(switchFactor*ce->D2 *ce->dt )*gasdev()
				+ sqrt((1.0-switchFactor)*ce->D1*ce->dt)*gasdev();
				
		    ce->x[i+1] = ce->x[i] + delta_x;
	    
		    // concentration
		    
		    ce->c[i+1] = ce->c[i]
		    + (1.0-switchFactor)*( ce->kp1 - ce->kd1 * ce->c[i] ) * ce->dt
		    + switchFactor*( ce->kp - ce->kd * ce->c[i] ) * ce->dt
		    - (delta_x / ce->x[i+1]) * ce->c[i]; // note that dt cancels with dt of delta_x
		    //+ sqrt( (ce->kp + ce->kd*ce->c[i])*ce->D2 *ce->dt )*gasdev();
		    
		}else{
		  if((i==0) &&(j==0) ) printf("NO FUNKTION\n");
		}

                      
                if(i == ce->N - 2){ 
		  t_hit = ce->N - 2;
		}	
	    }

	    // print hitting total_time
	    if(PRINT_HITTING_TIME){
	      int int_t_hit = (int)t_hit;
	      if(int_t_hit == (ce->N-2)){
		//fprintf(fileHittingTime,"%14.6f\n", 0.0/0);	
		fprintf(fileHittingTime,"\n");
	      }else{
		fprintf(fileHittingTime,"%14.6f\n",t_hit*ce->dt);
	      }
	    }
	    
	    // revert and average
	    int int_t_hit = (int)t_hit;
	    if(int_t_hit == (ce->N-2)){
		countNonTerm += 1.0;
		printf("nonTerm");
		printf("\n");
	    } else if(t_hit <= (int)(1.0*ce->TminLength /ce->dt)){
		countToShort += 1.0;
		printf("to short living");
	    } else{
		//int t_rev = ((int)(ce->t_rev_max - 1.0) < int_t_hit) ? (ce->t_rev_max -1.0) : int_t_hit;
		int t_rev = ((int)(ce->t_rev_max) < int_t_hit) ? ce->t_rev_max : int_t_hit;		
		// --1-- covariance
		for(int ii = 0.0; ii <= t_rev; ii+= ce->wfreq){
		  for(int jj = 0.0; jj <= t_rev; jj+= ce->wfreq){
		      // remember that norm is for state n-1
		      // norm zero means no partner at all -- so no cov contribution
		      //ce->x[int_t_hit-j] == 0.0 means either start !!! -- obsolet as only contributing parts are in foor loop
		      //if(!(ce->x[int_t_hit-i] == 0.0 || ce->x[int_t_hit-j] == 0.0) && !(ce->Norm[i] == 0.0 || ce->Norm[j] == 0.0)){
		      if(!(ce->Norm[ii] == 0.0 || ce->Norm[jj] == 0.0)){
		      // prefactor norm[j] has to be tce smaller of both
			  int norm = min(ce->Norm[ii],ce->Norm[jj]);
			  int iidx = (int)(ii/ce->wfreq);
			  int jidx = (int)(jj/ce->wfreq);
			  double XX = ( ce->x[int_t_hit-ii] - (ce->xmean[ii] / ce->Norm[ii]) );
			  double YY = ( ce->x[int_t_hit-jj] - (ce->xmean[jj] / ce->Norm[jj]) );
			  //printf("XX %f \n", XX);
			  //printf("YY %f \n", YY);
			  ce->cov[iidx][jidx]      += XX*YY*(norm) /(norm+1.0); 
			  //printf(" nrealiz: %i norm:  %i covariance:   %f\n",j, norm, ce->cov[iidx][jidx]);
		      }
		  }
		}
		
		// --2-- mean and variance

		// print forward trajectory
		if(PRINT_TRJ_FWD){
		  int tag_first_zero = 0;
		  // for(int ii=0; ii<= (int)(ce->t_rev_max/ce->traFreq); ii++){//(int)(ce->N/ce->traFreq); ii++){
		  for(int ii=0; ii<= (int)(ce->N/ce->traFreq); ii++){
		    if(ce->x[ii*ce->traFreq] > 0.0) fprintf(fileTrj_fwd,"%14.6f",ce->x[ii*ce->traFreq]);
		    else if(tag_first_zero == 0 ){ 
		      fprintf(fileTrj_fwd,"%14.6f", 0.0);
		      tag_first_zero = 1;
		    }
		    else fprintf(fileTrj_fwd,"%14.6f", 0.0/0);
		  }
		  fprintf(fileTrj_fwd,"\n");
		}				
		
		
		if(PRINT_TRJ){
		  for(int ii=0; ii<= (int)(ce->t_rev_max/ce->traFreq); ii++){
		    if( (int)(ii*ce->traFreq)<= t_rev ) fprintf(fileTrj,"%14.6f",ce->x[int_t_hit-ii*ce->traFreq]);
		    else fprintf(fileTrj,"%14.6f", 0.0/0);		    
		  }
		    if(PRINT_DELTA_T_SWITCH) fprintf(fileDeltaTswitch,"%14.6f\n",ce->deltaTswitch);
		}
		
		
		if(PRINT_CONC){
		  for(int ii=0; ii<= (int)(ce->t_rev_max/ce->traFreq); ii++){
		    if( (int)(ii*ce->traFreq)<= t_rev ) fprintf(fileConc,"%14.6f",ce->c[int_t_hit-ii*ce->traFreq]);
		    else fprintf(fileConc,"%14.6f", 0.0/0);
		  }
		}		
		
		for(int ii = 0.0; ii <= t_rev; ii++){
		  //if(PRINT_TRJ) fprintf(fileTrj,"%14.6f",ce->x[int_t_hit-ii]);
		  ce->xmean [ii]    += ce->x[int_t_hit-ii];
		  ce->xmeanConc[ii]    += ce->c[int_t_hit-ii];
		  ce->x2mean[ii]    += ce->x[int_t_hit-ii]*ce->x[int_t_hit-ii];
		  ce->x3mean[ii]    += ce->x[int_t_hit-ii]*ce->x[int_t_hit-ii]*ce->x[int_t_hit-ii];
		  ce->Norm  [ii]    += 1.0;
		}
	    }
	    if(PRINT_TRJ) fprintf(fileTrj,"\n");
	    free(ce->x);
	    if(PRINT_CONC) fprintf(fileConc,"\n");
	    free(ce->c);
	}
          
        // Normalization of mean and variance    
        for(int ii=0; ii< ce->t_rev_max; ii++){
	  if(ce->Norm[ii]  == 0.0){
	    ce->xmean[ii]   = 0.0;
	    ce->xmeanConc[ii]   = 0.0;
	    ce->x2mean[ii]  = 0.0;
	    ce->x3mean[ii]  = 0.0;
	  } else{
	    ce->xmean[ii]  /= ce->Norm[ii];
	    ce->xmeanConc[ii]  /= ce->Norm[ii];
	    ce->x2mean[ii] /= ce->Norm[ii];
	    ce->x3mean[ii] /= ce->Norm[ii];
	    
	    xmean2power     = ce->xmean[ii]*ce->xmean[ii];
	    var 	    = ce->x2mean[ii]-xmean2power;
	    std 	    = sqrt(var);
	    ce->x3mean[ii]  = (ce->x3mean[ii]-3.0*ce->x2mean[ii]*ce->xmean[ii] + 2.0*xmean2power*ce->xmean[ii]) / (std*std*std);    
	    
	    ce->x2mean[ii]  =  var;
	  }
	}
	
		  
	// Normalization of covariance
	//int t_rev_max_wfreq = (int)(ce->t_rev_max/ce->wfreq);
	for(int ii = 0.0; ii < ce->t_rev_max_wfreq; ii++){
	  for(int jj = 0.0; jj < ce->t_rev_max_wfreq; jj++){
	    int idx = ii * ce->wfreq;
	    int jdx = jj * ce->wfreq;
	    int norm = min(ce->Norm[idx],ce->Norm[jdx]);
	    if(norm == 0.0){
	      ce->cov[ii][jj] = 0.0;
	    } else{
	      ce->cov[ii][jj]/=norm;
	    }
	  }
	}
	
	// consitency check
	/*
	for(int ii = 0.0; ii < ce->t_rev_max_wfreq; ii++){
	  if (ii%ce->wfreq == 0.0){
	    int idx = (int) (ii * ce->wfreq);
	    printf("variance:  %f	cov:   %f \n", ce->x2mean[idx] , ce->cov[ii][ii]);
	  }
	}
	*/
	  
	  
	
	// Output
	
	for(int ii=0; ii< ce->t_rev_max; ii++){
	  //if( (ii>0) && (ce->xmean[ii]) <= 0.0 ) fprintf(fileMean,"%14.6f", 0.0/0);
	  //else fprintf(fileMean,"%f\n",ce->xmean[ii]);
	  fprintf(fileMean,"%f\n",ce->xmean[ii]);
	}
	for(int ii=0; ii< ce->t_rev_max; ii++){
	  //if( (ii>0) && (ce->xmeanConc[ii]) <= 0.0 ) fprintf(fileMeanConc,"%14.6f", 0.0/0);
	  //else fprintf(fileMeanConc,"%f\n",ce->xmeanConc[ii]);	  
	  fprintf(fileMeanConc,"%f\n",ce->xmeanConc[ii]);
	}		
	for(int ii=0; ii< ce->t_rev_max; ii++){
	  fprintf(fileVar,"%f\n",ce->x2mean[ii]);
	}
	for(int ii=0; ii< ce->t_rev_max; ii++){
	  fprintf(fileSkew,"%f\n",ce->x3mean[ii]);
	}	
	for(int ii = 0.0; ii < ce->t_rev_max_wfreq; ii++){
	  for(int jj = 0.0; jj < ce->t_rev_max_wfreq; jj++){    
	    fprintf(fileCov,"%f\n",ce->cov[ii][jj]);
	  }
	}
	
	fclose(fileMean);
	fclose(fileVar);
	fclose(fileCov);
	fclose(fileMeanConc);
	if(PRINT_TRJ )fclose(fileTrj);
	if(PRINT_CONC )fclose(fileConc);
	if(PRINT_TRJ_FWD)fclose(fileTrj_fwd);
	if(PRINT_HITTING_TIME)fclose(fileHittingTime);
	if(PRINT_DELTA_T_SWITCH)fclose(fileDeltaTswitch);
	fclose(fileSkew);
	if(READ_XMEAN)fclose(fileReadXmean);
}

int main(int argc, char** argv){
	CE* ce = CE_create();
        // defaults

	/*
	ce->alpha          = 0.3 ; 
	ce->beta	   = 0.0 ; 
	ce->D1 	       	   = 0.7 ;
	ce->D2		   = 0.7 ;
	ce->gamma1	   = 1.0 ;
	ce->gamma2         = 1.0 ;
	ce->nrealiz	   = 2000;
	ce->seed           = 123 ;  
	ce->dt 	       	   = 1e-2;
	ce->x_start_sim    = 5.2 ;
	ce->x0_for_reverse = 0.0 ;
	ce->total_time     = 200 ;   
	ce->Tswitch        = 0   ;
	ce->deltaTswitch   = 0	 ;
	ce->TminLength     = 0   ;
	ce->wfreq 	   = 5   ;
	ce->traFreq        = 1   ;
	ce->max_rtime      = ce->total_time ;
	ce->x_upper_bound  = 10000;
	ce->Tabsorb        = ce->Tswitch;
	ce->var_start_sim  = 0.0;
	*/
      
	switch (argc) {
	    default:
	    case 34: ce->gamma2pre	 = atof(argv[33]);
	    case 33: ce->gamma1pre	 = atof(argv[32]);
	    case 32: ce->kd1		 = atof(argv[31]);
	    case 31: ce->kp1		 = atof(argv[30]);	    
	      
	    case 30: ce->alpha		 = atof(argv[29]);
	    case 29: ce->gamma2		 = atof(argv[28]);
	    case 28: ce->gamma1 	 = atof(argv[27]);	      
	    case 27: ce->maxDeltaTswitch = atof(argv[26]);
	    case 26: ce->varDeltaTswitch = atof(argv[25]);
	    case 25: ce->var_start_sim   = atof(argv[24]);
	    case 24: ce->Tabsorb	 = atoi(argv[23]);
	    case 23: ce->x_upper_bound   = atof(argv[22]);
	    case 22: ce->max_rtime	 = atoi(argv[21]);            
	    case 21: ce->traFreq 	 = atoi(argv[20]); 
	    case 20: ce->wfreq	 	 = atoi(argv[19]);
	    case 19: ce->TminLength 	 = atoi(argv[18]);
	    case 18: ce->deltaTswitch	 = atof(argv[17]);
	    case 17: ce->Tswitch 	 = atoi(argv[16]);
	    case 16: ce->total_time 	 = atoi(argv[15]);
	    case 15: ce->x0_for_reverse  = atof(argv[14]);

	    case 14: ce->D2   		 = atof(argv[13]);
	    case 13: ce->D1   		 = atof(argv[12]);
	    

	    case 12 : ce->dt 		 = atof(argv[11]);
	    case 11: ce->seed 		 = atoi(argv[10]);
	    case 10: ce->nrealiz 	 = atoi(argv[9 ]);

	    case 9 : ce->c_start_sim	 = atof(argv[8 ]);	    
	    
	    case 8 : ce->kd		 = atof(argv[7 ]);
	    case 7 : ce->kp		 = atof(argv[6 ]);
	    
	    case 6 : ce->Nb		 = atoi(argv[5 ]);
	    case 5 : ce->A		 = atof(argv[4 ]);

	    case 4 : ce->x_start_sim 	 = atof(argv[3 ]);	    
	    
	    case 3 : ce->K 	 	 = atof(argv[2 ]);
	    case 2 : ce->xi 		 = atof(argv[1 ]);
	    case 1 : break;
	}

	
	    printf("# # # # P A R A M E T E R S # # # # # \n");
	    printf("gamma2pre        %f\n",ce->gamma2pre      	 );	    
	    printf("gamma1pre        %f\n",ce->gamma1pre      	 );	    
	    printf("kd1              %f\n",ce->kd1      	 );
	    printf("kp1              %f\n",ce->kp1      	 );	    
	    printf("alpha            %f\n",ce->alpha      	 );
	    printf("gamma2           %f\n",ce->gamma2      	 );
	    printf("gamma1           %f\n",ce->gamma1         	 );
	    printf("varDeltaTswitch  %f\n",ce->varDeltaTswitch   );
	    printf("maxDeltaTswitch  %f\n",ce->maxDeltaTswitch   );
	    printf("var_start_sim    %f\n",ce->var_start_sim     );
	    printf("Tabsorb          %i\n",ce->Tabsorb	         );
	    printf("x_upper_bound    %f\n",ce->x_upper_bound     );
	    printf("max_rtime        %i\n",ce->max_rtime         );
	    printf("traFreq          %i\n",ce->traFreq 	         );
	    printf("wfreq            %i\n",ce->wfreq	         );
	    printf("TminLength       %i\n",ce->TminLength        );
	    printf("deltaTswitch     %f\n",ce->deltaTswitch      );
	    printf("Tswitch          %i\n",ce->Tswitch 	         );
	    printf("total_time       %i\n",ce->total_time        );
	    printf("x0_for_reverse   %f\n",ce->x0_for_reverse    );
	    printf("D2               %f\n",ce->D2   	         );
	    printf("D1               %f\n",ce->D1   	         );
	    printf("dt               %f\n", ce->dt 	         );
	    printf("seed             %i\n",ce->seed 	         );
	    printf("nrealiz          %i\n",ce->nrealiz 	         );
	    printf("c_start_sim      %f\n",ce->c_start_sim       );
	    printf("kd               %f\n",ce->kd		 );
	    printf("kp               %f\n",ce->kp		 );
	    printf("Nb               %i\n",ce->Nb	         );
	    printf("A                %f\n",ce->A		 );
	    printf("x_start_sim      %f\n",ce->x_start_sim       );
	    printf("K                %f\n",ce->K 	 	 );
	    printf("xi               %f\n",ce->xi 		 );
	    printf("============================================== \n");

	
	CE_init(ce);
	run(ce);
	            
        // free mean, var and cov    
        free(ce->xmean);
        free(ce->xmeanConc);	
	free(ce->x2mean);
	free(ce->x3mean);
	free(ce->Norm);
	
	/*FIXME memory to big !!!!!!1*/
	//int t_rev_max_wfreq = (int)(ce->t_rev_max/ce->wfreq);
	for(int i=0; i<= ce->t_rev_max_wfreq; i++) free(ce->cov[i]);
	free(ce->cov);
	
	CE_destroy(ce);
}




