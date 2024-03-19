#include<stdio.h>
#include<stdlib.h>
#include "types.h"
#include "randomnumbers.h"
#include "init.h"
#include <math.h>

static void run(CE *ce);

CE* CE_create() {
	printf("Hallo");
	CE* tmp = (CE*) calloc((size_t)1,sizeof(CE));
	return tmp;
}

void CE_destroy(CE* ce) {
	if (ce != NULL) free(ce);
}

int main(){
	printf("Hallo"); 
	CE* ce = CE_create(); 
  	CE_init(ce);
 	 run(ce);
	CE_destroy(ce);
}





static void run(CE *ce){
        double t_hit = 0;
        double countNonTerm = 0.0;
        double countToShort = 0.0;
        //ce->ttt = np.linspace(0.0, ce->N*ce->dt, ce->N)

	for(int j=0; j < ce->nrealiz; j++){
	    for(int i=0; i < (ce->N-1); i++){
                // init distribution
		printf("Hallo");
                if (i== 0){
                    ce->x = (double *) calloc(ce->N,sizeof(double));
                    ce->x[i]= ce->sigma*gasdev() + ce->x_start_sim; 
		}
                // boundaries
                if (ce->x[i]<= ce->x0_for_reverse + ((i < (int)(1.0*ce->Tswitch/ce->dt)) ? 0.0 : 0.2 ) ){
                    // refelcting                    
                    if(i < (int)(1.0+1.0*ce->Tabsorb/ce->dt)){
                        ce->x[i] = ce->x0_for_reverse + fabs(ce->x0_for_reverse - ce->x[i]);
		    }
                    // absorbing    
                    else{
                        ce->x[i]= ce->x0_for_reverse;
                        t_hit = i;
                        if (ce->t_hit_min > t_hit) ce->t_hit_min = t_hit;            
                        if (t_hit > ce->t_hit_max) ce->t_hit_max = t_hit;   
                        break;
		    }
		}
                else{
                    // effectively value ice->last one due to +1  
                    if(i < (int)(1.0*ce->Tswitch/ce->dt)){
                        ce->x[i+1] = ce->x[i] - ce->gamma*(ce->x[i]- ce->x_start_sim)*ce->dt 
				      + 0.5*ce->D*ce->beta* pow(ce->x[i],(2.0*ce->beta-1.0)) 
				      + sqrt(ce->D*ce->dt)*pow(ce->x[i],ce->beta)*gasdev();
		    } else{
                        ce->x[i+1] = ce->x[i] - ce->gamma*pow(ce->x[i],- ce->alpha)*ce->dt 
				      + 0.5*ce->D*ce->beta* pow(ce->x[i],(2.0*ce->beta-1.0))
				      + sqrt(ce->D*ce->dt)*pow(ce->x[i],ce->beta)*gasdev();
		    }
		}

                      
                if(i == ce->N - 2){ 
		  t_hit = ce->N - 2;
		}

                  
            //"""figure(1)    
            //if(j <= 200):
            //    plot(ce->ttt,x)
            //    show()"""
		
		
	    free(ce->x);
	    }
	}
}
/*            
              # revert
            for i in range(int(t_hit)+1):
                # take care for never terminating ones
                if(t_hit == (self.N-2)):
                    y[i] = x[t_hit-i]*np.nan
                    if(i == t_hit): 
                        countNonTerm += 1.0
                        print "nonTerm"
                # take care of to short living oneself.     
                elif(t_hit <= int(1.0*self.TminLength /self.dt)):
                    y[i] = x[t_hit-i]*np.nan
                    if(i == t_hit): 
                        countToShort += 1.0
                        print "to short living"           
                else:
                    y[i] = x[t_hit-i]
                
                
            # self.m up mean and variances
                #take care of never terminating ones
                #take care of to short living ones
            if(~np.isnan(y[0])):
                for ii in range(t_hit+1):
                    ymean[ii] += y[ii]
                    yvar[ii] += y[ii]*y[ii]
                    yNorm[ii] += 1.0
            
            if(boolCov==1):
                # consider only trajectorieself.living longer than TminLength for covariance
                if(t_hit > int(1.0*self.TminLength /self.dt)):
                    z[j] = y[0:self.N/self.traFreq:self.wfreq]
                else:
                    z[j]= y[0:self.N/self.traFreq:self.wfreq] * np.nan
*/

/*                    
	}
	
	
	
	
	
	
	
        # calculate masked coovarinace matrix
        if(boolCov==1):
            zcovar = maCov.maskedCov(z)

        
        # output to screen
        print "number of non terminating trajectories"
        print countNonTerm
        print "T_min, fraction ------(t_hit_min*dt) (t_hit_min/N)"
        print t_hit_min*self.dt
        print 1.0*t_hit_min/self.N   
        print "T_max, fraction: ------ (t_hit_max*dt) (t_hit max / s.N)"
        print t_hit_max*self.dt
        maxLengthofTraj = (1.0/self.N)*t_hit_max
        print maxLengthofTraj
          
        # return evaluation of mean and variance
        for ii in range(self.N):
            ymean[ii] /= yNorm[ii]
            yvar[ii] = yvar[ii] / yNorm[ii] - ymean[ii] * ymean[ii]

        # reload plots            
        #self.ax2.relim()
        #self.ax2.autoscale_view()
        #self.fig2.canvas.draw()
        #plt.pause(1)  
        
        if(boolCov==1):
            return (ymean, yvar, t_hit_min, t_hit_max, zcovar )
        else:
            return (ymean, yvar, t_hit_min, t_hit_max)




*/
