/**************************************************************************************************************** 
C++ version - of SDEquation RMA model, to verify is structural sensitivity still works in a stochastic context
FB 11/04/2020 
*****************************************************************************************************************/ 

//libraries
#include<iostream>
#include<fstream>
#include<math.h>
#include<cassert>
#include<cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

using namespace std;

//Parameters
#define PI 3.14159265
#define SEED 17
#define INTERVALS 20
#define EVERYX 25

const float r=1.0;
const float m=0.1;
const float epsilon=0.3;
// Holling FR
const float a_H = 3.05;
const float b_H = 2.68;
// Ivlev FR
const float a_I = 1.0;
const float b_I = 2.0;
// Tangent FR
const float a_T = 0.99;
const float b_T = 1.48;
// Stochasticity 
const float sigmaE = 0.25;


//GENERAL STORAGE PARAMETERS
const float TimeMax=800.0; //maximal size of array for densities
const float dt=1.0/INTERVALS; //time increment to monitor the size of the populations
const int length_stored=500; //length of time (as integer) stored
const int klength=50; // length of vector for carrying capacity
const float K_min = 0.1;
const float K_max = 1.5;

// DATA STRUCTURES TO STORE RESULTS
// no need

const gsl_rng_type * T;
gsl_rng *rgsl;  		// random number generator - better be global TO IMPROVE

// FUNCTIONS (if any)


int main(){
// INITIALIZATION OF GNU SCIENTIFIC LIBRARY

    gsl_rng_env_setup();
    T=gsl_rng_default;
    rgsl = gsl_rng_alloc (T);     // pick random number generator
    int seed = SEED;
    gsl_rng_set (rgsl, seed);                  // set seed


    float t,K; // variables
    float n,p,n2,p2,n_I,p_I,n2_I,p2_I,n_T,p_T,n2_T,p2_T; 
    /*float N[length_stored],P[length_stored],N2[length_stored],P2[length_stored];//not needed for now. 
    //float N_I[length_stored],P_I[length_stored],N2_I[length_stored],P2_I[length_stored];
    //float N_T[length_stored],P_T[length_stored],N2_T[length_stored],P2_T[length_stored]; */
    int ti,k_index; 
    float PreyGrowth,PredGrowth,Predation,PreyGrowth2,PredGrowth2,Predation2; // variables
    float PreyGrowth_I,PredGrowth_I,Predation_I,PreyGrowth2_I,PredGrowth2_I,Predation2_I; // variables
    float PreyGrowth_T,PredGrowth_T,Predation_T,PreyGrowth2_T,PredGrowth2_T,Predation2_T; // variables
    float StochCompPrey,StochCompPred,StochCompPrey_I,StochCompPred_I,StochCompPrey_T,StochCompPred_T; 
    float time_for_saving_density; 

    // INPUT - OUPUT
    ofstream f1,f2,f3; // declare files -- had trouble with f[3] for some reason
    f1.open("RMA.txt");
    f2.open("RMA_I.txt");
    f3.open("RMA_T.txt");

for (k_index=1;k_index<klength;k_index++)
	{// loop on K
          K = K_min + (K_max - K_min)*k_index/klength;  //Defining K the carrying capacity
	
	  //Initial values
          n = 0.5;
          p = 0.5;
          n2 = 0.5;
          p2 = 0.5;

	  n_I = 0.5;
	  p_I = 0.5;
	  n2_I = 0.5;
	  p2_I = 0.5;

	  n_T = 0.5;
	  p_T = 0.5;
	  n2_T = 0.5;
	  p2_T = 0.5;
	  
        t=0.0;ti = 1; time_for_saving_density = TimeMax*INTERVALS-length_stored*EVERYX+1;
	while (t<TimeMax){
		
		// Record densities for a number of times and all K values
		if ((ti>time_for_saving_density) &(ti % EVERYX == 0) ){
			f1<<n<<" "<<p<<" "<<n2<<" "<<p2<<" "<<t<<" "<<K<<endl; // these would just be over time
			f2<<n_I<<" "<<p_I<<" "<<n2_I<<" "<<p2_I<<" "<<t<<" "<<K<<endl; 
			f3<<n_T<<" "<<p_T<<" "<<n2_T<<" "<<p2_T<<" "<<t<<" "<<K<<endl; 	
			
		}

		// Compute growth and predation processes for the various models
		PreyGrowth = r*n*(1.0-n/K);
		Predation = (a_H*n*p)/(b_H*n+1.0);
		PredGrowth = epsilon*Predation-m*p;
	
		PreyGrowth2 = r*n2*(1.0-n2/K);
		Predation2 = (a_H*n2*p2)/(b_H*n2+1.0);
		PredGrowth2 = epsilon*Predation2-m*p2;

		PreyGrowth_I = r*n_I*(1.0-n_I/K);
		Predation_I = a_I*(1.0-exp(-b_I*n_I))*p_I;
		PredGrowth_I = epsilon*Predation_I-m*p_I;

		PreyGrowth2_I = r*n2_I*(1.0-n2_I/K);
		Predation2_I = a_I*(1.0-exp(-b_I*n2_I))*p2_I;
		PredGrowth2_I = epsilon*Predation2_I-m*p2_I;

		PreyGrowth_T = r*n_T*(1.0-n_T/K);
		Predation_T = a_T*tanh(b_T*n_T)*p_T;
		PredGrowth_T = epsilon*Predation_T-m*p_T;

		PreyGrowth2_T = r*n2_T*(1.0-n2_T/K);
		Predation2_T =  a_T*tanh(b_T*n2_T)*p2_T;
		PredGrowth2_T = epsilon*Predation2_T-m*p2_T;

		// Add stochasticity
		// This is Euler-Maruyama (simplest) scheme 
		// Multiplication by abundance for environmental stochasticity 
		StochCompPrey = sigmaE*n*sqrt(dt)*gsl_ran_gaussian(rgsl,1.0); 
		StochCompPred = sigmaE*p*sqrt(dt)*gsl_ran_gaussian(rgsl,1.0); 
		StochCompPrey_I = sigmaE*n_I*sqrt(dt)*gsl_ran_gaussian(rgsl,1.0); 
		StochCompPred_I = sigmaE*p_I*sqrt(dt)*gsl_ran_gaussian(rgsl,1.0);
		StochCompPrey_T = sigmaE*n_T*sqrt(dt)*gsl_ran_gaussian(rgsl,1.0); 
		StochCompPred_T = sigmaE*p_T*sqrt(dt)*gsl_ran_gaussian(rgsl,1.0);
		
		// Pop growth
		n = n + dt*(PreyGrowth - Predation) + StochCompPrey; //+ DemogCompPrey; // Only env. stoch. on the prey
		p = p + dt*PredGrowth + StochCompPred ; // + DemogCompPred Put only dem. stoch. on the predator, ideally
		//deterministic equivalent
		n2 = n2 + dt*(PreyGrowth2 - Predation2); 
		p2 = p2 + dt*PredGrowth2; 
		//Ivlev
		n_I = n_I + dt*(PreyGrowth_I - Predation_I) + StochCompPrey_I;
		p_I = p_I + dt*PredGrowth_I + StochCompPred_I ; 
		n2_I = n2_I + dt*(PreyGrowth2_I - Predation2_I); 
		p2_I = p2_I + dt*PredGrowth2_I; 
		// Arctan
		n_T = n_T + dt*(PreyGrowth_T - Predation_T) + StochCompPrey_T;
		p_T = p_T + dt*PredGrowth_T + StochCompPred_T; 
		n2_T = n2_T + dt*(PreyGrowth2_T - Predation2_T);  
		p2_T = p2_T + dt*PredGrowth2_T; 
			
		//if ( (ti % 100 ) == 0){
		//cout<<" t =  "<<t<<" Iteration "<<ti<<endl; // to screen
		//}
		t=t+dt;ti++; 

		} // end of loop on time
cout<<" K=  "<<K<<endl; // to screen -- say where we are in the loop. 

}// end of loop on K

//for (int i=1; i<=3; i++){f[i].close();};//closing files, for some reason f[i] seems to be causing trouble
f1.close(); f2.close(); f3.close(); 

return 0;

}
