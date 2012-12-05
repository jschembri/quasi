//File for the 1D Quasi Code
//Jeremy Schembri on Nov 8, 2012
#include <iostream>
#include <math.h>
#include <string>
#include <cstdlib> 
#include <fstream>
#include "constants.h"

using namespace std;

int main(int argc, char **argv){

   double endtime = atof(argv[1]);
   double time = 0;
   double row[x_spaces+1];
   double velocity[x_spaces+1];
   double energy[x_spaces+1];
   double epsilon[x_spaces+1];
   double pressure[x_spaces+1];
   double Mach[x_spaces+1];
   double areas[x_spaces+1];
   double alpha = 0.1;
   double e0, FplusHalf,FminusHalf;

   e0 = P0/(fluid_gamma-1.0) + row0*pow(u0,2)/2.0;


   // Creating x value
   double x_value[x_spaces+1];
   for (int i=0; i<=x_spaces; i++){
      x_value[i] = x_lower + delta_x*i;
   }

   // U Vector
	double Uplus1[x_spaces+1][3];
   double  U[x_spaces+1][3];
   double  F[x_spaces+1][3];

   for (int i=0;i<=x_spaces;i++){
         U[i][0] = row0;
         U[i][1] = row0*u0;
         U[i][2] = e0; 

         F[i][0] = row0*u0;
         F[i][1] = row0*pow(u0,2) + P0;
         F[i][2] = (e0+P0)*u0;

   }


// loop and work

	double volumes[x_spaces+1];
	for(int i=0;i<=x_spaces;i++){
		if (i == 0) {
			volumes[i] = delta_x*(area(x_value[i])+ delta_x/2.0);
		}else if (i==x_spaces){
			volumes[i] = delta_x*(area(x_value[i])- delta_x/2.0);
		}else{
			volumes[i] = delta_x/2.0*(area(x_value[i] + delta_x/2.0) + area(x_value[i] - delta_x/2.0));
		}
	}


  // Q Vector Initialization
   double  Q[x_spaces+1][3];
   for (int i=0;i<=x_spaces;i++){
         Q[i][0] = 0;
			if (i==0){
		   	Q[i][1] = P0*(area(x_value[i] + delta_x/2.0) - area(x_value[i])); 
			}else if (i==x_spaces){
		   	Q[i][1] = P0*(area(x_value[i]) - area(x_value[i]-delta_x/2.0)); 
			}else{
		   	Q[i][1] = P0*(area(x_value[i] + delta_x/2.0) - area(x_value[i] - delta_x/2.0)); 
			}

         Q[i][2] = 0;
   }

	//Setting up if time ==0
	for (int i=0; i<=x_spaces; i++){
		row[i]= U[i][0];
		velocity[i] = U[i][1]/row[i];
		energy[i] = U[i][2];
		epsilon[i] = energy[i]/row[i] - pow(velocity[i],2)/2.0;
		pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
	}



	while(time < endtime){	
	// Finite Volume analysis
		for (int i=1; i<=x_spaces; i++){
		   for (int j=0; j<=2; j++){
				if (i==x_spaces){
					FplusHalf = 0.5*(F[i][j] + F[i][j])-0.5*alpha*(velocity[i]-velocity[i]);
				}else{
					FplusHalf = 0.5*(F[i+1][j] + F[i][j])-0.5*alpha*(velocity[i+1]-velocity[i]);
				}
				FminusHalf = 0.5*(F[i][j] + F[i-1][j])-0.5*(velocity[i]-velocity[i-1]);
				if (i==x_spaces){
 					Uplus1[i][j] = U[i][j] - delta_t/volumes[i]*(FplusHalf*area(x_value[i])- FplusHalf*area(x_value[i] -delta_x/2.0))+delta_t/volumes[i]*Q[i][j];
				}else{
 					Uplus1[i][j] = U[i][j] - delta_t/volumes[i]*(FplusHalf*area(x_value[i]+delta_x/2.0)- FminusHalf*area(x_value[i] -delta_x/2.0))+delta_t/volumes[i]*Q[i][j];
				}
		          
		   }
		}

		for (int i=1; i<=x_spaces; i++){
			row[i]= Uplus1[i][0];
			velocity[i] = Uplus1[i][1]/row[i];
			energy[i] = Uplus1[i][2];
			epsilon[i] = energy[i]/row[i] - pow(velocity[i],2)/2.0;
			pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
		}

		for(int i=1;i<=x_spaces;i++){
			for(int j=0; j<=2;j++){
				U[i][j] = Uplus1[i][j];
			}
			F[i][0] = row[i]*velocity[i];
			F[i][1] = row[i]*pow(velocity[i],2) + pressure[i];
			F[i][2] = (energy[i]+pressure[i])*velocity[i];

			Q[i][0] = 0;
			if (i==x_spaces){
				Q[i][1] = pressure[i]*(area(x_value[i]) - area(x_value[i] - delta_x/2.0)); 
			}else{
				Q[i][1] = pressure[i]*(area(x_value[i] + delta_x/2.0) - area(x_value[i] - delta_x/2.0)); 
			}
			Q[i][2] = 0;

		}
		time += delta_t;
	}


for (int i=0; i<=x_spaces; i++){
	Mach[i] = velocity[i] / pow(fluid_gamma*pressure[i]/row[i],0.5);
}
	for(int i=0;i<=x_spaces;i++){
		areas[i] = area(x_value[i]);
	}

 printarray (x_value,x_spaces+1, "X Value");
 printarray (row,x_spaces+1, "Y Value");
 printarray (areas,x_spaces+1, "Area");
 printarray (velocity,x_spaces+1, "Velocity");
 printarray (pressure,x_spaces+1, "Pressure");
 printarray (energy,x_spaces+1, "Energy");
 printarray (volumes,x_spaces+1, "Volumes");
 printarray (Mach,x_spaces+1, "Mach");




	return 0; 


}



