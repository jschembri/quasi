
//File for the 1D Quasi Code
//Jeremy Schembri on Nov 8, 2012
#include <iostream>
#include <math.h>
#include <cmath>
#include <string>
#include <cstdlib> 
#include <fstream>
#include "constants.h"

using namespace std;

double max_residual( double guess[],double answer[], int length){
	double max;
	for(int i=0;i<length;i++){
		if (i==0){
			max = abs(answer[i]-guess[i] );
		}

		if (abs(answer[i]-guess[i]) > max){
			max = abs(answer[i]-guess[i] );
		}
	}
	return max;
}



int main(int argc, char **argv){

	//Setting up values from 0 to i+1
   double endtime = atof(argv[1]);
   double time = 0;
   double row[x_spaces+3];
   double velocity[x_spaces+3];
   double energy[x_spaces+3];
   double epsilon[x_spaces+3];
   double pressure[x_spaces+3];
   double Mach[x_spaces+1];
   double areas[x_spaces+3];
   double alpha = 80;
   double e0, FplusHalf,FminusHalf;

   double temp_Mach[x_spaces+3];
	int iteration = 0;
	//int max_iterations = endtime/delta_t;
	//double iteration_list[max_iterations+1];
	//double residual_list[max_iterations+1];

	double uend =  Mstart*pow(fluid_gamma*Pend/row0,0.5);
   e0 = Pend/(fluid_gamma-1.0) + row0*pow(uend,2)/2.0;
	

   // Creating x value
   double x_value[x_spaces+3];
   for (int i=0; i<=x_spaces+2; i++){
		if (i==x_spaces+2){
			x_value[i] = x_value[i-1];
		}else if (i==0){
			x_value[i] = x_lower;
		}else{
      	x_value[i] = x_lower + delta_x*(i-1);
		}
   }

   // U Vector
	double Uplus1[x_spaces+1][3];
   double  U[x_spaces+3][3];
   double  F[x_spaces+3][3];

   for (int i=0;i<=x_spaces+2;i++){
         U[i][0] = row0;
         U[i][1] = row0*uend;
         U[i][2] = e0; 
         F[i][0] = row0*uend;
         F[i][1] = row0*pow(uend,2) + Pend;
         F[i][2] = (e0+Pend)*uend;

   }
for(int i=0;i<=x_spaces+2;i++){
	areas[i] = area(x_value[i]);
}
// Putting in the real Mach Equation
// ---------------------------------
	double real_Mach[x_spaces+3];
	double A_star = pow((pow(S0,2)/func(Mstart)),0.5);

	real_Mach[0] = Mstart;
   for (int i=1; i<=x_spaces+2; i++){
		if (i <= x_spaces/2.0+1){
			real_Mach[i] = bisection_search(0, 1.0, areas[i]/A_star);
		}else{
			real_Mach[i] = pbisection_search(1.0, 10, areas[i]/A_star);
		}
	}
//------------------
//end of added script
// loop and work

	double volumes[x_spaces+3];
	for(int i=0;i<=x_spaces+3;i++){
		volumes[i] = delta_x/2.0*(area(x_value[i] + delta_x/2.0) + area(x_value[i] - delta_x/2.0));
	}

//  // Q Vector Initialization
   double  Q[x_spaces+3][3];
   for (int i=0;i<=x_spaces+2;i++){
         Q[i][0] = 0;
			Q[i][1] = Pend*(area(x_value[i] + delta_x/2.0) - area(x_value[i] - delta_x/2.0)); 
         Q[i][2] = 0;
   }

//	//Setting up if time ==0
	for (int i=0; i<=x_spaces+2; i++){
		row[i]= U[i][0];
		velocity[i] = U[i][1]/row[i];
		energy[i] = U[i][2];
		epsilon[i] = energy[i]/row[i] - pow(velocity[i],2)/2.0;
		pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
	}
for (int i=0; i<=x_spaces; i++){
	temp_Mach[i] = velocity[i] / pow(fluid_gamma*pressure[i]/row[i],0.5);
}
//iteration_list[iteration] = iteration;
//residual_list[iteration] = max_residual(temp_Mach, real_Mach,x_spaces+1);


	while(time < endtime){	
	// Finite Volume analysis
		for (int i=1; i<=x_spaces+1; i++){
		   for (int j=0; j<=2; j++){
				FplusHalf = 0.5*(F[i][j] + F[i+1][j])-0.5*alpha*(U[i+1][j]-U[i][j]);
				FminusHalf = 0.5*(F[i-1][j] + F[i][j])-0.5*alpha*(U[i][j]-U[i-1][j]);
				Uplus1[i][j] = U[i][j] - delta_t/volumes[i]*(FplusHalf*area(x_value[i]+delta_x/2.0)- FminusHalf*area(x_value[i] -delta_x/2.0))+delta_t/volumes[i]*Q[i][j];
			}
		}

		for (int i=0; i<=x_spaces+2; i++){
			if (i==x_spaces+2){
				row[i]= row[i-1];
				velocity[i] = velocity[i-1];
				pressure[i] = Pend;
				energy[i] = energy[i-1];
				epsilon[i] = epsilon[i-1];
			}else if (i==x_spaces+1){
				row[i]= Uplus1[i][0];
				velocity[i] = Uplus1[i][1]/row[i];
				pressure[i] = Pend;
				energy[i] = (pressure[i]+row[i]*pow(velocity[i],2)/2.0)/(fluid_gamma-1.0);
				epsilon[i] = energy[i]/row[i] - pow(velocity[i],2)/2.0;
			}else if (i==1){
				row[i]= row0;
				velocity[i] = Uplus1[i][1]/row0;
				energy[i] = Uplus1[i][2];
				epsilon[i] = energy[i]/row0 - pow(velocity[i],2)/2.0;
				pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
				velocity[i] = Mstart*pow(fluid_gamma*pressure[i]/row[i],0.5);
				epsilon[i] = energy[i]/row0 - pow(velocity[i],2)/2.0;
				pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
			}else if (i==0){
				row[i]= row0;
				velocity[i] = Uplus1[i+1][1]/row0;
				energy[i] = Uplus1[i+1][2];
				epsilon[i] = energy[i]/row0 - pow(velocity[i],2)/2.0;
				pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
				velocity[i] = Mstart*pow(fluid_gamma*pressure[i]/row[i],0.5);
				epsilon[i] = energy[i]/row0 - pow(velocity[i],2)/2.0;
				pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
			}else{
				row[i]= Uplus1[i][0];
				velocity[i] = Uplus1[i][1]/row[i];
				energy[i] = Uplus1[i][2];
				epsilon[i] = energy[i]/row[i] - pow(velocity[i],2)/2.0;
				pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
			}
		}

		for(int i=0;i<=x_spaces+2;i++){
			for(int j=0; j<=2;j++){
				if (i == x_spaces+2){
					U[i][j] = U[i-1][j];
				}else if (i==0){
					U[i][j] = Uplus1[i+1][j];
				}else{
					U[i][j] = Uplus1[i][j];	
				}
			}
			if (i == x_spaces+2){
				F[i][0] = F[i-1][0];
				F[i][1] = F[i-1][1];
				F[i][2] = F[i-1][2];
				Q[i][0] = 0;
				Q[i][1] = Q[i-1][1]; 
				Q[i][2] = 0;
			}else if (i==0){
				F[i][0] = row[i+1]*velocity[i+1];
				F[i][1] = row[i+1]*pow(velocity[i+1],2) + pressure[i+1];
				F[i][2] = (energy[i+1]+pressure[i+1])*velocity[i+1];
				Q[i][0] = 0;
				Q[i][1] = pressure[i+1]*(area(x_value[i+1] + delta_x/2.0) - area(x_value[i+1] - delta_x/2.0)); 
				Q[i][2] = 0;
			}else{
				F[i][0] = row[i]*velocity[i];
				F[i][1] = row[i]*pow(velocity[i],2) + pressure[i];
				F[i][2] = (energy[i]+pressure[i])*velocity[i];
				Q[i][0] = 0;
				Q[i][1] = pressure[i]*(area(x_value[i] + delta_x/2.0) - area(x_value[i] - delta_x/2.0)); 
				Q[i][2] = 0;

			}

		}
	//added code
	for (int k=0; k<=x_spaces+2; k++){
		temp_Mach[k] = velocity[k] / pow(fluid_gamma*pressure[k]/row[k],0.5);
	}
	//iteration_list[iteration] = iteration;
	//residual_list[iteration] = max_residual( temp_Mach,real_Mach, x_spaces+3);
	//iteration += 1;
	//end of added code
		time += delta_t;
	}


for (int i=0; i<=x_spaces+2; i++){
	Mach[i] = velocity[i] / pow(fluid_gamma*pressure[i]/row[i],0.5);
}
 //cout << "double uend = pbisection_search(0.001, 1, 3): " << bisection_search(0.001, 1, 3) <<endl;

 printarray (x_value,x_spaces+3, "X Value");
 printarray (row,x_spaces+3, "Y Value");
 printarray (areas,x_spaces+3, "Area");
 printarray (velocity,x_spaces+3, "Velocity");
 printarray (pressure,x_spaces+3, "Pressure");
 printarray (energy,x_spaces+3, "Energy");
 printarray (volumes,x_spaces+3, "Volumes");
 printarray (Mach,x_spaces+3, "Mach"); 
 printarray (velocity,x_spaces+3, "Velocity"); 
 //printarray (real_Mach,x_spaces+1, "Real Mach");
 //printarray (temp_Mach,x_spaces+1, "temp Mach");
 //printarray (residual_list,max_iterations, "residual_list");
 //printarray (iteration_list,max_iterations, "iteration_list");




	return 0; 


}



