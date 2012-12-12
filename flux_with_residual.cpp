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
		//cout << "Max : " << abs(answer[i]-guess[i]) << endl;
	}
	return max;
}



int main(int argc, char **argv){

	//Setting up values from 0 to i+1
   double endtime = atof(argv[1]);
   double time = 0;
   double row[x_spaces+2];
   double velocity[x_spaces+2];
   double energy[x_spaces+2];
   double epsilon[x_spaces+2];
   double pressure[x_spaces+2];
   double Mach[x_spaces+2];
   double areas[x_spaces+2];
   double alpha = 100;
   double e0, FplusHalf,FminusHalf;





   double temp_Mach[x_spaces+1];
	int iteration = 0;
	int max_iterations = endtime/delta_t;
	double iteration_list[max_iterations+1];
	double residual_list[max_iterations+1];
	//cout << "max_iterations: " <<max_iterations<<endl;
   e0 = P0/(fluid_gamma-1.0) + row0*pow(u0,2)/2.0;

	double temp1[x_spaces];
	double temp2[x_spaces];
	double temp3[x_spaces];

	double residual_R1[max_iterations+1];
	double residual_R2[max_iterations+1];
	double residual_R3[max_iterations+1];

   // Creating x value
   double x_value[x_spaces+2];
   for (int i=0; i<=x_spaces+1; i++){
      x_value[i] = x_lower + delta_x*i;
		if (i==x_spaces+1){
			x_value[i] = x_value[i-1];
		}
   }

   // U Vector
	double Uplus1[x_spaces+1][3];
   double  U[x_spaces+2][3];
   double  F[x_spaces+2][3];

   for (int i=0;i<=x_spaces+1;i++){
         U[i][0] = row0;
         U[i][1] = row0*u0;
         U[i][2] = e0; 
         F[i][0] = row0*u0;
         F[i][1] = row0*pow(u0,2) + P0;
         F[i][2] = (e0+P0)*u0;

   }
for(int i=0;i<=x_spaces+1;i++){
	areas[i] = area(x_value[i]);
}
// Putting in the real Mach Equation
// ---------------------------------
   double M0,c;
	double real_Mach[x_spaces+1];
	c = pow(fluid_gamma*P0/row0,0.5); 
	M0 = u0 /c;
	double A_star = pow((pow(S0,2)/func(M0)),0.5);


	real_Mach[0] = M0;

   for (int i=1; i<=x_spaces; i++){
		if (i <= x_spaces/2.0){
			real_Mach[i] = pbisection_search(1.1, 10, areas[i]/A_star);
		}else{
			real_Mach[i] = pbisection_search(1.1, 10, areas[i]/A_star);
		}
	}
//------------------
//end of added script
// loop and work

	double volumes[x_spaces+1];
	for(int i=0;i<=x_spaces;i++){
			volumes[i] = delta_x/2.0*(area(x_value[i] + delta_x/2.0) + area(x_value[i] - delta_x/2.0));
//			cout << "delta_x/2.0: " << delta_x/2.0<<endl;
//			cout << "area(x_value[i] + delta_x/2.0): " << area(x_value[i] + delta_x/2.0) << endl;
//			cout << "area(x_value[i] - delta_x/2.0): " << area(x_value[i] - delta_x/2.0) << endl;
	}
//	cout << "area(-1): " << area(-1) <<endl;

//  // Q Vector Initialization
   double  Q[x_spaces+2][3];
   for (int i=0;i<=x_spaces+1;i++){
         Q[i][0] = 0;
			Q[i][1] = P0*(area(x_value[i] + delta_x/2.0) - area(x_value[i] - delta_x/2.0)); 
         Q[i][2] = 0;
   }

//	//Setting up if time ==0
	for (int i=0; i<=x_spaces+1; i++){
		row[i]= U[i][0];
		velocity[i] = U[i][1]/row[i];
		energy[i] = U[i][2];
		epsilon[i] = energy[i]/row[i] - pow(velocity[i],2)/2.0;
		pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
	}
for (int i=0; i<=x_spaces; i++){
	temp_Mach[i] = velocity[i] / pow(fluid_gamma*pressure[i]/row[i],0.5);
}
iteration_list[iteration] = iteration;
residual_list[iteration] = max_residual(temp_Mach, real_Mach,x_spaces+1);


	while(time < endtime){	
	// Finite Volume analysis
		for (int i=1; i<=x_spaces; i++){
		   for (int j=0; j<=2; j++){
				FplusHalf = 0.5*(F[i][j] + F[i+1][j])-0.5*alpha*(U[i+1][j]-U[i][j]);
				FminusHalf = 0.5*(F[i-1][j] + F[i][j])-0.5*alpha*(U[i][j]-U[i-1][j]);
				Uplus1[i][j] = U[i][j] - delta_t/volumes[i]*(FplusHalf*area(x_value[i]+delta_x/2.0)- FminusHalf*area(x_value[i] -delta_x/2.0))+delta_t/volumes[i]*Q[i][j];
			}
		}

		for (int i=1; i<=x_spaces+1; i++){
			if (i==x_spaces+1){
				row[i]= row[i-1];
				velocity[i] = velocity[i-1];
				energy[i] = energy[i-1];
				epsilon[i] = epsilon[i-1];
				pressure[i] = pressure[i-1];
			}else{
				row[i]= Uplus1[i][0];
				velocity[i] = Uplus1[i][1]/row[i];
				energy[i] = Uplus1[i][2];
				epsilon[i] = energy[i]/row[i] - pow(velocity[i],2)/2.0;
				pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
			}
		}

		for(int i=1;i<=x_spaces+1;i++){
			for(int j=0; j<=2;j++){
				if (i == x_spaces+1){
					U[i][j] = U[i-1][j];
				}else{
					U[i][j] = Uplus1[i][j];	
				}
			}
			if (i == x_spaces+1){
				F[i][0] = F[i-1][0];
				F[i][1] = F[i-1][1];
				F[i][2] = F[i-1][2];
				Q[i][0] = 0;
				Q[i][1] = Q[i-1][1]; 
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
	for (int k=0; k<=x_spaces; k++){
		temp_Mach[k] = velocity[k] / pow(fluid_gamma*pressure[k]/row[k],0.5);
	}
	iteration_list[iteration] = iteration;
	residual_list[iteration] = max_residual( temp_Mach,real_Mach, x_spaces+1);


	double zeroes[x_spaces];
	for (int i=1;i<=x_spaces;i++){
		for (int j=0;j<=2;j++){
			FplusHalf = 0.5*(F[i][j] + F[i+1][j])-0.5*alpha*(U[i+1][j]-U[i][j]);
			FminusHalf = 0.5*(F[i-1][j] + F[i][j])-0.5*alpha*(U[i][j]-U[i-1][j]);
			if (j==0){
				temp1[i-1] = (delta_t/volumes[i]*abs(FplusHalf*area(x_value[i]+delta_x/2.0)- FminusHalf*area(x_value[i] -delta_x/2.0))); 
			}else if(j==1){
 				temp2[i-1]  = 0.001*(delta_t/volumes[i]*abs(FplusHalf*area(x_value[i]+delta_x/2.0)- FminusHalf*area(x_value[i] -delta_x/2.0)-Q[i][j]));
			}else if (j==2){
 				temp3[i-1] = 0.00001*(delta_t/volumes[i]*abs(FplusHalf*area(x_value[i]+delta_x/2.0)- FminusHalf*area(x_value[i] -delta_x/2.0)));
			}
		}
		zeroes[i-1] = 0;
	}
	residual_R1[iteration] = max_residual(zeroes,temp1,x_spaces); 
	residual_R2[iteration] = max_residual(zeroes,temp2,x_spaces); 
	residual_R3[iteration] = max_residual(zeroes,temp3,x_spaces); 

	iteration += 1;
	//end of added code
		time += delta_t;
	}


for (int i=0; i<=x_spaces; i++){


	Mach[i] = velocity[i] / pow(fluid_gamma*pressure[i]/row[i],0.5);
}





 printarray (x_value,x_spaces+1, "X Value");
 printarray (row,x_spaces+1, "Y Value");
 printarray (areas,x_spaces+1, "Area");
 printarray (velocity,x_spaces+1, "Velocity");
 printarray (pressure,x_spaces+1, "Pressure");
 printarray (energy,x_spaces+1, "Energy");
 printarray (volumes,x_spaces+1, "Volumes");
 printarray (Mach,x_spaces+1, "Mach"); 
 //printarray (real_Mach,x_spaces+1, "Real Mach");
 //printarray (temp_Mach,x_spaces+1, "temp Mach");
 printarray (residual_list,max_iterations, "residual_list");
 printarray (iteration_list,max_iterations, "iteration_list"); 
 printarray (residual_R1,max_iterations, "R1");
 printarray (residual_R2,max_iterations, "R2");
 printarray (residual_R3,max_iterations, "R3");

	//for (int i=0;i<=x_spaces;i++){
		//cout << "i: " << abs(FplusHalf*area(x_value[i]+delta_x/2.0)- FminusHalf*area(x_value[i] -delta_x/2.0)) << endl;
	//}





	return 0; 


}



