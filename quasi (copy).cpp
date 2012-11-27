//File for the 1D Quasi Code
//Jeremy Schembri on Nov 8, 2012
#include <iostream>
#include <math.h>
#include <string>
#include <cstdlib> 
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
   

   // Creating x value
   double x_value[x_spaces+1];
   for (int i=0; i<=x_spaces; i++){
      x_value[i] = x_lower + delta_x*i;
   }

   // Creating Volume
   double volume[x_spaces+1];
   volume[0] = delta_x*area(x_value[0]);
   volume[x_spaces] = delta_x*area(x_value[x_spaces]);

   for (int i=1; i<=x_spaces-1;i++){
      volume[i] = delta_x*(area(x_value[i] + delta_x/2.0)+area(x_value[i]-delta_x/2.0))/2.0;
   }

  double Uplus1[x_spaces+1][3];
   // U Vector
   double  U[x_spaces+1][3];
	U[0][0] = row0;
	U[0][1] = row0*u0;
	U[0][2]= e0;
	Uplus1[0][0] = row0;
	Uplus1[0][1] = row0*u0;
	Uplus1[0][2]= e0;
   
   for (int i=1;i<=x_spaces;i++){
         U[i][0] = rowinf;
         U[i][1] = rowinf*uinf;
         U[i][2] = einf; 
   }


   // F Vector
  double  F[x_spaces+1][3];
        F[0][0] = row0*u0;
        F[0][1] = row0*pow(u0,2) + P0;
        F[0][2]= (e0+P0)*u0;

   for (int i=1;i<=x_spaces;i++){
         F[i][0] = rowinf*uinf;
         F[i][1] = rowinf*pow(uinf,2) + Pinf;
         F[i][2] = (einf+Pinf)*uinf;
   }

  // Q Vector Initialization
  double  Q[x_spaces+1][3];
        Q[0][0] = 0;
        Q[0][1] = P0/area(x_value[0])*der_area(x_value[0]);
        Q[0][2]= 0;

   for (int i=1;i<=x_spaces;i++){
         Q[i][0] = 0;
         Q[i][1] = Pinf/area(x_value[i])*der_area(x_value[i]);;
         Q[i][2] = 0;
   }

while(time <= endtime){	
// Finite Volume analysis
   for (int i=1; i<=x_spaces; i++){
      for (int j=0; j<=2; j++){
         Uplus1[i][j] = U[i][j] - delta_t/volume[i]*(F[i][j]*area(x_value[i]+delta_x/2.0) - F[i-1][j]*area(x_value[i]-delta_x/2.0)) + delta_t/volume[i]*Q[i][j];     
      }
   }

// Unpacking the converservation vectors into density, velocity, energy and pressure
/*
double row[x_spaces+1];
double velocity[x_spaces+1];
double energy[x_spaces+1];
double epsilon[x_spaces+1];
double pressure[x_spaces+1];
*/


for (int i=0; i<=x_spaces; i++){
   row[i]= Uplus1[i][0];
   velocity[i] = Uplus1[i][1]/row[i];
   energy[i] = Uplus1[i][2];
   epsilon[i] = energy[i]/row[i] - pow(velocity[i],2)/2.0;
   pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
}
   //row[0]= U[0][0];
// Resetting all variables for next cycle

for(int i=1;i<=x_spaces;i++){
   for(int j=0; j<=2;j++){
      U[i][j] = Uplus1[i][j];
   }

   
      F[i][0] = row[i]*velocity[i];
      F[i][1] = row[i]*pow(velocity[i],2) + pressure[i];
      F[i][2] = (energy[i]+pressure[i])*velocity[i];

      Q[i][0] = 0;
      Q[i][1] = pressure[i]/area(x_value[i])*der_area(x_value[i]);
      Q[i][2] = 0;

}

   time += delta_t;
}


// loop and work
 double areas[x_spaces+1];
 for(int i=0;i<=x_spaces;i++){
 	areas[i] = area(x_value[i]);
 }


 printarray (x_value,x_spaces+1, "X Value");
 printarray (row,x_spaces+1, "Y Value");
 printarray (areas,x_spaces+1, "Area");
 printarray (velocity,x_spaces+1, "Velocity");
 return 0; 


}



