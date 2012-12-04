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

	//initalize filename

   double time = 0;

   double Mach[x_spaces+1];
   double areas[x_spaces+1];
   double row[x_spaces+1];
   double velocity[x_spaces+1];
   double energy[x_spaces+1];
   double epsilon[x_spaces+1];
   double pressure[x_spaces+1];
   double temperature[x_spaces+1];
   double e0, M0;
   double x_value[x_spaces+1];

   e0 = P0/(fluid_gamma-1) + row0*pow(u0,2)/2.0;
	M0 = u0 / pow(fluid_gamma*(fluid_gamma-1.0)*e0,0.5);
	double A_star = pow((pow(S0,2)/func(M0)),0.5);

   for (int i=0; i<=x_spaces; i++){
      x_value[i] = x_lower + delta_x*i;
		areas[i] = area(x_value[i]);
   }
	
	Mach[0] = M0;
   for (int i=1; i<=x_spaces; i++){
		Mach[i] = bisection_search(0.001, 0.9, areas[i]/A_star);
	}

   // Creating x value
//	cout << "Test case: "<< bisection_search(0.001, 0.9, 2) <<endl;
//	cout << "func(0.31): "<< func(0.31) <<endl;
//	cout << "func(2.2) "<< func(2.2) <<endl;

//	cout << "E0: " << e0 <<endl;
//	cout << "S0: " << S0 <<endl;
//	cout << "A_star: " << A_star <<endl;

	printarray (x_value,x_spaces+1, "X Value");
	printarray (Mach,x_spaces+1, "Mach");

 return 0; 


}



