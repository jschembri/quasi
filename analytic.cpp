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
   double e0, M0,c;
   double x_value[x_spaces+1];

   e0 = P0/(fluid_gamma-1.0) + row0*pow(u0,2)/2.0;
	c = pow(fluid_gamma*P0/row0,0.5); 
	M0 = u0 /c;
	double A_star = pow((pow(S0,2)/func(M0)),0.5);

   for (int i=0; i<=x_spaces; i++){
      x_value[i] = x_lower + delta_x*i;
		areas[i] = area(x_value[i]);
   }
	
	Mach[0] = M0;
   for (int i=1; i<=x_spaces; i++){
		if (i <= x_spaces/2.0){
			Mach[i] = pbisection_search(1.1, 10, areas[i]/A_star);
		}else{
			Mach[i] = pbisection_search(1.1, 10, areas[i]/A_star);
		}
	}
	double A, exponent;
   for (int i=0; i<=x_spaces; i++){
      row[i] = row0*pow(1+0.5*(fluid_gamma-1.0)*pow(Mach[i],2),(fluid_gamma-1.0));
		A = 1+0.5*(fluid_gamma-1.0)*pow(Mach[i],2);
		exponent = -1.0*fluid_gamma / (fluid_gamma-1.0);
		pressure[i] = P0*pow(A,exponent);
   }
   // Creating x value
//	cout << "Test case: "<< bisection_search(0.001, 0.9, 2) <<endl;
//	cout << "func(0.31): "<< func(0.31) <<endl;
//	cout << "func(2.2) "<< func(2.2) <<endl;


// 	cout << "func(0.377) "<< func(0.377) <<endl;
//	cout << "binary_func(0.377, 100.0/60.0): "<<binary_func(0.377, 100/60.0)<<endl;
// 	cout << "area(5): " << area(5) <<endl;
//	cout << "u_crit: "<<bisection_search(0.001, 0.9, 100/60.0)<<endl;
//	cout << "S0: " << S0 <<endl;
//	cout << "A_star: " << A_star <<endl;

	printarray (x_value,x_spaces+1, "X Value");
	printarray (areas,x_spaces+1, "Areas");
	printarray (Mach,x_spaces+1, "Mach");
//	printarray (row,x_spaces+1, "Row");
//	printarray (pressure,x_spaces+1, "Pressure");

 return 0; 


}



