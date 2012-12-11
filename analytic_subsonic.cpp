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

   double Mach[x_spaces+1];
   double areas[x_spaces+1];
   double x_value[x_spaces+1];
	double row[x_spaces+1];
	double pressure[x_spaces+1];

	double A_star = pow((pow(S0,2)/func(Mstart)),0.5);

   for (int i=0; i<=x_spaces; i++){
      x_value[i] = x_lower + delta_x*i;
		areas[i] = area(x_value[i]);
   }
	
	Mach[0] = Mstart;
   for (int i=1; i<=x_spaces; i++){
		if (i <= x_spaces/2.0){
			Mach[i] = bisection_search(0, 1.0, areas[i]/A_star);
		}else{
			Mach[i] = pbisection_search(1.0, 10, areas[i]/A_star);
		}
	}

	double A, exponent;
	A = (1.0 + 0.5*(fluid_gamma-1.0)*pow(Mach[x_spaces],2));
	exponent = 1.0/(fluid_gamma-1.0);
	row[0] = row0;
	double P0_subsonic = Pend*pow(A,fluid_gamma*exponent);
	//pressure[0] = P0;
	for (int i=0;i<=x_spaces;i++){
		A = (1.0 + 0.5*(fluid_gamma-1.0)*pow(Mach[i],2));
		row[i] = row0 / pow(A,exponent);
		pressure[i] = P0_subsonic*pow(A,-1.0*fluid_gamma*exponent);
	}


	//cout << "P0_subsonic" << P0_subsonic <<endl;
	//cout << "Mach[x_spaces]" << Mach[x_spaces]<<endl;

	printarray (x_value,x_spaces+1, "X Value");
	printarray (areas,x_spaces+1, "Areas");
	printarray (Mach,x_spaces+1, "Mach");
	printarray (row,x_spaces+1, "Density");
	printarray (pressure,x_spaces+1, "Pressure");
 return 0; 


}



