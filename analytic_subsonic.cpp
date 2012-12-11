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


	double A_star = pow((pow(S0,2)/func(Mstar)),0.5);

   for (int i=0; i<=x_spaces; i++){
      x_value[i] = x_lower + delta_x*i;
		areas[i] = area(x_value[i]);
   }
	
	Mach[0] = Mstar;
   for (int i=1; i<=x_spaces; i++){
		if (i <= x_spaces/2.0){
			Mach[i] = bisection_search(0, 1.0, areas[i]/A_star);
		}else{
			Mach[i] = pbisection_search(1.0, 10, areas[i]/A_star);
		}
	}


	printarray (x_value,x_spaces+1, "X Value");
	printarray (areas,x_spaces+1, "Areas");
	printarray (Mach,x_spaces+1, "Mach");
 return 0; 


}



