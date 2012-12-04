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


   double row[x_spaces+1];
   double velocity[x_spaces+1];
   double energy[x_spaces+1];
   double epsilon[x_spaces+1];
	double areas[x_spaces+1];
   double e0;
   e0 = P0/(fluid_gamma-1) + row0*pow(u0,2)/2.0;

	double knownAreas[x_spaces+1];

	row[0] = row0;
	velocity[0] = u0;
	energy[0] = e0;
	areas[0] = area(0);
	energy[0] = e0;

	long double i;
	int j=0;
   char *inname = "Pressure.txt";
   char *pointer_Areas = "knownAreas.txt";

   ifstream infile(inname);

	double pressure[x_spaces+1];
   if (!infile) {
	   cout << "There was a problem opening file "
           << inname
           << " for reading."
            << endl;
      return 0;
   }
   //cout << "Opened " << inname << " for reading." << endl;
   while (infile >> i) {
		pressure[j] = i;
	   j +=1;
   }
   ifstream areafile(pointer_Areas);
	j=0;
   while (areafile >> i) {
		knownAreas[j] = i;
	   j +=1;
   }

   double x_value[x_spaces+1];
   for (int i=0; i<=x_spaces; i++){
      x_value[i] = x_lower + delta_x*i;
   }

	// the math required to move backwards
	double rowEpsilon;

	for (int i=1; i<=x_spaces; i++){
		rowEpsilon = pressure[i]/(fluid_gamma-1.0);
		energy[i] = rowEpsilon +0.5*(row[i-1]*pow(velocity[i-1],2));
		velocity[i] = (energy[i-1]+pressure[i-1])*velocity[i-1]/(energy[i]+pressure[i]);
		row[i] = (row[i-1]*pow(velocity[i-1],2))/(pow(velocity[i],2));
		areas[i] = (row[i-1]*velocity[i-1]*areas[i-1])/(row[i]*velocity[i]);

	}

	printarray (x_value,x_spaces+1, "X Value");
	printarray (pressure,x_spaces+1, "Pressure Value");
 printarray (velocity,x_spaces+1, "Velocity");
 printarray (pressure,x_spaces+1, "Pressure");
 printarray (energy,x_spaces+1, "Energy");
	printarray (areas,x_spaces+1, "Area");
	printarray (knownAreas,x_spaces+1, "Known Areas");
	return 0; 


}



