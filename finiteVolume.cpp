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
	ofstream myfile;
	ofstream pressureFile;
	myfile.open("finiteVolume.txt");
	pressureFile.open("Pressure.txt");

   double endtime = atof(argv[1]);
   myfile << "The END time is: "<<endtime << endl;
   myfile << "i|Density|Velocity|E "<<endtime << endl;
   double time = 0;

   double row[x_spaces+1];
   double velocity[x_spaces+1];
   double energy[x_spaces+1];
   double epsilon[x_spaces+1];
   double pressure[x_spaces+1];
   double e0;

   e0 = P0/(fluid_gamma-1) + row0*pow(u0,2)/2.0;

   // Creating x value
   double x_value[x_spaces+1];
   for (int i=0; i<=x_spaces; i++){
      x_value[i] = x_lower + delta_x*i;
   }

   // U Vector
	double Uplus1[x_spaces+1][3];
   double  U[x_spaces+1][3];

   for (int i=0;i<=x_spaces;i++){
         U[i][0] = row0;
         U[i][1] = row0*u0;
         U[i][2] = e0; 
   }


   // F Vector
   double  F[x_spaces+1][3];
   for (int i=0;i<=x_spaces;i++){
         F[i][0] = row0*u0;
         F[i][1] = row0*pow(u0,2) + P0;
         F[i][2] = (e0+P0)*u0;
   }

// loop and work
	double areas[x_spaces+1];
	for(int i=0;i<=x_spaces;i++){
		areas[i] = area(x_value[i]);
	}

	double volumes[x_spaces+1];
	for(int i=0;i<=x_spaces;i++){
		if (i == 0){
			volumes[i] = delta_x/2.0*(area(x_value[i]));
		}else{
			volumes[i] = delta_x/2.0*(area(x_value[i] + delta_x/2) + area(x_value[i] - delta_x/2.0));
		}
	}


  // Q Vector Initialization
   double  Q[x_spaces+1][3];
   for (int i=0;i<=x_spaces;i++){
         Q[i][0] = 0;
			if (i==0){
				Q[i][1] = P0*(area(x_value[i])); 
			}else{
		   	Q[i][1] = P0*(area(x_value[i] + delta_x/2) - area(x_value[i] - delta_x/2.0)); 
			}

         Q[i][2] = 0;
   }
	// May need to add Q[0][1] to deal with discontinuity

	//Setting up if time ==0
	for (int i=0; i<=x_spaces; i++){
		row[i]= U[i][0];
		velocity[i] = U[i][1]/row[i];
		energy[i] = U[i][2];
		epsilon[i] = energy[i]/row[i] - pow(velocity[i],2)/2.0;
		pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
	}



	while(time < endtime){	
		myfile << "The time is: "<<time << endl;
		myfile << "i|Density|Velocity|E|P|F(0)|F(1)|F(2)|+S.5|-S.5|Volume"<<endtime << endl;

	// Finite Volume analysis
		for (int i=1; i<=x_spaces; i++){
		   for (int j=0; j<=2; j++){
		      Uplus1[i][j] = U[i][j] - delta_t/volumes[i]*(F[i][j]*area(x_value[i]+delta_x/2.0)- F[i-1][j]*area(x_value[i] -delta_x/2.0))+delta_t/volumes[i]*Q[i][j];     
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
			Q[i][1] = pressure[i]*(area(x_value[i] + delta_x/2) - area(x_value[i] - delta_x/2.0)); 
			Q[i][2] = 0;

		}
		time += delta_t;
	}


double Fplus1[x_spaces+1];
double Fplus2[x_spaces+1];
double Fplus3[x_spaces+1];
for (int i=0; i<=x_spaces;i++){
	
	Fplus1[i] = (F[i][0]);
	Fplus2[i] = (F[i][1]);
	Fplus3[i] = (F[i][2]);
}



 printarray (x_value,x_spaces+1, "X Value");
 printarray (row,x_spaces+1, "Y Value");
 printarray (areas,x_spaces+1, "Area");
 printarray (velocity,x_spaces+1, "Velocity");
 printarray (pressure,x_spaces+1, "Pressure");
 printarray (energy,x_spaces+1, "Energy");
 printarray (volumes,x_spaces+1, "Volumes");

	cout << "\n";
 printarray (Fplus1,x_spaces+1, "Fplus1");
	cout << "\n";
 printarray (Fplus2,x_spaces+1, "Fplus2");
	cout << "\n";
 printarray (Fplus3,x_spaces+1, "Fplus3");
	cout << "\n";

/*
int i=0;
int j=1;
cout << "Term: " << (F[i][j]*areas[i] - F[i-1][j]*areas[i-1]) <<endl;
cout << "Area0: " << areas[0] <<endl;
cout << "Area1: " << areas[1] <<endl;
cout << "F[0][0]: " << F[2][1] <<endl;
cout << "F[1][0]: " << F[3][1] <<endl;
*/

 //printMatrix(U,x_spaces+1,3);
	for (int i=0; i<=x_spaces; i++){
		pressureFile << pressure[i] << endl;
	}

	myfile.close();
	pressureFile.close();
	return 0; 


}



