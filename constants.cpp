// Settings the constants for the 1D Quasi Code
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string>
#include <cstdlib> 

using namespace std;

double row0 = 1.2;
double u0 = 100;
//double e0 = 214000;
double P0 = 1000;

//double rowinf = 1.2;
//double uinf = 100;
//double einf = 214000;
//double Pinf = 1000;
double fluid_gamma =1.4;




double x_lower = 0;
double x_higher = 10;
int x_spaces = 100;
double delta_x = (x_higher - x_lower)/x_spaces;
double delta_t = 0.0001; //(in seconds)
double PI = 3.141592654;

int delta(float x, float a){
   if (x >= a){
      return 1;
   }
   else {
      return 0;
   }   

}

<<<<<<< HEAD

double area(double x){
   return 10*((-1.0/10.0*pow(x,2)+5)*(delta(x,0)-delta(x,4)) +(0.4*pow((x-5),2)+3)*(delta(x,4)-delta(x,6)) +  (-1.0/10.0*pow((x-10),2)+5)*(delta(x,6))); 
}
=======
/*
double area(double x){
   return 10*((-1.0/10.0*pow(x,2)+5)*(delta(x,0)-delta(x,4)) +(0.4*pow((x-5),2)+3)*(delta(x,4)-delta(x,6)) +  (-1.0/10.0*pow((x-10),2)+5)*(delta(x,6))); 
}
*/
>>>>>>> 6a7cfb924453230b4b24d5ec6a1d5ae1ef939766

//// the derivative of  the area
//double der_area(double x){
//   return (-2.0/10.0*x)*(delta(x,0)-delta(x,4)) +(0.8*(x-5))*(delta(x,4)-delta(x,6)) +  (-2.0/10.0*(x-10))*(delta(x,6));
//}



//The area profile
double area(double x){
	double R = 10 + 0.1/(x_higher-x_lower)*x;
   return PI*pow(R,2); 
}

//// the derivative of  the area
//double der_area(double x){
//   double A = 2.0/(x_higher-x_lower);
//   return 2*PI*(1 + 2*A*x)*(2*A);

//}

<<<<<<< HEAD
//double area(double x){
//	double R = pow(0.1*(x-5.0),2)+1;
//   return R; 
//}
=======
/*
double area(double x){
	double R = pow(0.1*(x-5.0),2)+1;
   return R; 
}
*/


>>>>>>> 6a7cfb924453230b4b24d5ec6a1d5ae1ef939766

// the derivative of  the area
//double der_area(double x){
//   double A = 2.0/(x_higher-x_lower);
//   return (0.2*x-1);

//}



void printarray (double arg[], int length, string input) {
  for (int i=0; i<length; i++)
      if (i==0){
         cout << input <<" Start" << "," << arg[i] << ",";
      }else if (i==length-1){
         cout << arg[i] << "," << input <<" End" << ",";
      }else{
         cout << arg[i] << ",";
      }
}

/*
void printMatrix(double U[],int length,int depth){
	for (int i=0;i<length;i++){
		for (int j=0; j<depth; j++){
			if (i != length) {
				cout << U[i][j] <<" ,";
			}else{
				cout << U[i][j] << endl;
			}
		}
	}
}
*/

long double S0 = area(x_lower);



