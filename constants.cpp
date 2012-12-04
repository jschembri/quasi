// Settings the constants for the 1D Quasi Code
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string>
#include <cstdlib> 

using namespace std;

double row0 = 1.2;
double u0 = 100;
double P0 = 100000;
double fluid_gamma =1.4;




double x_lower = 0;
double x_higher = 10;
int x_spaces = 10;
double delta_x = (x_higher - x_lower)/x_spaces;
double delta_t = 0.0001; //(in seconds)
double PI = 3.141592654;
double c_v = 718;

int delta(float x, float a){
   if (x >= a){
      return 1;
   }
   else {
      return 0;
   }   

}


//double area(double x){
//   return 10*((-1.0/10.0*pow(x,2)+5)*(delta(x,-1)-delta(x,4)) +(0.4*pow((x-5),2)+3)*(delta(x,4)-delta(x,6)) +  (-1.0/10.0*pow((x-10),2)+5)*(delta(x,6))); 
//}

/*
double area(double x){
   return 10*((-1.0/10.0*pow(x,2)+5)*(delta(x,0)-delta(x,4)) +(0.4*pow((x-5),2)+3)*(delta(x,4)-delta(x,6)) +  (-1.0/10.0*pow((x-10),2)+5)*(delta(x,6))); 
}
*/

//// the derivative of  the area
//double der_area(double x){
//   return (-2.0/10.0*x)*(delta(x,0)-delta(x,4)) +(0.8*(x-5))*(delta(x,4)-delta(x,6)) +  (-2.0/10.0*(x-10))*(delta(x,6));
//}



//The area profile
double area(double x){
	double R = 10 + 1/(x_higher-x_lower)*x;
   return PI*pow(R,2); 
}

//// the derivative of  the area
//double der_area(double x){
//   double A = 2.0/(x_higher-x_lower);
//   return 2*PI*(1 + 2*A*x)*(2*A);

//}


//double area(double x){
//	double R = pow(0.1*(x-5.0),2)+1;
//   return R; 
//}

/*
double area(double x){
	double R = pow(0.1*(x-5.0),2)+1;
   return R; 
}
*/




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

double func(double M){
	double place_holder = 2.0/(fluid_gamma+1.0)*(1+((fluid_gamma-1.0)/2.0*pow(M,2)));
	double exponent = (fluid_gamma+1.0)/(fluid_gamma-1.0);
	return 1.0/pow(M,2)*pow(place_holder,exponent);
}

double binary_func(double M, double A_ratio){
	return func(M) - pow(A_ratio,2);
}

double bisection_search(double umin, double umax, double A_ratio){
	double uguess, ymin, ymax, yguess;
	while (umax-umin>0.0001){
		uguess = (umin + umax)/2;

		ymin = -binary_func(umin,A_ratio);
		yguess = -binary_func(uguess,A_ratio);
		ymax = -binary_func(umax,A_ratio);

		if (yguess > 0){
			umax = uguess;
		} else if (yguess < 0){
			umin = uguess;
		} else if (yguess == 0){
			return uguess;

		}
	//		cout <<"ymin : " << ymin  << endl;
	//		cout <<"yGuess: " << yguess << endl;
	//		cout <<"ymax : " << ymax  << endl;
	}

	return uguess;
}




double S0 = area(x_lower);



