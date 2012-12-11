#ifndef CONSTANTS_H
#define CONSTANTS_H

using namespace std;

extern double delta_x ;
extern double delta_t;

extern double row0;
extern double u0;
extern double e0;
extern double P0;
extern double S0;
extern double Mstar;

extern double rowinf;
extern double uinf;
extern double einf;
extern double Pinf;
extern double Pend;
extern double fluid_gamma;


extern double x_lower;
extern double x_higher;
extern int x_spaces;
extern double PI;

extern int delta(float x, float a);

// Equations
double area(double x);
double der_area(double x);
double func(double M);
double binary_func(double M, double A_ratio);
double bisection_search(double umin, double umax, double A_ratio);
double pbisection_search(double umin, double umax, double A_ratio);
void printarray (double arg[], int length, string input);
//void printMatrix(double U[],int length,int depth);

#endif

