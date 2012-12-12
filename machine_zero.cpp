//File for the 1D Quasi Code
//Jeremy Schembri on Nov 8, 2012
#include <iostream>
#include <math.h>
#include <cmath>
#include <string>
#include <cstdlib> 
#include <fstream>
#include "constants.h"

using namespace std;


int main(int argc, char **argv){

	double machineZero = 1.0;
	int count = 1;
	while (1+ machineZero != 1.0 && count < 1000){
		machineZero = machineZero/2.0;
	printf("Machine Zero is: %e \n", machineZero);
		cout << count << endl;
	count += 1;
	}

	printf("Machine Zero is: %f \n", machineZero);
	return 0; 


}



