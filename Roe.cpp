//File for the 1D Quasi Code
//Jeremy Schembri on Nov 8, 2012
#include <iostream>
#include <math.h>
#include <cmath>
#include <string>
#include <cstdlib> 
#include <fstream>
#include "constants.h"
#include <Eigen/Dense>
using namespace Eigen;

using namespace std;

double max_residual( double guess[],double answer[], int length){
	double max;
	for(int i=0;i<length;i++){
		if (i==0){
			max = abs(answer[i]-guess[i] );
		}

		if (abs(answer[i]-guess[i]) > max){
			max = abs(answer[i]-guess[i] );
		}
		//cout << "Max : " << abs(answer[i]-guess[i]) << endl;
	}
	return max;
}

double max_array( double guess[],int length){
	double max;
	for(int i=0;i<length;i++){
		if (i==0){
			max = abs(guess[i] );
		}

		if (abs(guess[i]) > max){
			max = abs(guess[i] );
		}
	}
	return max;
}

MatrixXd A_flux_matrix(double row, double row_plus1, double velocity, double velocity_plus1, double energy, double energy_plus1, double pressure, double pressure_plus1){

	double row_half = pow(row,0.5)*pow(row_plus1,0.5);
	double velocity_half  = (pow(row,0.5)*velocity+pow(row_plus1,0.5)*velocity_plus1) / (pow(row,0.5)+pow(row_plus1,0.5));
	double enthalpy_half = (pow(row,0.5)*(energy+pressure)/row + pow(row_plus1,0.5)*(energy_plus1+pressure_plus1)/row_plus1) / (pow(row,0.5)+pow(row_plus1,0.5));
	double c_at_i = pow((fluid_gamma-1.0)*(enthalpy_half -0.5*pow(velocity_half,2)),0.5);

	//cout << "Row: " << row << endl;
	//cout << "Half Row: " << row_half << endl;
	//cout << "Velocity: " << velocity << endl;
	//cout << "Velocity Half: " << velocity_half << endl;
	//cout << "c_at_i: " << c_at_i << endl;


	MatrixXd S_matrix(3,3);
	MatrixXd S_inverse_matrix(3,3);
	MatrixXd CA_matrix(3,3);
	MatrixXd CA_inverse_matrix(3,3);
	MatrixXd lambda_matrix(3,3);
	MatrixXd lambda_matrix_positive(3,3);
	MatrixXd lambda_matrix_negative(3,3);

	double Alpha_for_matrix = 0.5*pow(velocity_half,2);
	double Beta = fluid_gamma-1.0;

	for (int m=0;m<=2;m++){
		for (int n=0;n<=2;n++){
			S_matrix(m,n) = 0;
			S_inverse_matrix(m,n) = 0;
			CA_matrix(m,n) = 0;
			CA_inverse_matrix(m,n) = 0;
			lambda_matrix(m,n) = 0;
			lambda_matrix_positive(m,n) = 0;
			lambda_matrix_negative(m,n) = 0;
		}
	}



	S_matrix(0,0) = 1; 
	S_matrix(1,0) = -velocity_half / row_half; S_matrix(1,1) = 1.0/row_half;
	S_matrix(2,0) = Alpha_for_matrix*Beta; S_matrix(2,1) = -velocity_half*Beta;S_matrix(2,2) = Beta;
	S_inverse_matrix(0,0) = 1; 
	S_inverse_matrix(1,0) = velocity_half; S_inverse_matrix(1,1) = row_half;
	S_inverse_matrix(2,0) = Alpha_for_matrix; S_inverse_matrix(2,1) = row_half*velocity_half;S_inverse_matrix(2,2) = 1.0/Beta;


	CA_matrix(0,0) = 1; CA_matrix(0,2) = -1.0/pow(c_at_i,2);
	CA_matrix(1,1) = row_half*c_at_i;CA_matrix(1,2) = 1;
	CA_matrix(2,1) = -row_half*c_at_i;CA_matrix(2,2) = 1;
	CA_inverse_matrix(0,0) = 1; CA_inverse_matrix(0,1) = 0.5/pow(c_at_i,2); CA_inverse_matrix(0,2) = 0.5/pow(c_at_i,2);
										 CA_inverse_matrix(1,1) = 0.5/(row_half*c_at_i); CA_inverse_matrix(1,2) = -0.5/(row_half*c_at_i);
										 CA_inverse_matrix(2,1) = 0.5; CA_inverse_matrix(2,2) = 0.5;

	lambda_matrix(0,0) = velocity_half;
	lambda_matrix(1,1) = velocity_half+c_at_i;
	lambda_matrix(2,2) = velocity_half-c_at_i;

	lambda_matrix_positive(0,0) = velocity_half;
	if (velocity_half+c_at_i >=0){
		lambda_matrix_positive(1,1) = velocity_half+c_at_i;
	}
	if (velocity_half-c_at_i >=0){
		lambda_matrix_positive(2,2) = velocity_half-c_at_i;
	}

	lambda_matrix_negative(0,0) = velocity_half;
	if (velocity_half+c_at_i <=0){
		lambda_matrix_negative(1,1) = velocity_half+c_at_i;
	}
	if (velocity_half-c_at_i <=0){
		lambda_matrix_negative(2,2) = velocity_half-c_at_i;
	}

	MatrixXd A(3,3);
	MatrixXd Aplus(3,3);
	MatrixXd Aminus(3,3);
	MatrixXd fluxA(3,3);
	//cout << "S_inverse_matrix :\n" << S_inverse_matrix <<endl;
	//cout << "SCA_inverse_matrix :\n" << CA_inverse_matrix <<endl;
	//cout << "lambda_matrix_positi :\n" << lambda_matrix_positive <<endl;
	//cout << "CA_matrix :\n" << CA_matrix <<endl;
	//cout << "S_matrix :\n" << S_matrix<<endl;

	Aplus = S_inverse_matrix*CA_inverse_matrix*lambda_matrix_positive*CA_matrix*S_matrix;
	A = S_inverse_matrix*CA_inverse_matrix*lambda_matrix*CA_matrix*S_matrix;
	Aminus = A - Aplus;
	fluxA = Aplus - Aminus;
	//cout << "Aplus: \n" << Aplus << endl;
	//cout << "Aminus: \n" << Aminus << endl;
	//cout << "A: \n" << fluxA << endl;
	return fluxA;
}




int main(int argc, char **argv){

	//Setting up values from 0 to i+1
   double endtime = atof(argv[1]);
   double time = 0;
   double row[x_spaces+1];
   double velocity[x_spaces+1];
   double energy[x_spaces+1];
   double epsilon[x_spaces+1];
   double pressure[x_spaces+1];
   double Mach[x_spaces+1];
   double temperature[x_spaces+1];
   double areas[x_spaces+1];
   double alpha;// = 500;
   double FplusHalf,FminusHalf;

	//Adding info regarding Roe Scheme
	MatrixXd AplusHalf(3,3);
	MatrixXd AminusHalf(3,3);


	//ending info for Roe Scheme
	
	//solving for u0
	double denominator = 1.0/(Mach0*Mach0*fluid_gamma*R*Tt) - (fluid_gamma-1.0)/(fluid_gamma+1.0)*1.0/(aStar*aStar); 
	double u0 = pow(1.0/denominator,0.5);
	double P0 = static_pressure(u0);
	//double P0 = Pt;
	double T0 = static_temperature(u0);
	//double T0 = Tt;
	double row0 = P0 / (R*T0);
   double e0 = P0/(fluid_gamma-1.0) + row0*pow(u0,2)/2.0;
	


   double temp_Mach[x_spaces+1];
	int count = 0;	
	int iteration = 0;
	int skip = 100;
	int max_count = (endtime/delta_t)/skip;
	double iteration_list[max_count+1];
	double residual_list[max_count+1];



	double temp1[x_spaces];
	double temp2[x_spaces];
	double temp3[x_spaces];

	double residual_R1[max_count+1];
	double residual_R2[max_count+1];
	double residual_R3[max_count+1];

   // Creating x value
   double x_value[x_spaces+1];
   for (int i=0; i<=x_spaces; i++){
      x_value[i] = x_lower + delta_x*i;
   }

   // U Vector
	double Uplus1[x_spaces+1][3];
   double  U[x_spaces+1][3];
   double  F[x_spaces+1][3];

   for (int i=0;i<=x_spaces;i++){
         U[i][0] = row0;
         U[i][1] = row0*u0;
         U[i][2] = e0; 
         F[i][0] = row0*u0;
			if (i==x_spaces){
		      F[i][1] = row0*pow(u0,2) + Pend;
		      F[i][2] = (e0+Pend)*u0;
			}else{
		      F[i][1] = row0*pow(u0,2) + P0;
		      F[i][2] = (e0+P0)*u0;
			}

   }
for(int i=0;i<=x_spaces;i++){
	areas[i] = area(x_value[i]);
}
// Putting in the real Mach Equation
// ---------------------------------
   double M_chock_start,c_chocke_start;
	double real_Mach[x_spaces+1];
	c_chocke_start = pow(fluid_gamma*P0/row0,0.5); 
	M_chock_start= u0 /c_chocke_start;;
	double A_star = pow((pow(S0,2)/func(M_chock_start)),0.5);


	real_Mach[0] = M_chock_start;

   for (int i=1; i<=x_spaces; i++){
		if (i <= x_spaces/2.0){
			real_Mach[i] = pbisection_search(1.1, 10, areas[i]/A_star);
		}else{
			real_Mach[i] = pbisection_search(1.1, 10, areas[i]/A_star);
		}
	}
//------------------
//end of added script
// loop and work

	double volumes[x_spaces+1];
	for(int i=0;i<=x_spaces;i++){
			volumes[i] = delta_x/2.0*(area(x_value[i] + delta_x/2.0) + area(x_value[i] - delta_x/2.0));
	}

//  // Q Vector Initialization
   double  Q[x_spaces+2][3];
   for (int i=0;i<=x_spaces;i++){
         Q[i][0] = 0;
			if (i==x_spaces){
				Q[i][1] = Pend*(area(x_value[i] + delta_x/2.0) - area(x_value[i] - delta_x/2.0)); 
			}else{
				Q[i][1] = P0*(area(x_value[i] + delta_x/2.0) - area(x_value[i] - delta_x/2.0)); 
			}
         Q[i][2] = 0;
   }

//	//Setting up if time ==0
	for (int i=0; i<=x_spaces; i++){
		row[i]= U[i][0];
		velocity[i] = U[i][1]/row[i];
		energy[i] = U[i][2];
		epsilon[i] = energy[i]/row[i] - pow(velocity[i],2)/2.0;
		if (i==x_spaces){
			pressure[i] = Pend;
		}else{
			pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
		}
		temperature[i] = T0;
	}
for (int i=0; i<=x_spaces; i++){
	temp_Mach[i] = velocity[i] / pow(fluid_gamma*pressure[i]/row[i],0.5);
}
iteration_list[count] = iteration;
residual_list[count] = max_residual(temp_Mach, real_Mach,x_spaces+1);

	VectorXf deltaU(3);

	while(time < endtime){	
	// Finite Volume analysis
		for (int i=1; i<=x_spaces-1; i++){
			//added code for Roe scheme included for i+0.5

			//ended code for Roe Scheme


		   for (int j=0; j<=2; j++){
//				if (i <=9){
//					alpha = 10*(i);
//					//alpha = 200;
//				}else if (i<=x_spaces-1 && i>=x_spaces-10){
//					alpha = 10*((x_spaces-1)-i+1);
//					//alpha = 200;
//				}else{
//					alpha = 200;
//				}
//				FplusHalf = 0.5*(F[i][j] + F[i+1][j])-0.5*alpha*(U[i+1][j]-U[i][j]);
//				FminusHalf = 0.5*(F[i-1][j] + F[i][j])-0.5*alpha*(U[i][j]-U[i-1][j]);

				AplusHalf = A_flux_matrix(row[i], row[i+1], velocity[i], velocity[i+1], energy[i], energy[i+1], pressure[i], pressure[i+1]);
				AminusHalf = A_flux_matrix(row[i-1], row[i], velocity[i-1], velocity[i], energy[i-1], energy[i], pressure[i-1], pressure[i]);

				deltaU(0) = U[i+1][0] - U[i][0];
				deltaU(1) = U[i+1][1] - U[i][1];
				deltaU(2) = U[i+1][2] - U[i][2];
				FplusHalf = 0.5*(F[i][j] + F[i+1][j]) -0.5*(AplusHalf(j,0)*deltaU(0) +AplusHalf(j,1)*deltaU(1)+AplusHalf(j,2)*deltaU(2));  					
			
				deltaU(0) = U[i][0] - U[i-1][0];
				deltaU(1) = U[i][1] - U[i-1][1];
				deltaU(2) = U[i][2] - U[i-1][2];
				FminusHalf = 0.5*(F[i-1][j] + F[i][j]) -0.5*(AminusHalf(j,0)*deltaU(0) +AminusHalf(j,1)*deltaU(1)+AminusHalf(j,2)*deltaU(2)); 

				Uplus1[i][j] = U[i][j] - delta_t/volumes[i]*(FplusHalf*area(x_value[i]+delta_x/2.0)- FminusHalf*area(x_value[i] -delta_x/2.0))+delta_t/volumes[i]*Q[i][j];

			/*
			if (i ==50){
				cout << "\n Aplus Half: " << AplusHalf <<endl;
				cout << "1,2" << AplusHalf(1,2) <<endl;
				cout << "Aminus Half "<< AminusHalf << endl; 
				cout << "Alpha1: " << (AplusHalf(j,0)*deltaU(0) +AplusHalf(j,1)*deltaU(1)+AplusHalf(j,2)*deltaU(2)) << endl;
				cout << "Alpha2: " << (AminusHalf(j,0)*deltaU(0) +AminusHalf(j,1)*deltaU(1)+AminusHalf(j,2)*deltaU(2)) << endl;
			}
			*/

			}
		}

		for (int i=1; i<=x_spaces-1; i++){
				row[i]= Uplus1[i][0];
				velocity[i] = Uplus1[i][1]/row[i];
				energy[i] = Uplus1[i][2];
				epsilon[i] = energy[i]/row[i] - pow(velocity[i],2)/2.0;
				pressure[i] = (fluid_gamma -1.0)*row[i]*epsilon[i];
				temperature[i] = pressure[i] / (R*row[i]);
		}

		//special condition when i ==0, the first cell (only for supersonic)
		temperature[0] = T0;
		double body = 1.0 -(fluid_gamma-1.0)/(fluid_gamma+1.0)*pow(velocity[0]/aStar,2);
		double exponent = 1.0/(fluid_gamma-1.0); 
		double deltaP_over_delta_u = Pt*(fluid_gamma/(fluid_gamma-1.0))*pow(body,exponent)*(-2.0*(fluid_gamma-1.0)/(fluid_gamma+1.0)*velocity[0]/pow(aStar,2));
		double CFL = velocity[0]*delta_t/delta_x;
		double c0 = pow(fluid_gamma*pressure[0]/row[0],0.5);
		double c1 = pow(fluid_gamma*pressure[1]/row[1],0.5);
		double delta_t1 = CFL*delta_x /(velocity[0] + c0);
		double lambda = delta_t1/delta_x*(0.5*(velocity[1]+velocity[0]) - 0.5*(c1+c0));
		double delta_u = (-lambda*(pressure[1] - pressure[0] - row[0]*c0*(velocity[1] - velocity[0])))/(deltaP_over_delta_u-row[0]*c0);
		velocity[0] = velocity[0] + delta_u;
		temperature[0] = static_temperature(velocity[0]);
		//temperature[0] = Tt;
		pressure[0] = Pt*pow(temperature[0]/Tt,fluid_gamma/(fluid_gamma-1.0));
		//pressure[0] = Pt;
		row[0] = pressure[0] / (R*temperature[0]);
		energy[0] = row[0]*(c_v*temperature[0] +0.5*pow(velocity[0],2));
		c0 = pow(fluid_gamma*pressure[0]/row[0],0.5);
		Mach[0] = velocity[0] / c0;
		//added code by Jeremy
//		Mach[0] = Mach0;
//		velocity[0] = c0 * Mach0;


		//exit boundary conditions
		CFL = velocity[x_spaces]*delta_t/delta_x;
		double c_xspaces =pow(fluid_gamma*pressure[x_spaces]/row[x_spaces],0.5); 
		double c_maxMinusOne =pow(fluid_gamma*pressure[x_spaces-1]/row[x_spaces-1],0.5);  
		double delta_tMAX = CFL*delta_x/(velocity[x_spaces] + c_xspaces);
		double lambda1 = 0.5*(delta_tMAX/delta_x)*(velocity[x_spaces] + velocity[x_spaces-1]);
		double lambda2 =  0.5*(delta_tMAX/delta_x)*(velocity[x_spaces] + velocity[x_spaces-1]+c_xspaces+c_maxMinusOne);
		double lambda3 =  0.5*(delta_tMAX/delta_x)*(velocity[x_spaces] + velocity[x_spaces-1]-c_xspaces-c_maxMinusOne);

		double R1_lambda1 = -lambda1*(row[x_spaces] - row[x_spaces-1] - 1.0/pow(c_xspaces,2)*(pressure[x_spaces] - pressure[x_spaces-1]));
		double R2_lambda2 = -lambda2*(pressure[x_spaces]-pressure[x_spaces-1] + row[x_spaces]*c_xspaces*(velocity[x_spaces] - velocity[x_spaces-1]));
		double R3_lambda3 = -lambda3*(pressure[x_spaces]-pressure[x_spaces-1] - row[x_spaces]*c_xspaces*(velocity[x_spaces] - velocity[x_spaces-1]));
		double delta_P, delta_row;


		Mach[x_spaces] = (velocity[x_spaces] + velocity[x_spaces-1])/(c_xspaces+c_maxMinusOne);
		if (Mach[x_spaces] > 1){
			delta_P = 0.5*(R2_lambda2+R3_lambda3);
		}else{
			delta_P = 0;
		}

		delta_row = R1_lambda1 + delta_P/pow(c_xspaces,2);
		delta_u = (R2_lambda2 -delta_P)/(row[x_spaces]*c_xspaces);

		//updating properties
		row[x_spaces] = row[x_spaces] +delta_row;
		velocity[x_spaces] = velocity[x_spaces] +delta_u;
		pressure[x_spaces] = pressure[x_spaces] + delta_P;
		temperature[x_spaces] =  pressure[x_spaces]/(row[x_spaces]*R);
		energy[x_spaces] = row[x_spaces]*(c_v*temperature[x_spaces]+0.5*pow(velocity[x_spaces],2));
		c_xspaces =pow(fluid_gamma*pressure[x_spaces]/row[x_spaces],0.5); 
		Mach[x_spaces] = velocity[x_spaces] / c_xspaces;

		//end of exit boundary conditions


		for(int i=0;i<=x_spaces;i++){
			for(int j=0; j<=2;j++){
				U[i][j] = Uplus1[i][j];	

				F[i][0] = row[i]*velocity[i];
				F[i][1] = row[i]*pow(velocity[i],2) + pressure[i];
				F[i][2] = (energy[i]+pressure[i])*velocity[i];
				Q[i][0] = 0;
				Q[i][1] = pressure[i]*(area(x_value[i] + delta_x/2.0) - area(x_value[i] - delta_x/2.0)); 
				Q[i][2] = 0;
			}
		}

	//added code
	for (int k=0; k<=x_spaces; k++){
		temp_Mach[k] = velocity[k] / pow(fluid_gamma*pressure[k]/row[k],0.5);
	}
	if (iteration % skip ==0){
		iteration_list[count] = iteration;
		residual_list[count] = max_residual( temp_Mach,real_Mach, x_spaces+1);
		residual_R1[count] = max_array(temp1,x_spaces); 
		residual_R2[count] = max_array(temp2,x_spaces); 
		residual_R3[count] = max_array(temp3,x_spaces); 
		count += 1;
	}

	double zeroes[x_spaces];
	for (int i=1;i<=x_spaces-1;i++){
		for (int j=0;j<=2;j++){
				/*

				if (i <=9){
					alpha = 10*(i);
					//alpha = 200;
				}else if (i<=x_spaces-1 && i>=x_spaces-10){
					alpha = 10*((x_spaces-1)-i+1);
					//alpha = 200;
				}else{
					alpha = 200;
				}
			*/
				AplusHalf = A_flux_matrix(row[i], row[i+1], velocity[i], velocity[i+1], energy[i], energy[i+1], pressure[i], pressure[i+1]);
				AminusHalf = A_flux_matrix(row[i-1], row[i], velocity[i-1], velocity[i], energy[i-1], energy[i], pressure[i-1], pressure[i]);

				deltaU(0) = U[i+1][0] - U[i][0];
				deltaU(1) = U[i+1][1] - U[i][1];
				deltaU(2) = U[i+1][2] - U[i][2];
				FplusHalf = 0.5*(F[i][j] + F[i+1][j]) -0.5*(AplusHalf(j,0)*deltaU(0) +AplusHalf(j,1)*deltaU(1)+AplusHalf(j,2)*deltaU(2));  					
			
				deltaU(0) = U[i][0] - U[i-1][0];
				deltaU(1) = U[i][1] - U[i-1][1];
				deltaU(2) = U[i][2] - U[i-1][2];
				FminusHalf = 0.5*(F[i-1][j] + F[i][j]) -0.5*(AminusHalf(j,0)*deltaU(0) +AminusHalf(j,1)*deltaU(1)+AminusHalf(j,2)*deltaU(2)); 
		

			//FplusHalf = 0.5*(F[i][j] + F[i+1][j])-0.5*alpha*(U[i+1][j]-U[i][j]);
			//FminusHalf = 0.5*(F[i-1][j] + F[i][j])-0.5*alpha*(U[i][j]-U[i-1][j]);
			if (j==0){
				temp1[i-1] = (delta_t/volumes[i]*abs(FplusHalf*area(x_value[i]+delta_x/2.0)- FminusHalf*area(x_value[i] -delta_x/2.0))); 
			}else if(j==1){
 				temp2[i-1]  =(delta_t/volumes[i]*abs(FplusHalf*area(x_value[i]+delta_x/2.0)- FminusHalf*area(x_value[i] -delta_x/2.0)-Q[i][j]));
			}else if (j==2){
 				temp3[i-1] = (delta_t/volumes[i]*abs(FplusHalf*area(x_value[i]+delta_x/2.0)- FminusHalf*area(x_value[i] -delta_x/2.0)));
			}
		}
		zeroes[i-1] = 0;
	}

	//residual_R1[iteration/skip] = max_array(temp1,x_spaces); 
	//residual_R1[iteration/skip] = max_array(temp1,x_spaces); 
	//residual_R1[iteration/skip] = max_array(temp2,x_spaces); 
	//residual_R2[iteration/skip] = max_array(temp3,x_spaces); 


	iteration += 1;
	
	//end of added code
		time += delta_t;
	}

		for (int i=0; i<=x_spaces; i++){
			Mach[i] = velocity[i] / pow(fluid_gamma*pressure[i]/row[i],0.5);
		}

		for (int i=1; i<=x_spaces-1; i++){
			//cout << "i: temp1: " << i<<": " << temp1[i]<<endl;
		}
		//cout << "Max array: " << max_array(temp1,x_spaces)<<endl;

 printarray (x_value,x_spaces+1, "X Value");
 printarray (row,x_spaces+1, "Y Value");
 printarray (areas,x_spaces+1, "Area");
 printarray (velocity,x_spaces+1, "Velocity");
 printarray (pressure,x_spaces+1, "Pressure");
 printarray (energy,x_spaces+1, "Energy");
 printarray (temperature,x_spaces+1, "Temperature");
 printarray (volumes,x_spaces+1, "Volumes");
 printarray (Mach,x_spaces+1, "Mach"); 
 //printarray (real_Mach,x_spaces+1, "Real Mach");
 //printarray (temp_Mach,x_spaces+1, "temp Mach");
// printarray (residual_list,max_count, "residual_list");
 printarray (iteration_list,max_count, "iteration_list"); 
 printarray (residual_R1,max_count, "R1");
 printarray (residual_R2,max_count, "R2");
 printarray (residual_R3,max_count, "R3");

	//cout << "\n Aplus Half: \n" << AplusHalf <<endl;
	//cout << "Aminus Half \n"<< AminusHalf << endl; 
	//cout << "deltaU \n" << deltaU << endl;




	return 0; 


}



