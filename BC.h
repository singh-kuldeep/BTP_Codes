#include "iostream"
#include "math.h"
#include <vector>

#define GAMMA 1.4

// Function Implements the boundary condition
void BC(vector<vector<vector<vector<float> > > > & U,
	vector<vector<vector<vector<float> > > > & y_face_area,
	vector<vector<vector<vector<float> > > > & z_face_area,
	 int Nx, int Ny, int Nz)
{
	// inilet conditions (user given data)
	// one has to maintion only inlet Mach number totalpressure and
	// // the totaltemperature
	// float Static_pressure = 1e5 - 0.56 ; // Atmoshperic pressure
	// // float temperature = 273.15 ; // Atmoshperic pressure
	// // float density = 1.2918786 ;
	// // float Mach_free_stream = 0.40; // this is given by user 
	// // float free_Stream_Velocity = 331.368688 ; // M*u
	// // float total_energy = 324240.0 ;
	// float pressure ; 

	// // Inlet ghost cell updating usinig the totel 
	// //quantities(T_0, P_0) and flow direction  
	// for (int j =2; j < Ny-2; ++j)
	// {
	// 	for (int k =2; k < Nz-2; ++k)
	// 	{	
	// 		// these values are given based on the freestream conditions
	// 		pressure = 0.4*(U[3][j][k][4] - 0.5*( pow(U[3][j][k][1],2) + pow(
	// 			U[3][j][k][2],2) +  pow(U[3][j][k][1],2) )/U[3][j][k][0] ) ; 

	// 		U[0][j][k][0] = 1.2918786 ; //density ;
	// 		U[0][j][k][1] = 171.2352464; //inletdensity*inletvelocity ;
	// 		U[0][j][k][2] = 0 ;
	// 		U[0][j][k][3] = 0 ;
	// 		U[0][j][k][4] = (pressure/0.4)+11348.40 ; 
	// 		//inletpressure/(GAMMA-1) + 0.5*inletdensity*inletvelocity*
	// 		//inletvelocity ;

	// 		pressure = 0.4*(U[2][j][k][4] - 0.5*( pow(U[2][j][k][1],2) + pow(
	// 			U[2][j][k][2],2) + pow(U[2][j][k][1],2) )/U[2][j][k][0] ) ; 
			
	// 		U[1][j][k][0] = 1.2918786 ; //density ;
	// 		U[1][j][k][1] = 171.2352464 ;//inletdensity*inletvelocity ;
	// 		U[1][j][k][2] = 0 ;
	// 		U[1][j][k][3] = 0 ;
	// 		U[1][j][k][4] = (pressure/0.4)+11348.40 ; 
	// 		//inletpressure/(GAMMA-1) + 0.5*inletdensity*inletvelocity*
	// 		//inletvelocity ;
	// 	}
	// }




	// /* at exit updating the x ghost cell (Aim is to keep exit static pressure
	//  constant througout the iteration) and 
	//  P_live = P_ghost because the flow is subsonic*/ 

	// for (int j =2; j < Ny-2; ++j)
	// {
	// 	for (int k =2; k < Nz-2; ++k)
	// 	{
	// 		// p_live = p_ghost => p_live  =  Static_pressure 
	// 		U[Nx-3][j][k][4] = ( Static_pressure/(GAMMA-1) +  0.5*( pow(
	// 			U[Nx-3][j][k][1],2) + pow(U[Nx-3][j][k][2],2)+pow(
	// 			U[Nx-3][j][k][3],2) )/U[Nx-3][j][k][0] ) ;

	// 		U[Nx-4][j][k][4] = ( Static_pressure/(GAMMA-1) +  0.5*( pow(
	// 			U[Nx-4][j][k][1],2) + pow(U[Nx-4][j][k][2],2)+pow(
	// 			U[Nx-4][j][k][3],2) )/U[Nx-4][j][k][0] ) ;

	// 		for (int l = 0; l < 4; ++l) 
			//  rho is same and velocity is linearly extraploeted  
			// pressure is little less than the inlet pressure due to viscocity
			// so e will alos be different   
			
	// 		{
	// 			U[Nx-2][j][k][l] = U[Nx-3][j][k][l]  + (U[Nx-3][j][k][l] - U[Nx-4][j][k][l]) ;
	// 			U[Nx-1][j][k][l] = U[Nx-2][j][k][l]  + (U[Nx-2][j][k][l] - U[Nx-3][j][k][l]) ;
	// 		} 
	// 		U[Nx-2][j][k][4] = ( Static_pressure/(GAMMA-1) +  0.5*( pow(
	// 			U[Nx-2][j][k][1],2) + pow(U[Nx-2][j][k][2],2)+pow(
	// 			U[Nx-2][j][k][3],2) )/U[Nx-2][j][k][0] ) ;

	// 		U[Nx-1][j][k][4] = ( Static_pressure/(GAMMA-1) +  0.5*( pow(
	// 			U[Nx-1][j][k][1],2) + pow(U[Nx-1][j][k][2],2)+pow(
	// 			U[Nx-1][j][k][3],2) )/U[Nx-1][j][k][0] ) ;

	// 	}
	// }

// inilet conditions (user given data)
	// one has to maintion only inlet Mach number totalpressure and the totaltemperature
	float totaltemperature = 350 ;
	float totalpressure = 492766.8;
	float totaldensity = totalpressure / (gasconstant*totaltemperature) ;

	// Inlet ghost cell updating usinig the totel quantities(T_0, P_0) and flow direction  
	for (int j =2; j < Ny-2; ++j)
	{
		for (int k =2; k < Nz-2; ++k)
		{
			float inletpressure = (gamma-1)*(U[2][j][k][4] - 0.5*(pow(U[2][j][k][1],2) +
				pow(U[2][j][k][2],2)+pow(U[2][j][k][3],2))/U[2][j][k][0]) ;

			float Mach = sqrt( (2/(gamma-1)) * ( pow( (totalpressure/inletpressure),((gamma-1)/gamma) ) -1 ) ) ;

			float inlettemperature = totaltemperature / (1+((gamma-1)*Mach*Mach)/2);

			float inletvelocity = Mach * sqrt(gamma*gasconstant*inlettemperature) ;

			float inletdensity = totaldensity / pow((totalpressure/inletpressure),(1/gamma)) ;

			U[0][j][k][0] = inletdensity ;
			U[0][j][k][1] = inletdensity*inletvelocity ;
			U[0][j][k][2] = 0 ;
			U[0][j][k][3] = 0 ;
			U[0][j][k][4] = inletpressure/(gamma-1) + 0.5*inletdensity*inletvelocity*inletvelocity ;

			U[1][j][k][0] = inletdensity ;
			U[1][j][k][1] = inletdensity*inletvelocity ;
			U[1][j][k][2] = 0 ;
			U[1][j][k][3] = 0 ;
			U[1][j][k][4] = inletpressure/(gamma-1) + 0.5*inletdensity*inletvelocity*inletvelocity ;
		}
	}

// at exit updating the x ghost cell (this is true whene flow is supersonic)
	for (int j =2; j < Ny-2; ++j)
	{
		for (int k =2; k < Nz-2; ++k)
		{
			for (int l = 0; l < 5; ++l)
			{
				U[Nx-2][j][k][l] = U[Nx-3][j][k][l] ;
				U[Nx-1][j][k][l] = U[Nx-4][j][k][l] ;
			}
		}
	}


	// updating the wall ghost cell(Y - wall) 
	//[Here the bottom wall (y = 0 ) is viscous ]
	for (int i = 2; i < Nx-2; ++i)
	{
		for (int k = 2; k < Nz-2; ++k)
		{
			//using rightcell(j=2) filling j=1
			// using zero order extrapolation
			U[i][1][k][0] =  U[i][2][k][0] ;
			U[i][1][k][1] =  U[i][2][k][1] ;
			U[i][1][k][2] = -U[i][2][k][2] ;
			U[i][1][k][3] =  U[i][2][k][3] ;
			U[i][1][k][4] =  U[i][2][k][4] ;

			//using rightcell(j=3) filling j=0
			// using zero order extrapolation
			U[i][0][k][0] =  U[i][3][k][0] ;
			U[i][0][k][1] =  U[i][3][k][1] ;
			U[i][0][k][2] = -U[i][3][k][2] ;
			U[i][0][k][3] =  U[i][3][k][3] ;
			U[i][0][k][4] =  U[i][3][k][4] ;

			//using rightcell(j=Ny-3) filling j=Ny-2
			// here every term has been multiplied by -1 
			//because cell area vector is -A
			// using zero order extrapolation
			U[i][Ny-2][k][0] =  U[i][Ny-3][k][0] ;
			U[i][Ny-2][k][1] =  U[i][Ny-3][k][1] ;
			U[i][Ny-2][k][2] = -U[i][Ny-3][k][2] ;
			U[i][Ny-2][k][3] =  U[i][Ny-3][k][3] ;
			U[i][Ny-2][k][4] =  U[i][Ny-3][k][4] ;

			// using rightcell(j=Ny-4) filling j=Ny-1
			// here every term has been multiplied by -1 
			//because cell area vector is -A
			// using zero order extrapolation
			U[i][Ny-1][k][0] =  U[i][Ny-4][k][0] ;
			U[i][Ny-1][k][1] =  U[i][Ny-4][k][1] ;
			U[i][Ny-1][k][2] = -U[i][Ny-4][k][2] ;
			U[i][Ny-1][k][3] =  U[i][Ny-4][k][3] ;
			U[i][Ny-1][k][4] =  U[i][Ny-4][k][4] ;

		}
	}

	// // updating the z ghost cells
	for (int i = 2; i < Nx-2; ++i)
	{
		for (int j = 2; j < Ny-2; ++j)
		{
			// using cell (k=2) filling k=1
			// using zero order extrapolation
			U[i][j][1][0] =  U[i][j][2][0] ;
			U[i][j][1][1] =  U[i][j][2][1] ;
			U[i][j][1][2] =  U[i][j][2][2] ;
			U[i][j][1][3] = -U[i][j][2][3] ;
			U[i][j][1][4] =  U[i][j][2][4] ;

			// using cell (k=3) filling k=0
			// using zero order extrapolation
			U[i][j][0][0] =  U[i][j][3][0] ;
			U[i][j][0][1] =  U[i][j][3][1] ;
			U[i][j][0][2] =  U[i][j][3][2] ;
			U[i][j][0][3] = -U[i][j][3][3] ;
			U[i][j][0][4] =  U[i][j][3][4] ;

			// using cell (k=Nz-3) filling k=Nz-2
			// using zero order extrapolation
			U[i][j][Nz-2][0] =  U[i][j][Nz-3][0] ;
			U[i][j][Nz-2][1] =  U[i][j][Nz-3][1] ;
			U[i][j][Nz-2][2] =  U[i][j][Nz-3][2] ;
			U[i][j][Nz-2][3] = -U[i][j][Nz-3][3] ;
			U[i][j][Nz-2][4] =  U[i][j][Nz-3][4] ;

			//using cell(k=Nz-4) filling k=Nz-1
			// using zero order extrapolation
			U[i][j][Nz-1][0] =  U[i][j][Nz-4][0] ;
			U[i][j][Nz-1][1] =  U[i][j][Nz-4][1] ;
			U[i][j][Nz-1][2] =  U[i][j][Nz-4][2] ;
			U[i][j][Nz-1][3] = -U[i][j][Nz-4][3] ;
			U[i][j][Nz-1][4] =  U[i][j][Nz-4][4] ;
		}
	}
}

