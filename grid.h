/* Here the grid function gives the area vectors of the all faces, volume of
all the cells   and the barycenter distance(characterstic_length) for all the
cells */

#include "iostream"
#include "math.h"
#include <vector>

const double PI=3.14159265;
using namespace std ; 

// Function defines the area vector and cell volumes  
void grid( vector<vector<vector<vector<float> > > > & x_face_area,
	  vector<vector<vector<vector<float> > > > & y_face_area,
      vector<vector<vector<vector<float> > > > & z_face_area,
      vector<vector<vector<float> > > & V_cell, 
      int Nx, int Ny, int Nz, int N,
      vector<vector<vector<float> > > & characterstic_length)
{
// Grids for bump(this is the first case which will test the scheme)

	// Creating a 4D vector object for grid points
	typedef vector<float> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> matrix4D;

	// this store previous values of variables (density , three momentum, energy)
	matrix4D grid_point(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 

 	float lenght = 10 ; // domain length in meters
 	
 	int Nsb =  2 + 2 ; // stating point of the bump
 	int Neb =  Nsb + N/2 ; // stating point of the bump
 	float bumphight = 0.075 * lenght ; 

	float delta = lenght / N ; 
	float deltatx = delta ;
	float deltaty = lenght/(Ny) ;
	float deltatz = delta ;


	// this file is opend to store the grid poits
	ofstream grid ;
	grid.open("grid.csv");
	grid <<  "x" << "," << "y" <<   endl ;
	
// Grid points without any rotation
// NOTE : In each direction no. of grid points would be just one more 
// than the number of (actual cells + ghost cells)   

// First defining the grid points
	for (int i =0; i < Nx+1; ++i) 
	{
		for (int  j=0;  j < Ny+1; j++)
		{
			for (int  k=0;  k < Nz+1; k++)
			{	
				grid_point[i][j][k][0] = i*deltatx ;   
				grid_point[i][j][k][2] = k*deltatz ;
				
				if( i >= Nsb + 1 && i <= Neb + 1 )
				{
					float del_y_i = bumphight * (i-Nsb)/(Neb-Nsb) ; 
					grid_point[i][j][k][1] = del_y_i + j * (lenght-del_y_i)/(Ny); 
					// this dy should start from dy = delta and than should start decreasing 
				}
				else if(i>= Neb + 2 && i <= Neb + 2 + (Neb - Nsb) )
				{
					grid_point[i][j][k][1] = grid_point[2*Neb+2-i][j][k][1] ;
				}
				else
				{
					grid_point[i][j][k][1] = j*deltaty ;
				}

			}
			// grid << grid_point[i][j][4][0] << "," << grid_point[i][j][4][1] << endl ; 
		}	
	}

	// here comes the area vectors
	for (int i = 0; i  < Nx; ++i)
	{
		for (int  j = 0;  j < Ny; ++j)
		{
			for (int  k = 0;  k < Nz; ++k)
			{
				x_face_area[i][j][k][0] = (grid_point[i][j+1][k][1]-grid_point[i][j][k][1])*deltatz ;
				x_face_area[i][j][k][1] = 0 ;
				x_face_area[i][j][k][2] = 0 ;

				y_face_area[i][j][k][0] = -deltatz*(grid_point[i+1][j][k][1]-grid_point[i][j][k][1]) ;
				y_face_area[i][j][k][1] =  deltatz*(grid_point[i+1][j][k][0]-grid_point[i][j][k][0]) ;
				y_face_area[i][j][k][2] = 0 ;

				z_face_area[i][j][k][0] = 0 ; 
				z_face_area[i][j][k][1] = 0 ;
				z_face_area[i][j][k][2] = 0.5*deltatx*( (grid_point[i][j+1][k][1]-grid_point[i][j][k][1]) + 
					(grid_point[i+1][j+1][k][1]-grid_point[i+1][j][k][1]) ); 
			}
		}	
	}



	// cell volume 
	for (int i = 0; i  < Nx; ++i)
	{
		for (int  j= 0;  j < Ny; ++j)
		{
			for (int  k= 0;  k < Nz; ++k)
			{
				if (i < N+2)
				{
					V_cell[i][j][k] = deltatx*deltaty*deltatz ; 
				}
				else if( i >= N+2 && i <= floor(3*N/2) + 2 )
				{
					V_cell[i][j][k] = 0.5 * deltatx * ((grid_point[i][j+1][k][1] - grid_point[i][j][k][1])  + 
					(grid_point[i+1][j+1][k][1]-grid_point[i+1][j][k][1])) * deltatz ; 
				}
				else 
				{
					V_cell[i][j][k] = V_cell[3*N+4-i][j][k]  ;
				}
			}
		}	
	}
// here are the grid points and area vector with theta degree rotation 
// and the cell volume will not change after the rotation so no neet to worry abuout that

	float theta = PI * 45 / 180 ; // theta is radtion	about the y axis
	for (int i = 0; i <  Nx ; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			for (int k = 0; k < Nz; ++k)
			{
				grid_point[i][j][k][0] = cos(theta)*grid_point[i][j][k][0] + sin(theta)*grid_point[i][j][k][2] ;  
				grid_point[i][j][k][1] = grid_point[i][j][k][1] ; 
				grid_point[i][j][k][2] = -sin(theta)*grid_point[i][j][k][0] + cos(theta)*grid_point[i][j][k][2] ;

				x_face_area[i][j][k][0] = cos(theta)*x_face_area[i][j][k][0] + sin(theta)*x_face_area[i][j][k][2] ;
				x_face_area[i][j][k][1] = x_face_area[i][j][k][1] ;
				x_face_area[i][j][k][2] = -sin(theta)*x_face_area[i][j][k][0] + cos(theta)*x_face_area[i][j][k][2] ;	 
				
				y_face_area[i][j][k][0] = cos(theta)*y_face_area[i][j][k][0] + sin(theta)*y_face_area[i][j][k][2] ;
				y_face_area[i][j][k][1] = y_face_area[i][j][k][1] ;
				y_face_area[i][j][k][2] = -sin(theta)*y_face_area[i][j][k][0] + cos(theta)*y_face_area[i][j][k][2] ;

				z_face_area[i][j][k][0] = cos(theta)*z_face_area[i][j][k][0] + sin(theta)*z_face_area[i][j][k][2] ;
				z_face_area[i][j][k][1] = z_face_area[i][j][k][1] ;
				z_face_area[i][j][k][2] = -sin(theta)*z_face_area[i][j][k][0] + cos(theta)*z_face_area[i][j][k][2] ;
			}
		}
	}

	for (int i = 0; i < Nx; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			grid << grid_point[i][j][4][0] << "," << grid_point[i][j][4][1] << endl ; 
		}
	}

}






























// Flate plate grid

// {
// 	// Grids for straight duct(just to check wheter the basic code working)
// 	float deltax = 0.00012 ; 
// 	float deltay = 0.000045 ;
// 	float deltaz = 0.00014 ;
// 	// this daltax and other conditions will make Re = 1.0e5 

// 	// int nx = 10 + 4; 
// 	// int ny = 20 + 4; 
// 	// int nz = 1 + 4; 

// 	for (int i = 0; i  < nx+1; ++i)
// 	{
// 		for (int  j = 0;  j < ny+1; ++j)
// 		{
// 			for (int  k = 0;  k < nz+1; ++k)
// 			{
// 				// area vector for x direction cells
// 				x_face_area[i][j][k][0] = deltay*deltaz*pow(1.1,j) ;   
// 				x_face_area[i][j][k][1] = 0 ;
// 				x_face_area[i][j][k][2] = 0 ;

// 				// area vector for y direction cells
// 				y_face_area[i][j][k][0] = 0 ;   
// 				y_face_area[i][j][k][1] = deltaz*deltax ;
// 				y_face_area[i][j][k][2] = 0 ;

// 				// area vector for z direction cells
// 				z_face_area[i][j][k][0] = 0 ;   
// 				z_face_area[i][j][k][1] = 0 ;
// 				z_face_area[i][j][k][2] = deltay*deltax*pow(1.1,j) ;
// 			}
// 		}	
// 	}

// 	// cell volume 
// 	for (int i = 0; i  < nx; ++i)
// 	{
// 		for (int  j= 0;  j < ny; ++j)
// 		{
// 			for (int  k= 0;  k < nz; ++k)
// 			{
// 				V_cell[i][j][k] = deltax*deltay*deltaz*pow(1.1,j) ;

// 				characterstic_length[i][j][k] = 10000000 ; 
// 				if (characterstic_length[i][j][k] > deltax)
// 				{
// 					characterstic_length[i][j][k] = deltax ;
// 				}
// 				if (characterstic_length[i][j][k] > deltay*pow(1.1,j))
// 				{
// 					characterstic_length[i][j][k] = deltay*pow(1.1,j) ;
// 				}
// 				if (characterstic_length[i][j][k] > deltaz)
// 				{
// 					characterstic_length[i][j][k] = deltaz ;
// 				}
					
// 			}
// 		}	
// 		// cout << "characterstic_length hai ye " <<
// 		// characterstic_length[i][j][k] << endl ;
// 	}

// // this file is opend to store the grid poits
// 	ofstream grid ;
// 	grid.open("Flat_plate.dat");
// 	// grid <<  "x" << "," << "y" <<   endl ;
// 	for (int i = 2; i < nx-2 ;  ++i)
// 	{
// 		float y = deltay*(1+1.1+(1.1*1.1)/2) ;
// 		for (int j = 2; j < ny-2; ++j)
// 		{
// 			grid << (i+0.5)*deltax << "," << y << endl ; 
// 			y = y + deltay*( pow(1.1,j) + pow(1.1,j+1) )/2 ;
// 		}
// 	}	
// }