/* Here the grid function gives the area vectors of the all faces, volume of
all the cells   and the barycenter distance(characterstic_length) for all the
cells */

#include "iostream"
#include "math.h"
#include <vector>

using namespace std ; 

// Function defines the area vector and cell volumes  
void grid( vector<vector<vector<vector<float> > > > & x_face_area,
	  vector<vector<vector<vector<float> > > > & y_face_area,
      vector<vector<vector<vector<float> > > > & z_face_area,
      vector<vector<vector<float> > > & V_cell, 
      int nx, int ny, int nz, 
      vector<vector<vector<float> > > & characterstic_length)
{

	// Grids for straight duct(just to check wheter the basic code working)
	float deltax = 0.00012 ; 
	float deltay = 0.000045 ;
	float deltaz = 0.00014 ;
	// this daltax and other conditions will make Re = 1.0e5 

	// int nx = 10 + 4; 
	// int ny = 20 + 4; 
	// int nz = 1 + 4; 

	for (int i = 0; i  < nx+1; ++i)
	{
		for (int  j = 0;  j < ny+1; ++j)
		{
			for (int  k = 0;  k < nz+1; ++k)
			{
				// area vector for x direction cells
				x_face_area[i][j][k][0] = deltay*deltaz*pow(1.1,j) ;   
				x_face_area[i][j][k][1] = 0 ;
				x_face_area[i][j][k][2] = 0 ;

				// area vector for y direction cells
				y_face_area[i][j][k][0] = 0 ;   
				y_face_area[i][j][k][1] = deltaz*deltax ;
				y_face_area[i][j][k][2] = 0 ;

				// area vector for z direction cells
				z_face_area[i][j][k][0] = 0 ;   
				z_face_area[i][j][k][1] = 0 ;
				z_face_area[i][j][k][2] = deltay*deltax*pow(1.1,j) ;
			}
		}	
	}

	// cell volume 
	for (int i = 0; i  < nx; ++i)
	{
		for (int  j= 0;  j < ny; ++j)
		{
			for (int  k= 0;  k < nz; ++k)
			{
				V_cell[i][j][k] = deltax*deltay*deltaz*pow(1.1,j) ;

				characterstic_length[i][j][k] = 10000000 ; 
				if (characterstic_length[i][j][k] > deltax)
				{
					characterstic_length[i][j][k] = deltax ;
				}
				if (characterstic_length[i][j][k] > deltay*pow(1.1,j))
				{
					characterstic_length[i][j][k] = deltay*pow(1.1,j) ;
				}
				if (characterstic_length[i][j][k] > deltaz)
				{
					characterstic_length[i][j][k] = deltaz ;
				}
					
			}
		}	
		// cout << "characterstic_length hai ye " <<
		// characterstic_length[i][j][k] << endl ;
	}

// this file is opend to store the grid poits
	ofstream kullu_grid ;
	kullu_grid.open("Flat_plate.dat");
	// kullu_grid <<  "x" << "," << "y" <<   endl ;
	for (int i = 2; i < nx-2 ;  ++i)
	{
		float y = deltay*(1+1.1+(1.1*1.1)/2) ;
		for (int j = 2; j < ny-2; ++j)
		{
			kullu_grid << (i+0.5)*deltax << "," << y << endl ; 
			y = y + deltay*( pow(1.1,j) + pow(1.1,j+1) )/2 ;
		}
	}	
}