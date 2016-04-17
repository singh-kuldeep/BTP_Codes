/*1. Lower the exit pressure ?(the pressure can be calculated using the analytical formula) 
2. implement the exit BC using the dp/dx at n-1 and n cell
3. charls harsh
*/
/*! 
@mainpage Main Solver 
* -Few points regarding the solver: 
*   1. One has to give the CFL number and the total number of time steps 
*	2. Detla_t is locally chosen for all the cells by keeping CFL constant
*   3. There is separate class which calculates the delta_t 
*   4. Steady flow 
*  More text here.
*/

/*
subscript _l is for "Left"
subscript _l_m is for "Left_minus"
subscript _r is for "Right"
subscript _r_p is for "Right_plus"
subscript _A is for "Area"
subscript _V is for "Volume"
*/
#include "iostream"
#include <vector>
#include <fstream>
#include "math.h"
#include "time.h"
#include "netflux.h"
#include "BC.h"
#include "grid.h"
#include "delta_t.h" 

#define GAMMA 1.4
#define R 287.14
#define CP 717.5

using namespace std ;

// BC function implements the boundary condition
void BC(vector<vector<vector<vector<float> > > > & U, 
	vector<vector<vector<vector<float> > > > & y_face_area, 
	vector<vector<vector<vector<float> > > > & z_face_area, 
	int nx, int ny, int nz) ;

/*
grid function generates the area vector and volume for the all cells in the
domain the function definition is in grid.h header file and the working is 
explained in the same file
*/
void grid(vector<vector<vector<vector<float> > > > & x_face_area, 
		  vector<vector<vector<vector<float> > > > & y_face_area, 
		  vector<vector<vector<vector<float> > > > & z_face_area, 
		  vector<vector<vector<float> > > & V_cell, 
		  int nx, int ny, int nz, 
		  vector<vector<vector<float> > > & characterstic_length) ;

int main ()
{
	time_t start, end ; /*These are the start and end point of the timer*/
	time(&start); /*Here timer starts counting the time*/ 

	int nt = 100000; /* Total time steps which will be input 
	given from the user*/
	float CFL = 0.1 ; /* CFL number give and stability condition and it 
	should be less than 1 always and the delta_t is directly depends on 
	CFL number*/

	// float delta_x = 0.0011 ; 
	// float delta_y = 0.00015 ;
	// float delta_z = 0.001 ;
	
	int nx = (10 + 4);
	int ny = 20 + 4 ;
	int nz = 1 + 4 ; 
	//Only one live cell in z direction because of 2D flow
	/* These are the grid specification 
		1. Structured grid has been used 
		2. If user wants denser grid than he can increase the nx, ny and nz
		3. Please remember that nx, ny and nz also includes the ghost cell so 
		actual cell no will be 4 less in count, in each direction 
	*/
	
	typedef vector<float> Dim1;
	typedef vector<Dim1>  Dim2;
	typedef vector<Dim2>  Dim3;
	typedef vector<Dim3>  matrix4D;
	// Creating a 1D, 2D, 3D and 4D vector object

	// this store previous time step values of conserved variables
	//conserved variables (rho, Momentum, energy)
	matrix4D U(nx,Dim3(ny,Dim2(nz,Dim1(5)))); 

	// this store new time step values of conserved variables 
	matrix4D U_new(nx,Dim3(ny,Dim2(nz,Dim1(5)))); 
	
	matrix4D x_face_area(nx+1,Dim3(ny+1,Dim2(nz+1,Dim1(3)))); 
	matrix4D y_face_area(nx+1,Dim3(ny+1,Dim2(nz+1,Dim1(3)))); 
	matrix4D z_face_area(nx+1,Dim3(ny+1,Dim2(nz+1,Dim1(3)))); 

	Dim3 V_cell(nx,Dim2(ny,Dim1(nz)));
	Dim3 characterstic_length(nx,Dim2(ny,Dim1(nz))); // Characteristic length
	//These variables stores grid information(Area vectors and cell volume)
	/*x_face_area, y_face_area, z_face_area store the area vectors of "Right
	 faces" of each cell in the 3D*/  
	
	grid(x_face_area,y_face_area,z_face_area,V_cell,nx,ny,nz,
		characterstic_length);

	/* By calling the grid function grids information is generated and grid 
	points are stored in a file */

	for (int i = 0; i < nx; ++i)
	{
		for (int j = 0; j < ny; ++j)
		{
			for (int k = 0; k < nz; ++k)
			{
				U[i][j][k][0] = 1.16;
				U[i][j][k][1] = 0;
				U[i][j][k][2] = 0;
				U[i][j][k][3] = 0;
				U[i][j][k][4] = 271767;

				U_new[i][j][k][0] = 1.16;
				U_new[i][j][k][1] = 0;
				U_new[i][j][k][2] = 0;
				U_new[i][j][k][3] = 0;
				U_new[i][j][k][4] = 271767;
			}
		}
	}
	/*Initial condition(these are just random values), just to initialize the 
	primitive variables so at the very beginning(first time step) solver can 
	use these to calculate the next time step(second time step) value of
	primitive variables	*/

	ofstream extra_data_file ;
	extra_data_file.open("Extra_data_file.dat");
	extra_data_file << nx-4 << "," << ny-4 << "," << nz-4 << "," << 1.2918786 
	<<"," << (171.2352464/1.2918786) << endl ; 
	/*Live cells grid points, free-stream density and free-stream velocity 
	is being stored in the Extra_data_file.dat file, values are stored in the
	CSV format. These values are required while plotting the pressure, 
	velocity etc. in the MATLAB */

	ofstream Residual_file ;
	Residual_file.open("Residual.dat");
	//Residual_file <<  "t(secs)" << "," << "rho_residual"  << "," << "
	//x_momentum_residual" << "," << "y_momentum_residual" << "," << "
	//z_momentum_residual" << "," << "energy_residual" << endl ;

	/*Residual.dat file contains all the primitive variable residuals of each
	 time step, in CSV format */

	for (int t = 0; t < nt; ++t)
	// time marching starts here 	
	{

		BC(U,y_face_area,z_face_area,nx,ny,nz) ; 
		/*To apply proper boundary condition, before every time step we need to
		 have proper primitive variable value in the ghost cells, which is 
		 done by calling the BC function. This takes care of all Inlet, Exit,
		 y-wall, Z-wall boundary condition*/

		/*Here three for loop is for go to each and every cell in the 3D 
		domain*/
		for (int i = 1; i < nx-2; ++i)
		{
			for (int j = 1; j < ny-2; ++j)
			{
				for (int k = 1; k < nz-2; ++k)
				{
					float x_interface_volume = 0.5*(V_cell[i][j][k] + 
						V_cell[i+1][j][k]) ;
					//x right interface volume

					float y_interface_volume = 0.5*(V_cell[i][j][k] + 
						V_cell[i][j+1][k]) ;
					//y right interface volume

					float z_interface_volume = 0.5*(V_cell[i][j][k] + 
						V_cell[i][j][k+1]) ;
					//z right interface volume

					delta_t delta_t(U[i][j][k],
						characterstic_length[i][j][k],CFL) ;
					/*Because we are interested in the steady flow solution so
					 there is no need to have a global delta_t, we can vary 
					 delta_t locally by keeping the	CFL number constant 
					 throughout the simulation, that is what class delta_t
					 does, here delta_t may change from point to point that 
					 will depend on the local velocity and the grid size*/

					netflux x_face(U,x_face_area,
						y_face_area,z_face_area,V_cell,delta_t.dt,i-1,
						j,k,i,j,k,i+1,j,k,i+2,j,k);

					netflux y_face(U,x_face_area,
						y_face_area,z_face_area,V_cell,delta_t.dt,i,j-1
						,k,i,j,k,i,j+1,k,i,j+2,k);

					netflux z_face(U,x_face_area,
						y_face_area,z_face_area,V_cell,delta_t.dt,i,j,
						k-1,i,j,k,i,j,k+1,i,j,k+2);
					/*For each cell location net flux vector is calculated at
					 all three "Right faces", using class net flux
					 net_flux contains three fluxes Euler, viscus and
					 numerical diffusion flux*/ 
					
					for (int l = 0; l < 5; ++l)
					{
						U_new[i][j][k][l] -=(delta_t.
							dt/x_interface_volume)*(x_face.net_flux[l]);
						U_new[i+1][j][k][l]+=(delta_t
							.dt/x_interface_volume)*(x_face.net_flux[l]);

						U_new[i][j][k][l] -=(delta_t.
							dt/y_interface_volume)*(y_face.net_flux[l]);
						U_new[i][j+1][k][l]+=(delta_t
							.dt/y_interface_volume)*(y_face.net_flux[l]);

						U_new[i][j][k][l] -=(delta_t.
							dt/z_interface_volume)*(z_face.net_flux[l]);
						U_new[i][j][k+1][l]+=(delta_t
							.dt/z_interface_volume)*(z_face.net_flux[l]);

						/*each cell has 3 Right face x_face, y_face and
						z_face so in the above calculation net flux(which is
						going through the right faces of the cell) is 
						subtracted from the main cell and added to it's three
						Right adjusent cells by multiplying the appropriate 
						delta_t and the volume of the interface cell*/
					}
				}
			}
		}

		float rho_residual = 0.0; 
		float x_momentum_residual = 0.0;
		float y_momentum_residual = 0.0;
		float z_momentum_residual = 0.0;
		float energy_residual = 0.0;
		/*Five above maintained parameters are the residual. Residuals help
		analyzing the stability and the convergence. Here the second norm was
		used for the residual calculation*/
		for (int x = 2; x < nx-2; ++x)
			{
				for (int y = 2; y < ny-2; ++y)
				{
					rho_residual   += pow((U_new
					[x][y][2][0]-U[x][y][2][0]),2);

					x_momentum_residual+= pow((U_new
					[x][y][2][1]-U[x][y][2][1]),2);

					y_momentum_residual+= pow((U_new
					[x][y][2][2]-U[x][y][2][2]),2);

					z_momentum_residual+= pow((U_new
					[x][y][2][3]-U[x][y][2][3]),2);

					energy_residual    += pow((U_new
					[x][y][2][4]-U[x][y][2][4]),2);     
				}
			}

		Residual_file << t << "," << sqrt(rho_residual/((nx-4)*(ny-4)*(nz-4))) << "," << 
		sqrt(x_momentum_residual/((nx-4)*(ny-4)*(nz-4))) << "," << sqrt(y_momentum_residual/((nx-4)*(ny-4)*(nz-4)))
		<< "," << sqrt(z_momentum_residual/((nx-4)*(ny-4)*(nz-4))) << "," << 
		sqrt(energy_residual/((nx-4)*(ny-4)*(nz-4))) << endl ;
		//Writing the residuals of each time step in the Residual_file 
		// float dt = delta_t.dt ;
		cout << "timestep   " << t << "   delta t   " << 0 << "   Residual  " <<rho_residual <<endl ; 

		for (int i = 2; i < nx-2; ++i)
		{
			for (int j = 2; j < ny-2; ++j)
			{
				for (int k = 2; k < nz-2; ++k)
				{
					for (int l = 0; l < 5; ++l)
					{
						U[i][j][k][l] = U_new[i][j][k][l] ;
					}				
				}
			}
		}
		/*before going to the new time-step update U
		by U_new*/

		if (t%250 == 0)
		{
			ofstream kullu_2D ;
			kullu_2D.open("conserved_variables.dat");
			//kullu_2D << "rho" << "," << "rho*u" << "," << "rho*
			//v" << "," << "rho*w" << "," << "energy" << endl ;

			for (int i = 2; i < nx-2; ++i)
			{
				for (int j = 2; j < ny-2; ++j)
				{
					kullu_2D << U[i][j][2][0] << "," << U[i][j][2][1] << "," 
					<< U[i][j][2][2] << "," << U[i][j][2][3] << "," <<
					U[i][j][2][4] << endl ;
				}
			}
			/*storing all the conserved_variables at XY plane(2D plane) after
			some time step in the CSV format. so that we can plot them 
			in between the simulation or without getting it completed*/
		}
	} 
	// time marching ends here

	time(&end) ;
	double diff = difftime (end,start); //gives the time taken by the solver
	cout << "Time taken by the solver in secs = " << diff << endl ;
	return 0;
}

