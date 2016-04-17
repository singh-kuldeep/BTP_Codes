/*This class computes the net flux. Here net flux contains three fluxes, Euler
flux, Viscus flux and the numerical diffusion flux. No big calculation is
involved here. This class uses the there other classes and with help of them
it just sums up the three fluxes given by these classes*/

#include "math.h"
#include "iostream"
#include "eulerflux.h"
#include "viscusflux.h"
#include "diffusion_flux.h"

using namespace std ; 

class netflux
{
	public:
	float net_flux[5]; //Stores the net flux vector 
	
	netflux( vector<vector<vector<vector<float> > > > & U,
	vector<vector<vector<vector<float> > > > & x_face_area,
	vector<vector<vector<vector<float> > > > & y_face_area,
	vector<vector<vector<vector<float> > > > & z_face_area,
	vector<vector<vector<float> > >& V_cell, float dt, int i_l_m, int
	j_l_m, int k_l_m, int i_l, int j_l, int k_l, int i_r, int j_r, int k_r,
	int i_r_p, int j_r_p, int k_r_p, int viscus)
{
	/*U is conserved variable vector

	x_face_area is face area vectors which are in i direction
	y_face_area is face area vectors which are in j direction
	z_face_area is face area vectors which are in k direction

    (i_l, j_l, k_l) is the left cell location and (i_r, j_r, k_r) are the
	right cell location (i_l_m, j_l_m,k_l_m) is the cell location which is
	at left to the left cell (i_r_p,j_r_p,k_r_p) is the cell location which 
	is at right to the right cell*/

	eulerflux left ( U[i_l][j_l][k_l] );
	eulerflux right( U[i_r][j_r][k_r] );
	/*left and right Euler flux is calculated, using the eulerflux class,
	using the left and right cell conserved variables*/

    // Numerical Diffusion flux

    diffusion_flux diffusion( U[i_l_m][j_l_m][k_l_m], U[i_l][j_l][k_l], 
    	U[i_r][j_r][k_r],U[i_r_p][j_r_p][k_r_p], x_face_area[i_l][j_l][k_l],
    	x_face_area[i_r][j_r][k_r], x_face_area[i_r_p][j_r_p][k_r_p], 
    	V_cell[i_l_m][j_l_m][k_l_m], V_cell[i_l][j_l][k_l],
    	V_cell[i_r][j_r][k_r], V_cell[i_r_p][j_r_p][k_r_p],dt);

	float hx = x_face_area[i_r][j_r][k_r][0]/V_cell[i_r][j_r][k_r];
	float hy = x_face_area[i_r][j_r][k_r][1]/V_cell[i_r][j_r][k_r];
	float hz = x_face_area[i_r][j_r][k_r][2]/V_cell[i_r][j_r][k_r];

	if (i_r - i_l == 1) 
	{
		diffusion_flux diffusion( U[i_l_m][j_l_m][k_l_m], U[i_l][j_l][k_l],
		U[i_r][j_r][k_r], U[i_r_p][j_r_p][k_r_p], x_face_area[i_l][j_l][k_l],
		x_face_area[i_r][j_r][k_r], x_face_area[i_r_p][j_r_p][k_r_p],
		V_cell[i_l_m][j_l_m][k_l_m], V_cell[i_l][j_l][k_l],
		V_cell[i_r][j_r][k_r], V_cell[i_r_p][j_r_p][k_r_p], dt);

	    hx = x_face_area[i_r][j_r][k_r][0]/V_cell[i_r][j_r][k_r] ;
	    hy = x_face_area[i_r][j_r][k_r][1]/V_cell[i_r][j_r][k_r] ;
	    hz = x_face_area[i_r][j_r][k_r][2]/V_cell[i_r][j_r][k_r] ; 
	} 

	if (j_r - j_l == 1)
	{
		diffusion_flux diffusion( U[i_l_m][j_l_m][k_l_m], U[i_l][j_l][k_l],
		U[i_r][j_r][k_r], U[i_r_p][j_r_p][k_r_p], y_face_area[i_l][j_l][k_l],
		y_face_area[i_r][j_r][k_r], y_face_area[i_r_p][j_r_p][k_r_p], 
		V_cell[i_l_m][j_l_m][k_l_m], V_cell[i_l][j_l][k_l],
		V_cell[i_r][j_r][k_r],	V_cell[i_r_p][j_r_p][k_r_p], dt);

	    hx = y_face_area[i_r][j_r][k_r][0]/V_cell[i_r][j_r][k_r] ;
	    hy = y_face_area[i_r][j_r][k_r][1]/V_cell[i_r][j_r][k_r] ;
	    hz = y_face_area[i_r][j_r][k_r][2]/V_cell[i_r][j_r][k_r] ;
	}

	if (k_r - k_l == 1)
	{
		diffusion_flux diffusion(U[i_l_m][j_l_m][k_l_m], U[i_l][j_l][k_l],
		U[i_r][j_r][k_r], U[i_r_p][j_r_p][k_r_p], z_face_area[i_l][j_l][k_l],
		z_face_area[i_r][j_r][k_r], z_face_area[i_r_p][j_r_p][k_r_p],
		V_cell[i_l_m][j_l_m][k_l_m], V_cell[i_l][j_l][k_l],
		V_cell[i_r][j_r][k_r], V_cell[i_r_p][j_r_p][k_r_p], dt);

	    hx = z_face_area[i_r][j_r][k_r][0]/V_cell[i_r][j_r][k_r] ;
	    hy = z_face_area[i_r][j_r][k_r][1]/V_cell[i_r][j_r][k_r] ;
	    hz = z_face_area[i_r][j_r][k_r][2]/V_cell[i_r][j_r][k_r] ;
	}

	// Viscus flux
	viscusflux leftviscus( U, x_face_area, y_face_area, z_face_area, 
	V_cell, i_l, j_l, k_l, viscus) ; 
	viscusflux rightviscus(U, x_face_area, y_face_area, z_face_area,
	V_cell, i_r, j_r, k_r, viscus) ;

	// Averaged interface Euler flux 
	float E_int[5] ;  
	float F_int[5] ; 
	float G_int[5] ;

	float V_int = ( V_cell[i_l][j_l][k_l] +
		V_cell[i_r][j_r][k_r])/2 ;
	for (int i = 0; i < 5; ++i)
	{
	    E_int[i] = 0.5*(left.xeulerflux[i]+right.xeulerflux[i] -
	     	leftviscus.x_viscus_flux[i] - rightviscus.x_viscus_flux[i]) ;

	    F_int[i] = 0.5*(left.yeulerflux[i]+right.yeulerflux[i] -
	     	leftviscus.y_viscus_flux[i] - rightviscus.y_viscus_flux[i]) ;

		G_int[i] = 0.5*(left.zeulerflux[i]+right.zeulerflux[i] - 
			leftviscus.z_viscus_flux[i] - rightviscus.z_viscus_flux[i]) ;

		net_flux[i] = ( E_int[i]*hx + F_int[i]*hy + G_int[i]
			*hz)*V_int + 0.5*diffusion.diffusionfluxvector[i] ;
	}
	
};

	// ~netflux();
	
};