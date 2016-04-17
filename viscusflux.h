/*This calss calculates the viscus flux at the interface of the two cells*/

#include "math.h"
#include "iostream"
#include "partial_derivative.h"

#define GAMMA 1.4
#define R 287.14
#define CP 717.5
#define b 1.458e-6 
#define S 110.4 
#define Prd 0.72

// using namespace std ;

class viscusflux
{
	public:
	float x_viscus_flux[5] ;
	float y_viscus_flux[5] ;
	float z_viscus_flux[5] ;

	/*x_vidcus_flux is viscus flux at the cell interface which is in the i 
	direction of the structured grid*/
	
	viscusflux(vector<vector<vector<vector<float> > > > & U, 
		vector<vector<vector<vector<float> > > > & x_face_area,
		vector<vector<vector<vector<float> > > > & y_face_area,
		vector<vector<vector<vector<float> > > > & z_face_area,
		vector<vector<vector<float> > >& cell_volume, int i, int j, int k, int viscus)
	{
		partial_derivative del(U, x_face_area, y_face_area, z_face_area, 
			cell_volume, i, j, k) ;

		float pressure = (GAMMA -1)*(U[i][j][k][4]-0.5*(pow(U[i][j][k][1],2)+
			pow(U[i][j][k][2],2)+pow(U[i][j][k][3],2))/U[i][j][k][0]) ;

		float temperature = pressure / (R*U[i][j][k][0]) ;

		float mu_M = b*sqrt(temperature)/(1+S/temperature) ;
	if (viscus == 1)
		{
		// viscus flux
		x_viscus_flux[0] = 0 ; 
		x_viscus_flux[1] = 2*mu_M*del.dubydx - 2*mu_M*(del.dubydx + del.dvbydy
			+ del.dwbydz)/3 ; // tau_xx
		x_viscus_flux[2] = mu_M*(del.dubydy + del.dvbydx) ; // tau_xy
		x_viscus_flux[3] = mu_M*(del.dubydz + del.dwbydx) ; // tau_xz
		x_viscus_flux[4] = (x_viscus_flux[1]*U[i][j][k][1] + x_viscus_flux[2]*
			U[i][j][k][2] + x_viscus_flux[3]*U[i][j][k][3])/U[i][j][k][0] + 
			(GAMMA*mu_M/Prd)*del.debydx ;

		y_viscus_flux[0] = 0 ; 
		y_viscus_flux[1] = x_viscus_flux[2] ; // tau_yx = tau_xy
		y_viscus_flux[2] = (2*mu_M*del.dvbydy-2*mu_M*(del.dubydx + del.dvbydy 
			+ del.dwbydz)/3) ; // tau_yy
		y_viscus_flux[3] = (mu_M*(del.dvbydz + del.dwbydy)) ; // tau_yz
		y_viscus_flux[4] =(y_viscus_flux[1]*U[i][j][k][1] + y_viscus_flux[2]*
			U[i][j][k][2] + y_viscus_flux[3]*U[i][j][k][3])/U[i][j][k][0] + 
			(GAMMA*mu_M/Prd)*del.debydy;

		z_viscus_flux[0] = 0 ; 
		z_viscus_flux[1] = x_viscus_flux[3] ; // tau_zx = tau_xz
		z_viscus_flux[2] = y_viscus_flux[3] ; // tau_zy
		z_viscus_flux[2] = 2*mu_M*del.dwbydz - 2*mu_M*(del.dubydx + del.dvbydy
			+ del.dwbydz)/3 ; // tau_zz
		z_viscus_flux[4] = (z_viscus_flux[1]*U[i][j][k][1] + z_viscus_flux[2]*
			U[i][j][k][2] + z_viscus_flux[3]*U[i][j][k][3])/U[i][j][k][0] +
			(GAMMA*mu_M/Prd)*del.debydz ;		
		}	

	if(viscus == 0)
		{
		// INVISID CASE
		x_viscus_flux[0] = 0 ;
		x_viscus_flux[1] = 0 ;
		x_viscus_flux[2] = 0 ;
		x_viscus_flux[3] = 0 ;
		x_viscus_flux[4] = 0 ;

		y_viscus_flux[0] = 0 ;
		y_viscus_flux[1] = 0 ;
		y_viscus_flux[2] = 0 ;
		y_viscus_flux[3] = 0 ;
		y_viscus_flux[4] = 0 ;

		z_viscus_flux[0] = 0 ;
		z_viscus_flux[1] = 0 ;
		z_viscus_flux[2] = 0 ;
		z_viscus_flux[3] = 0 ;
		z_viscus_flux[4] = 0 ;
		}
	
	};
	 // ~viscusflux();	
};