// it gives del(u,v,w,e)/ del(x,y,z) , In total 12 derivatives :)
// green's theorem is used to calculate the partial derivatives 

#include "math.h"
#include "iostream"

using namespace std ;

class partial_derivative
{
	public:
	float dubydx ;
	float dubydy ;
	float dubydz ;
	float dvbydx ;
	float dvbydy ;
	float dvbydz ;
	float dwbydx ;
	float dwbydy ;
	float dwbydz ;
	float debydx ;
	float debydy ;
	float debydz ;

	partial_derivative(vector<vector<vector<vector<float> > > > & U, 
		vector<vector<vector<vector<float> > > > & x_face_area,
		vector<vector<vector<vector<float> > > > & y_face_area,
		vector<vector<vector<vector<float> > > > & z_face_area,
		vector<vector<vector<float> > >& V_cell, int i, int j, int k)
	{
		float pressure = (gamma -1)*( U[i][j][k][4] - 0.5*(
			pow(U[i][j][k][1],2)+pow(U[i][j][k][2],2)+
			pow(U[i][j][k][3],2))/U[i][j][k][0] ) ;  
		
		// u
		dubydx = 0.5*( 
		   -(U[i][j][k][1]/U[i][j][k][0] + U[i-1][j][k][1]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][0]  +
			(U[i][j][k][1]/U[i][j][k][0] + U[i+1][j][k][1]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][0]+
		   -(U[i][j][k][1]/U[i][j][k][0] + U[i][j-1][k][1]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][0]  +
			(U[i][j][k][1]/U[i][j][k][0] + U[i][j+1][k][1]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][0]+
		   -(U[i][j][k][1]/U[i][j][k][0] + U[i][j][k-1][1]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][0]  +
			(U[i][j][k][1]/U[i][j][k][0] + U[i][j][k+1][1]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][0])
		/V_cell[i][j][k];

		dubydy = 0.5*( 
		   -(U[i][j][k][1]/U[i][j][k][0] + U[i-1][j][k][1]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][1]  +
			(U[i][j][k][1]/U[i][j][k][0] + U[i+1][j][k][1]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][1]+
		   -(U[i][j][k][1]/U[i][j][k][0] + U[i][j-1][k][1]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][1]  +
			(U[i][j][k][1]/U[i][j][k][0] + U[i][j+1][k][1]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][1]+
		   -(U[i][j][k][1]/U[i][j][k][0] + U[i][j][k-1][1]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][1]  +
			(U[i][j][k][1]/U[i][j][k][0] + U[i][j][k+1][1]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][1])
		/V_cell[i][j][k];

		dubydz = 0.5*( 
		   -(U[i][j][k][1]/U[i][j][k][0] + U[i-1][j][k][1]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][2]  +
			(U[i][j][k][1]/U[i][j][k][0] + U[i+1][j][k][1]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][2]+
		   -(U[i][j][k][1]/U[i][j][k][0] + U[i][j-1][k][1]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][2]  +
			(U[i][j][k][1]/U[i][j][k][0] + U[i][j+1][k][1]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][2]+
		   -(U[i][j][k][1]/U[i][j][k][0] + U[i][j][k-1][1]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][2]  +
			(U[i][j][k][1]/U[i][j][k][0] + U[i][j][k+1][1]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][2])
		/V_cell[i][j][k];

		// v
		dvbydx = 0.5*( 
		   -(U[i][j][k][2]/U[i][j][k][0] + U[i-1][j][k][2]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][0]  +
			(U[i][j][k][2]/U[i][j][k][0] + U[i+1][j][k][2]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][0]+
		   -(U[i][j][k][2]/U[i][j][k][0] + U[i][j-1][k][2]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][0]  +
			(U[i][j][k][2]/U[i][j][k][0] + U[i][j+1][k][2]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][0]+
		   -(U[i][j][k][2]/U[i][j][k][0] + U[i][j][k-1][2]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][0]  +
			(U[i][j][k][2]/U[i][j][k][0] + U[i][j][k+1][2]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][0])
		/V_cell[i][j][k];

		dvbydy = 0.5*( 
		   -(U[i][j][k][2]/U[i][j][k][0] + U[i-1][j][k][2]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][1]  +
			(U[i][j][k][2]/U[i][j][k][0] + U[i+1][j][k][2]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][1]+
		   -(U[i][j][k][2]/U[i][j][k][0] + U[i][j-1][k][2]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][1]  +
			(U[i][j][k][2]/U[i][j][k][0] + U[i][j+1][k][2]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][1]+
		   -(U[i][j][k][2]/U[i][j][k][0] + U[i][j][k-1][2]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][1]  +
			(U[i][j][k][2]/U[i][j][k][0] + U[i][j][k+1][2]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][1])
		/V_cell[i][j][k];

		dvbydz = 0.5*( 
		   -(U[i][j][k][2]/U[i][j][k][0] + U[i-1][j][k][2]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][2]  +
			(U[i][j][k][2]/U[i][j][k][0] + U[i+1][j][k][2]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][2]+
		   -(U[i][j][k][2]/U[i][j][k][0] + U[i][j-1][k][2]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][2]  +
			(U[i][j][k][2]/U[i][j][k][0] + U[i][j+1][k][2]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][2]+
		   -(U[i][j][k][2]/U[i][j][k][0] + U[i][j][k-1][2]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][2]  +
			(U[i][j][k][2]/U[i][j][k][0] + U[i][j][k+1][2]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][2])
		/V_cell[i][j][k];
		
		//w
		dwbydx = 0.5*( 
		   -(U[i][j][k][3]/U[i][j][k][0] + U[i-1][j][k][3]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][0]  +
			(U[i][j][k][3]/U[i][j][k][0] + U[i+1][j][k][3]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][0]+
		   -(U[i][j][k][3]/U[i][j][k][0] + U[i][j-1][k][3]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][0]  +
			(U[i][j][k][3]/U[i][j][k][0] + U[i][j+1][k][3]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][0]+
		   -(U[i][j][k][3]/U[i][j][k][0] + U[i][j][k-1][3]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][0]  +
			(U[i][j][k][3]/U[i][j][k][0] + U[i][j][k+1][3]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][0])
		/V_cell[i][j][k];

		dwbydy = 0.5*( 
		   -(U[i][j][k][3]/U[i][j][k][0] + U[i-1][j][k][3]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][1]  +
			(U[i][j][k][3]/U[i][j][k][0] + U[i+1][j][k][3]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][1]+
		   -(U[i][j][k][3]/U[i][j][k][0] + U[i][j-1][k][3]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][1]  +
			(U[i][j][k][3]/U[i][j][k][0] + U[i][j+1][k][3]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][1]+
		   -(U[i][j][k][3]/U[i][j][k][0] + U[i][j][k-1][3]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][1]  +
			(U[i][j][k][3]/U[i][j][k][0] + U[i][j][k+1][3]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][1])
		/V_cell[i][j][k];

		dwbydz = 0.5*( 
		   -(U[i][j][k][3]/U[i][j][k][0] + U[i-1][j][k][3]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][2]  +
			(U[i][j][k][3]/U[i][j][k][0] + U[i+1][j][k][3]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][2]+
		   -(U[i][j][k][3]/U[i][j][k][0] + U[i][j-1][k][3]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][2]  +
			(U[i][j][k][3]/U[i][j][k][0] + U[i][j+1][k][3]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][2]+
		   -(U[i][j][k][3]/U[i][j][k][0] + U[i][j][k-1][3]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][2]  +
			(U[i][j][k][3]/U[i][j][k][0] + U[i][j][k+1][3]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][2])
		/V_cell[i][j][k];

		//e 
		debydx = 0.5*( 
		   -(U[i][j][k][4]/U[i][j][k][0] + U[i-1][j][k][4]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][0]  +
			(U[i][j][k][4]/U[i][j][k][0] + U[i+1][j][k][4]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][0]+
		   -(U[i][j][k][4]/U[i][j][k][0] + U[i][j-1][k][4]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][0]  +
			(U[i][j][k][4]/U[i][j][k][0] + U[i][j+1][k][4]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][0]+
		   -(U[i][j][k][4]/U[i][j][k][0] + U[i][j][k-1][4]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][0]  +
			(U[i][j][k][4]/U[i][j][k][0] + U[i][j][k+1][4]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][0])
		/V_cell[i][j][k];

		debydy = 0.5*( 
		   -(U[i][j][k][4]/U[i][j][k][0] + U[i-1][j][k][4]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][1]  +
			(U[i][j][k][4]/U[i][j][k][0] + U[i+1][j][k][4]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][1]+
		   -(U[i][j][k][4]/U[i][j][k][0] + U[i][j-1][k][4]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][1]  +
			(U[i][j][k][4]/U[i][j][k][0] + U[i][j+1][k][4]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][1]+
		   -(U[i][j][k][4]/U[i][j][k][0] + U[i][j][k-1][4]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][1]  +
			(U[i][j][k][4]/U[i][j][k][0] + U[i][j][k+1][4]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][1])
		/V_cell[i][j][k];

		debydz = 0.5*( 
		   -(U[i][j][k][4]/U[i][j][k][0] + U[i-1][j][k][4]/U[i-1][j][k][0])*
		   x_face_area[i][j][k][2]  +
			(U[i][j][k][4]/U[i][j][k][0] + U[i+1][j][k][4]/U[i+1][j][k][0])*
			x_face_area[i+1][j][k][2]+
		   -(U[i][j][k][4]/U[i][j][k][0] + U[i][j-1][k][4]/U[i][j-1][k][0])*
		   y_face_area[i][j][k][2]  +
			(U[i][j][k][4]/U[i][j][k][0] + U[i][j+1][k][4]/U[i+1][j][k][0])*
			y_face_area[i][j+1][k][2]+
		   -(U[i][j][k][4]/U[i][j][k][0] + U[i][j][k-1][4]/U[i][j][k-1][0])*
		   z_face_area[i][j][k][2]  +
			(U[i][j][k][4]/U[i][j][k][0] + U[i][j][k+1][4]/U[i][j][k+1][0])*
			z_face_area[i][j][k+1][2])
		/V_cell[i][j][k];
		
	 };
	 // ~eulerflux();	
};