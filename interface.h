// checked
#include "math.h"
#include "iostream"

#define GAMMA 1.4

using namespace std ;

class interface
{
	public:
	float density_interface ;
	float u_interface ;
	float v_interface ;
	float w_interface ;
	float enthalpy_interface ;

	float U_jump_interface[5] ;
	float eigenvalue[5] ;


	float eigenvectormatrix[5][5] ;
	float eigenvectormatrixinvers[5][5] ; 

	float alpha_interface[5] ;
	// eigenvectormatrixinvers[5][5]*U_jump_interface

	float mu_interface[5] ;
	// this is some parameter = delta t * eigenvalue

	float Z_interface[5] ; // this is dame as mu_interface

	float pshi_interface[5] ;
	float g_interface[5];


	interface(vector<float>& U_l,vector<float>& U_r, vector<float>& 
		A_interface, float V_l,	float V_r, float dt) 
	{
		float V_interface = 0.5*(V_l + V_r) ;

		float pressure_left = (GAMMA -1)*( U_l[4] - 0.5*(pow(U_l[1],2)+pow(
			U_l[2],2)+ pow(U_l[3],2))/U_l[0] ) ;

		float pressure_right = (GAMMA -1)*( U_r[4] - 0.5*(pow(U_r[1],2)+pow(
			U_r[2],2)+ pow(U_r[3],2))/U_r[0] ) ;
		
		float enthalpyleft = (U_l[4] + pressure_left)/U_l[0];
		float enthalpyright = (U_r[4] + pressure_right)/U_r[0];

		density_interface = sqrt(U_l[0]*U_r[0]) ; 

		u_interface = (U_l[1] + U_r[1]*sqrt(U_l[0]/U_r[0]))/(U_l[0] +sqrt(
			U_r[0]*U_l[0])) ;

		v_interface = (U_l[2] + U_r[2]*sqrt(U_l[0]/U_r[0]))/(U_l[0] +sqrt(
			U_r[0]*U_l[0])) ;

		w_interface = (U_l[3] + U_r[3]*sqrt(U_l[0]/U_r[0]))/(U_l[0]+sqrt(
			U_r[0]*U_l[0])) ;

		enthalpy_interface = (enthalpyleft + enthalpyright*sqrt(
			U_r[0]/U_l[0]))/(1+sqrt(U_r[0]/U_l[0])) ; 

		// eigenvalues 
		float hx = A_interface[0] / V_interface ;
		float hy = A_interface[1] / V_interface ;
		float hz = A_interface[2] / V_interface ;

		float hn = sqrt(pow(A_interface[0],2) + pow( A_interface[1],2) +
			pow(A_interface[2],2)) / V_interface ;

		float Ucont = u_interface*hx + v_interface*hy + w_interface*hz ;

		float a_interface = sqrt((GAMMA -1)*(enthalpy_interface - 0.5*(pow(
			u_interface,2) + pow(v_interface,2) + pow(w_interface,2)))) ;
		// a_interface is sound velocity at interface

		eigenvalue[0] = Ucont - a_interface*hn ;
		eigenvalue[1] = Ucont ;
		eigenvalue[2] = Ucont ;
		eigenvalue[3] = Ucont ;
		eigenvalue[4] = Ucont + a_interface*hn ;

		// Jump vector at interface 
		U_jump_interface[0] = V_interface*(U_r[0] - U_l[0]) ;
		U_jump_interface[1] = V_interface*(U_r[1] - U_l[1]) ;
		U_jump_interface[2] = V_interface*(U_r[2] - U_l[2]) ;
		U_jump_interface[3] = V_interface*(U_r[3] - U_l[3]) ;
		U_jump_interface[4] = V_interface*(U_r[4] - U_l[4]) ;

		// Defying the eigenvector matrix "R"
		float hdesx = hx/hn ; 
		float hdesy = hy/hn ; 
		float hdesz = hz/hn ; 
		float Phi = u_interface*hdesx + v_interface*hdesy + w_interface*hdesz;
		float q = sqrt(pow(u_interface,2)+pow(v_interface,2)+pow(
			w_interface,2)) ;

		eigenvectormatrix[0][0] = 1 ; 
		eigenvectormatrix[0][1] = 1 ;
		eigenvectormatrix[0][2] = 0 ;
		eigenvectormatrix[0][3] = 0 ; 
		eigenvectormatrix[0][4] = 1 ;

		eigenvectormatrix[1][0] = u_interface - hdesx * a_interface ; 
		eigenvectormatrix[1][1] = u_interface ;
		eigenvectormatrix[1][2] = hdesy ;
		eigenvectormatrix[1][3] = hdesz ; 
		eigenvectormatrix[1][4] = u_interface + hdesx*a_interface ;

		eigenvectormatrix[2][0] = v_interface - hdesy * a_interface ; 
		eigenvectormatrix[2][1] = v_interface ;
		eigenvectormatrix[2][2] = hdesz ;
		eigenvectormatrix[2][3] = hdesx ; 
		eigenvectormatrix[2][4] = v_interface + hdesy*a_interface ;

		eigenvectormatrix[3][0] = w_interface - hdesz * a_interface ; 
		eigenvectormatrix[3][1] = w_interface ;
		eigenvectormatrix[3][2] = hdesx ;
		eigenvectormatrix[3][3] = hdesy ; 
		eigenvectormatrix[3][4] = w_interface + hdesz*a_interface ;

		eigenvectormatrix[4][0] = enthalpy_interface - (hdesx*u_interface + 
			hdesy*v_interface + hdesz*w_interface) * a_interface ; 
		eigenvectormatrix[4][1] = 0.5 * pow(q,2);
		eigenvectormatrix[4][2] = hdesx*w_interface + hdesz*v_interface + 
			hdesy*u_interface;
		eigenvectormatrix[4][3] = hdesy*w_interface + hdesx*v_interface + 
			hdesz*u_interface;
		eigenvectormatrix[4][4] = enthalpy_interface + (hdesx*u_interface +
			hdesy*v_interface + hdesz*w_interface) * a_interface ;

		// eigenvectors matrix defying
		eigenvectormatrixinvers[0][0] = 0.5*(0.5*(pow(q,2))*((GAMMA-1)/pow(
			a_interface,2))+(Phi/a_interface));
		eigenvectormatrixinvers[0][1] = -0.5*(u_interface*((GAMMA-1)/pow(
			a_interface,2))+(hdesx/a_interface));
		eigenvectormatrixinvers[0][2] = -0.5*(v_interface*((GAMMA-1)/pow(
			a_interface,2))+(hdesy/a_interface));
		eigenvectormatrixinvers[0][3] = -0.5*(w_interface*((GAMMA-1)/pow(
			a_interface,2))+(hdesz/a_interface));
		eigenvectormatrixinvers[0][4] = 0.5*(GAMMA -1)/ pow(a_interface,2) ;

		eigenvectormatrixinvers[1][0] = 1 - 0.5*pow(q,2)*((GAMMA-1)/pow(
			a_interface,2)) ; 
		eigenvectormatrixinvers[1][1] = u_interface*(GAMMA-1)/pow(
			a_interface,2) ; 
		eigenvectormatrixinvers[1][2] = v_interface*(GAMMA-1)/pow(
			a_interface,2) ;   
		eigenvectormatrixinvers[1][3] = w_interface*(GAMMA-1)/pow(
			a_interface,2) ; 
		eigenvectormatrixinvers[1][4] = -(GAMMA-1)/pow(a_interface,2) ;

		eigenvectormatrixinvers[2][0] = -(hdesy*u_interface + hdesz*
			v_interface + hdesx*w_interface);
		eigenvectormatrixinvers[2][1] = hdesy ; 
		eigenvectormatrixinvers[2][2] = hdesz ;
		eigenvectormatrixinvers[2][3] = hdesx ; 
		eigenvectormatrixinvers[2][4] = 0 ;

		eigenvectormatrixinvers[3][0] = -(hdesz*u_interface + hdesx*
			v_interface + hdesy*w_interface);
		eigenvectormatrixinvers[3][1] = hdesz ; 
		eigenvectormatrixinvers[3][2] = hdesx ;
		eigenvectormatrixinvers[3][3] = hdesy ; 
		eigenvectormatrixinvers[3][4] = 0 ;

		eigenvectormatrixinvers[4][0] = 0.5*(0.5*(pow(q,2))*((GAMMA-1)/pow(
			a_interface,2)) - (Phi/a_interface));
		eigenvectormatrixinvers[4][1] = 0.5*( -u_interface*((GAMMA-1)/pow(
			a_interface,2))+ (hdesx/a_interface));
		eigenvectormatrixinvers[4][2] = 0.5*( -v_interface*((GAMMA-1)/pow(
			a_interface,2))+ (hdesy/a_interface));
		eigenvectormatrixinvers[4][3] = 0.5*( -w_interface*((GAMMA-1)/pow(
			a_interface,2))+ (hdesz/a_interface));
		eigenvectormatrixinvers[4][4] = 0.5*(GAMMA -1)/ pow(a_interface,2) ;

		//alpha_interface defying 
 		for (int i = 0; i < 5; ++i)
		{
			alpha_interface[i] = 0.0 ;
			for (int l = 0; l < 5; ++l)
			{
				alpha_interface[i] += eigenvectormatrixinvers[i][l]*
				U_jump_interface[l] ;
			}
			// cout <<   alpha_interface[i] << endl ;
		}

		// Defying the mu_interface or Z_interface
		for (int i = 0; i < 5; ++i)
		{
			mu_interface[i] = dt * eigenvalue[i] ;
			Z_interface[i] = dt * eigenvalue[i] ;
		}
			
// // pshi_interface defying (for invisid flow)
// 		for (int i = 0; i < 5; ++i)
// 		{
// 			pshi_interface[i] = pow(Z_interface[i],2) + 0.25 ;
// 		}
// // pshi_interface defying (for viscus flow)
		float deltaf = 0.2 ; 
		// this is given constant value (it is between 0.1 to 0.5)
		
		for (int i = 0; i < 5; ++i)
		{
			if(fabs(Z_interface[i]) >= deltaf ){
				pshi_interface[i] = fabs(Z_interface[i]);
			}
			else
			{
				pshi_interface[i] = 0.5*(pow(Z_interface[i],2)+pow(deltaf,2)) 
				/ deltaf ;
			}
		}

		// g_interface defying
		for (int i = 0; i < 5; ++i)
		{
			g_interface[i] = 0.5*(pshi_interface[i] - 
				pow(Z_interface[i],2))*alpha_interface[i];
		}
	 };
	// ~interface();
};